/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019,2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements common routines for PME tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "pmetestcommon.h"

#include <cstring>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme_gather.h"
#include "gromacs/ewald/pme_gpu_internal.h"
#include "gromacs/ewald/pme_gpu_staging.h"
#include "gromacs/ewald/pme_grid.h"
#include "gromacs/ewald/pme_internal.h"
#include "gromacs/ewald/pme_redistribute.h"
#include "gromacs/ewald/pme_solve.h"
#include "gromacs/ewald/pme_spread.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

bool pmeSupportsInputForMode(const gmx_hw_info_t& hwinfo, const t_inputrec* inputRec, CodePath mode)
{
    bool       implemented;
    gmx_mtop_t mtop;
    switch (mode)
    {
        case CodePath::CPU: implemented = true; break;

        case CodePath::GPU:
            implemented = (pme_gpu_supports_build(nullptr) && pme_gpu_supports_hardware(hwinfo, nullptr)
                           && pme_gpu_supports_input(*inputRec, mtop, nullptr));
            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
    return implemented;
}

uint64_t getSplineModuliDoublePrecisionUlps(int splineOrder)
{
    /* Arbitrary ulp tolerance for sine/cosine implementation. It's
     * hard to know what to pick without testing lots of
     * implementations. */
    const uint64_t sineUlps = 10;
    return 4 * (splineOrder - 2) + 2 * sineUlps * splineOrder;
}

//! PME initialization
PmeSafePointer pmeInitWrapper(const t_inputrec*        inputRec,
                              const CodePath           mode,
                              const gmx_device_info_t* gpuInfo,
                              const PmeGpuProgram*     pmeGpuProgram,
                              const Matrix3x3&         box,
                              const real               ewaldCoeff_q,
                              const real               ewaldCoeff_lj)
{
    const MDLogger dummyLogger;
    const auto     runMode       = (mode == CodePath::CPU) ? PmeRunMode::CPU : PmeRunMode::Mixed;
    t_commrec      dummyCommrec  = { 0 };
    NumPmeDomains  numPmeDomains = { 1, 1 };
    gmx_pme_t*     pmeDataRaw =
            gmx_pme_init(&dummyCommrec, numPmeDomains, inputRec, false, false, true, ewaldCoeff_q,
                         ewaldCoeff_lj, 1, runMode, nullptr, gpuInfo, pmeGpuProgram, dummyLogger);
    PmeSafePointer pme(pmeDataRaw); // taking ownership

    // TODO get rid of this with proper matrix type
    matrix boxTemp;
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            boxTemp[i][j] = box[i * DIM + j];
        }
    }
    const char* boxError = check_box(PbcType::Unset, boxTemp);
    GMX_RELEASE_ASSERT(boxError == nullptr, boxError);

    switch (mode)
    {
        case CodePath::CPU: invertBoxMatrix(boxTemp, pme->recipbox); break;

        case CodePath::GPU:
            pme_gpu_set_testing(pme->gpu, true);
            pme_gpu_update_input_box(pme->gpu, boxTemp);
            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }

    return pme;
}

//! Simple PME initialization based on input, no atom data
PmeSafePointer pmeInitEmpty(const t_inputrec*        inputRec,
                            const CodePath           mode,
                            const gmx_device_info_t* gpuInfo,
                            const PmeGpuProgram*     pmeGpuProgram,
                            const Matrix3x3&         box,
                            real                     ewaldCoeff_q,
                            real                     ewaldCoeff_lj)
{
    return pmeInitWrapper(inputRec, mode, gpuInfo, pmeGpuProgram, box, ewaldCoeff_q, ewaldCoeff_lj);
    // hiding the fact that PME actually needs to know the number of atoms in advance
}

//! Make a GPU state-propagator manager
std::unique_ptr<StatePropagatorDataGpu> makeStatePropagatorDataGpu(const gmx_pme_t& pme)
{
    // TODO: Pin the host buffer and use async memory copies
    // TODO: Special constructor for PME-only rank / PME-tests is used here. There should be a mechanism to
    //       restrict one from using other constructor here.
    return std::make_unique<StatePropagatorDataGpu>(
            pme_gpu_get_device_stream(&pme), pme_gpu_get_device_context(&pme),
            GpuApiCallBehavior::Sync, pme_gpu_get_padding_size(&pme), nullptr);
}

//! PME initialization with atom data
void pmeInitAtoms(gmx_pme_t*               pme,
                  StatePropagatorDataGpu*  stateGpu,
                  const CodePath           mode,
                  const CoordinatesVector& coordinates,
                  const ChargesVector&     charges)
{
    const index atomCount = coordinates.size();
    GMX_RELEASE_ASSERT(atomCount == charges.ssize(), "Mismatch in atom data");
    PmeAtomComm* atc = nullptr;

    switch (mode)
    {
        case CodePath::CPU:
            atc              = &(pme->atc[0]);
            atc->x           = coordinates;
            atc->coefficient = charges;
            gmx_pme_reinit_atoms(pme, atomCount, charges.data());
            /* With decomposition there would be more boilerplate atc code here, e.g. do_redist_pos_coeffs */
            break;

        case CodePath::GPU:
            // TODO: Avoid use of atc in the GPU code path
            atc = &(pme->atc[0]);
            // We need to set atc->n for passing the size in the tests
            atc->setNumAtoms(atomCount);
            gmx_pme_reinit_atoms(pme, atomCount, charges.data());

            stateGpu->reinit(atomCount, atomCount);
            stateGpu->copyCoordinatesToGpu(arrayRefFromArray(coordinates.data(), coordinates.size()),
                                           gmx::AtomLocality::All);
            pme_gpu_set_kernelparam_coordinates(pme->gpu, stateGpu->getCoordinates());

            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Getting local PME real grid pointer for test I/O
static real* pmeGetRealGridInternal(const gmx_pme_t* pme)
{
    const size_t gridIndex = 0;
    return pme->fftgrid[gridIndex];
}

//! Getting local PME real grid dimensions
static void pmeGetRealGridSizesInternal(const gmx_pme_t* pme,
                                        CodePath         mode,
                                        IVec& gridSize,       //NOLINT(google-runtime-references)
                                        IVec& paddedGridSize) //NOLINT(google-runtime-references)
{
    const size_t gridIndex = 0;
    IVec         gridOffsetUnused;
    switch (mode)
    {
        case CodePath::CPU:
            gmx_parallel_3dfft_real_limits(pme->pfft_setup[gridIndex], gridSize, gridOffsetUnused,
                                           paddedGridSize);
            break;

        case CodePath::GPU:
            pme_gpu_get_real_grid_sizes(pme->gpu, &gridSize, &paddedGridSize);
            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Getting local PME complex grid pointer for test I/O
static t_complex* pmeGetComplexGridInternal(const gmx_pme_t* pme)
{
    const size_t gridIndex = 0;
    return pme->cfftgrid[gridIndex];
}

//! Getting local PME complex grid dimensions
static void pmeGetComplexGridSizesInternal(const gmx_pme_t* pme,
                                           IVec& gridSize,       //NOLINT(google-runtime-references)
                                           IVec& paddedGridSize) //NOLINT(google-runtime-references)
{
    const size_t gridIndex = 0;
    IVec         gridOffsetUnused, complexOrderUnused;
    gmx_parallel_3dfft_complex_limits(pme->pfft_setup[gridIndex], complexOrderUnused, gridSize,
                                      gridOffsetUnused, paddedGridSize); // TODO: what about YZX ordering?
}

//! Getting the PME grid memory buffer and its sizes - template definition
template<typename ValueType>
static void pmeGetGridAndSizesInternal(const gmx_pme_t* /*unused*/,
                                       CodePath /*unused*/,
                                       ValueType*& /*unused*/, //NOLINT(google-runtime-references)
                                       IVec& /*unused*/,       //NOLINT(google-runtime-references)
                                       IVec& /*unused*/)       //NOLINT(google-runtime-references)
{
    GMX_THROW(InternalError("Deleted function call"));
    // explicitly deleting general template does not compile in clang/icc, see https://llvm.org/bugs/show_bug.cgi?id=17537
}

//! Getting the PME real grid memory buffer and its sizes
template<>
void pmeGetGridAndSizesInternal<real>(const gmx_pme_t* pme, CodePath mode, real*& grid, IVec& gridSize, IVec& paddedGridSize)
{
    grid = pmeGetRealGridInternal(pme);
    pmeGetRealGridSizesInternal(pme, mode, gridSize, paddedGridSize);
}

//! Getting the PME complex grid memory buffer and its sizes
template<>
void pmeGetGridAndSizesInternal<t_complex>(const gmx_pme_t* pme,
                                           CodePath /*unused*/,
                                           t_complex*& grid,
                                           IVec&       gridSize,
                                           IVec&       paddedGridSize)
{
    grid = pmeGetComplexGridInternal(pme);
    pmeGetComplexGridSizesInternal(pme, gridSize, paddedGridSize);
}

//! PME spline calculation and charge spreading
void pmePerformSplineAndSpread(gmx_pme_t* pme,
                               CodePath   mode, // TODO const qualifiers elsewhere
                               bool       computeSplines,
                               bool       spreadCharges)
{
    GMX_RELEASE_ASSERT(pme != nullptr, "PME data is not initialized");
    PmeAtomComm* atc                          = &(pme->atc[0]);
    const size_t gridIndex                    = 0;
    const bool   computeSplinesForZeroCharges = true;
    real*        fftgrid                      = spreadCharges ? pme->fftgrid[gridIndex] : nullptr;
    real*        pmegrid                      = pme->pmegrid[gridIndex].grid.grid;

    switch (mode)
    {
        case CodePath::CPU:
            spread_on_grid(pme, atc, &pme->pmegrid[gridIndex], computeSplines, spreadCharges,
                           fftgrid, computeSplinesForZeroCharges, gridIndex);
            if (spreadCharges && !pme->bUseThreads)
            {
                wrap_periodic_pmegrid(pme, pmegrid);
                copy_pmegrid_to_fftgrid(pme, pmegrid, fftgrid, gridIndex);
            }
            break;

        case CodePath::GPU:
        {
            // no synchronization needed as x is transferred in the PME stream
            GpuEventSynchronizer* xReadyOnDevice = nullptr;
            pme_gpu_spread(pme->gpu, xReadyOnDevice, gridIndex, fftgrid, computeSplines, spreadCharges);
        }
        break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Getting the internal spline data buffer pointer
static real* pmeGetSplineDataInternal(const gmx_pme_t* pme, PmeSplineDataType type, int dimIndex)
{
    GMX_ASSERT((0 <= dimIndex) && (dimIndex < DIM), "Invalid dimension index");
    const PmeAtomComm* atc          = &(pme->atc[0]);
    const size_t       threadIndex  = 0;
    real*              splineBuffer = nullptr;
    switch (type)
    {
        case PmeSplineDataType::Values:
            splineBuffer = atc->spline[threadIndex].theta.coefficients[dimIndex];
            break;

        case PmeSplineDataType::Derivatives:
            splineBuffer = atc->spline[threadIndex].dtheta.coefficients[dimIndex];
            break;

        default: GMX_THROW(InternalError("Unknown spline data type"));
    }
    return splineBuffer;
}

//! PME solving
void pmePerformSolve(const gmx_pme_t*  pme,
                     CodePath          mode,
                     PmeSolveAlgorithm method,
                     real              cellVolume,
                     GridOrdering      gridOrdering,
                     bool              computeEnergyAndVirial)
{
    t_complex*   h_grid              = pmeGetComplexGridInternal(pme);
    const bool   useLorentzBerthelot = false;
    const size_t threadIndex         = 0;
    switch (mode)
    {
        case CodePath::CPU:
            if (gridOrdering != GridOrdering::YZX)
            {
                GMX_THROW(InternalError("Test not implemented for this mode"));
            }
            switch (method)
            {
                case PmeSolveAlgorithm::Coulomb:
                    solve_pme_yzx(pme, h_grid, cellVolume, computeEnergyAndVirial, pme->nthread, threadIndex);
                    break;

                case PmeSolveAlgorithm::LennardJones:
                    solve_pme_lj_yzx(pme, &h_grid, useLorentzBerthelot, cellVolume,
                                     computeEnergyAndVirial, pme->nthread, threadIndex);
                    break;

                default: GMX_THROW(InternalError("Test not implemented for this mode"));
            }
            break;

        case CodePath::GPU:
            switch (method)
            {
                case PmeSolveAlgorithm::Coulomb:
                    pme_gpu_solve(pme->gpu, h_grid, gridOrdering, computeEnergyAndVirial);
                    break;

                default: GMX_THROW(InternalError("Test not implemented for this mode"));
            }
            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! PME force gathering
void pmePerformGather(gmx_pme_t* pme, CodePath mode, PmeForceOutputHandling inputTreatment, ForcesVector& forces)
{
    PmeAtomComm* atc       = &(pme->atc[0]);
    const index  atomCount = atc->numAtoms();
    GMX_RELEASE_ASSERT(forces.ssize() == atomCount, "Invalid force buffer size");
    const bool forceReductionWithInput = (inputTreatment == PmeForceOutputHandling::ReduceWithInput);
    const real   scale                 = 1.0;
    const size_t threadIndex           = 0;
    const size_t gridIndex             = 0;
    real*        pmegrid               = pme->pmegrid[gridIndex].grid.grid;
    real*        fftgrid               = pme->fftgrid[gridIndex];

    switch (mode)
    {
        case CodePath::CPU:
            atc->f = forces;
            if (atc->nthread == 1)
            {
                // something which is normally done in serial spline computation (make_thread_local_ind())
                atc->spline[threadIndex].n = atomCount;
            }
            copy_fftgrid_to_pmegrid(pme, fftgrid, pmegrid, gridIndex, pme->nthread, threadIndex);
            unwrap_periodic_pmegrid(pme, pmegrid);
            gather_f_bsplines(pme, pmegrid, !forceReductionWithInput, atc, &atc->spline[threadIndex], scale);
            break;

        case CodePath::GPU:
        {
            // Variable initialization needs a non-switch scope
            PmeOutput output = pme_gpu_getOutput(*pme, GMX_PME_CALC_F);
            GMX_ASSERT(forces.size() == output.forces_.size(),
                       "Size of force buffers did not match");
            if (forceReductionWithInput)
            {
                std::copy(std::begin(forces), std::end(forces), std::begin(output.forces_));
            }
            pme_gpu_gather(pme->gpu, inputTreatment, reinterpret_cast<float*>(fftgrid));
            std::copy(std::begin(output.forces_), std::end(output.forces_), std::begin(forces));
        }
        break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! PME test finalization before fetching the outputs
void pmeFinalizeTest(const gmx_pme_t* pme, CodePath mode)
{
    switch (mode)
    {
        case CodePath::CPU: break;

        case CodePath::GPU: pme_gpu_synchronize(pme->gpu); break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Setting atom spline values/derivatives to be used in spread/gather
void pmeSetSplineData(const gmx_pme_t*             pme,
                      CodePath                     mode,
                      const SplineParamsDimVector& splineValues,
                      PmeSplineDataType            type,
                      int                          dimIndex)
{
    const PmeAtomComm* atc       = &(pme->atc[0]);
    const index        atomCount = atc->numAtoms();
    const index        pmeOrder  = pme->pme_order;
    const index        dimSize   = pmeOrder * atomCount;
    GMX_RELEASE_ASSERT(dimSize == splineValues.ssize(), "Mismatch in spline data");
    real* splineBuffer = pmeGetSplineDataInternal(pme, type, dimIndex);

    switch (mode)
    {
        case CodePath::CPU:
            std::copy(splineValues.begin(), splineValues.end(), splineBuffer);
            break;

        case CodePath::GPU:
            std::copy(splineValues.begin(), splineValues.end(), splineBuffer);
            pme_gpu_transform_spline_atom_data(pme->gpu, atc, type, dimIndex, PmeLayoutTransform::HostToGpu);
            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Setting gridline indices to be used in spread/gather
void pmeSetGridLineIndices(gmx_pme_t* pme, CodePath mode, const GridLineIndicesVector& gridLineIndices)
{
    PmeAtomComm* atc       = &(pme->atc[0]);
    const index  atomCount = atc->numAtoms();
    GMX_RELEASE_ASSERT(atomCount == gridLineIndices.ssize(), "Mismatch in gridline indices size");

    IVec paddedGridSizeUnused, gridSize(0, 0, 0);
    pmeGetRealGridSizesInternal(pme, mode, gridSize, paddedGridSizeUnused);

    for (const auto& index : gridLineIndices)
    {
        for (int i = 0; i < DIM; i++)
        {
            GMX_RELEASE_ASSERT((0 <= index[i]) && (index[i] < gridSize[i]),
                               "Invalid gridline index");
        }
    }

    switch (mode)
    {
        case CodePath::GPU:
            memcpy(pme_gpu_staging(pme->gpu).h_gridlineIndices, gridLineIndices.data(),
                   atomCount * sizeof(gridLineIndices[0]));
            break;

        case CodePath::CPU:
            atc->idx.resize(gridLineIndices.size());
            std::copy(gridLineIndices.begin(), gridLineIndices.end(), atc->idx.begin());
            break;
        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Getting plain index into the complex 3d grid
inline size_t pmeGetGridPlainIndexInternal(const IVec& index, const IVec& paddedGridSize, GridOrdering gridOrdering)
{
    size_t result;
    switch (gridOrdering)
    {
        case GridOrdering::YZX:
            result = (index[YY] * paddedGridSize[ZZ] + index[ZZ]) * paddedGridSize[XX] + index[XX];
            break;

        case GridOrdering::XYZ:
            result = (index[XX] * paddedGridSize[YY] + index[YY]) * paddedGridSize[ZZ] + index[ZZ];
            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
    return result;
}

//! Setting real or complex grid
template<typename ValueType>
static void pmeSetGridInternal(const gmx_pme_t*                        pme,
                               CodePath                                mode,
                               GridOrdering                            gridOrdering,
                               const SparseGridValuesInput<ValueType>& gridValues)
{
    IVec       gridSize(0, 0, 0), paddedGridSize(0, 0, 0);
    ValueType* grid;
    pmeGetGridAndSizesInternal<ValueType>(pme, mode, grid, gridSize, paddedGridSize);

    switch (mode)
    {
        case CodePath::GPU: // intentional absence of break, the grid will be copied from the host buffer in testing mode
        case CodePath::CPU:
            std::memset(grid, 0,
                        paddedGridSize[XX] * paddedGridSize[YY] * paddedGridSize[ZZ] * sizeof(ValueType));
            for (const auto& gridValue : gridValues)
            {
                for (int i = 0; i < DIM; i++)
                {
                    GMX_RELEASE_ASSERT((0 <= gridValue.first[i]) && (gridValue.first[i] < gridSize[i]),
                                       "Invalid grid value index");
                }
                const size_t gridValueIndex =
                        pmeGetGridPlainIndexInternal(gridValue.first, paddedGridSize, gridOrdering);
                grid[gridValueIndex] = gridValue.second;
            }
            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Setting real grid to be used in gather
void pmeSetRealGrid(const gmx_pme_t* pme, CodePath mode, const SparseRealGridValuesInput& gridValues)
{
    pmeSetGridInternal<real>(pme, mode, GridOrdering::XYZ, gridValues);
}

//! Setting complex grid to be used in solve
void pmeSetComplexGrid(const gmx_pme_t*                    pme,
                       CodePath                            mode,
                       GridOrdering                        gridOrdering,
                       const SparseComplexGridValuesInput& gridValues)
{
    pmeSetGridInternal<t_complex>(pme, mode, gridOrdering, gridValues);
}

//! Getting the single dimension's spline values or derivatives
SplineParamsDimVector pmeGetSplineData(const gmx_pme_t* pme, CodePath mode, PmeSplineDataType type, int dimIndex)
{
    GMX_RELEASE_ASSERT(pme != nullptr, "PME data is not initialized");
    const PmeAtomComm* atc       = &(pme->atc[0]);
    const size_t       atomCount = atc->numAtoms();
    const size_t       pmeOrder  = pme->pme_order;
    const size_t       dimSize   = pmeOrder * atomCount;

    real*                 sourceBuffer = pmeGetSplineDataInternal(pme, type, dimIndex);
    SplineParamsDimVector result;
    switch (mode)
    {
        case CodePath::GPU:
            pme_gpu_transform_spline_atom_data(pme->gpu, atc, type, dimIndex, PmeLayoutTransform::GpuToHost);
            result = arrayRefFromArray(sourceBuffer, dimSize);
            break;

        case CodePath::CPU: result = arrayRefFromArray(sourceBuffer, dimSize); break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
    return result;
}

//! Getting the gridline indices
GridLineIndicesVector pmeGetGridlineIndices(const gmx_pme_t* pme, CodePath mode)
{
    GMX_RELEASE_ASSERT(pme != nullptr, "PME data is not initialized");
    const PmeAtomComm* atc       = &(pme->atc[0]);
    const size_t       atomCount = atc->numAtoms();

    GridLineIndicesVector gridLineIndices;
    switch (mode)
    {
        case CodePath::GPU:
            gridLineIndices = arrayRefFromArray(
                    reinterpret_cast<IVec*>(pme_gpu_staging(pme->gpu).h_gridlineIndices), atomCount);
            break;

        case CodePath::CPU: gridLineIndices = atc->idx; break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
    return gridLineIndices;
}

//! Getting real or complex grid - only non zero values
template<typename ValueType>
static SparseGridValuesOutput<ValueType> pmeGetGridInternal(const gmx_pme_t* pme,
                                                            CodePath         mode,
                                                            GridOrdering     gridOrdering)
{
    IVec       gridSize(0, 0, 0), paddedGridSize(0, 0, 0);
    ValueType* grid;
    pmeGetGridAndSizesInternal<ValueType>(pme, mode, grid, gridSize, paddedGridSize);
    SparseGridValuesOutput<ValueType> gridValues;
    switch (mode)
    {
        case CodePath::GPU: // intentional absence of break
        case CodePath::CPU:
            gridValues.clear();
            for (int ix = 0; ix < gridSize[XX]; ix++)
            {
                for (int iy = 0; iy < gridSize[YY]; iy++)
                {
                    for (int iz = 0; iz < gridSize[ZZ]; iz++)
                    {
                        IVec         temp(ix, iy, iz);
                        const size_t gridValueIndex =
                                pmeGetGridPlainIndexInternal(temp, paddedGridSize, gridOrdering);
                        const ValueType value = grid[gridValueIndex];
                        if (value != ValueType{})
                        {
                            auto key        = formatString("Cell %d %d %d", ix, iy, iz);
                            gridValues[key] = value;
                        }
                    }
                }
            }
            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
    return gridValues;
}

//! Getting the real grid (spreading output of pmePerformSplineAndSpread())
SparseRealGridValuesOutput pmeGetRealGrid(const gmx_pme_t* pme, CodePath mode)
{
    return pmeGetGridInternal<real>(pme, mode, GridOrdering::XYZ);
}

//! Getting the complex grid output of pmePerformSolve()
SparseComplexGridValuesOutput pmeGetComplexGrid(const gmx_pme_t* pme, CodePath mode, GridOrdering gridOrdering)
{
    return pmeGetGridInternal<t_complex>(pme, mode, gridOrdering);
}

//! Getting the reciprocal energy and virial
PmeOutput pmeGetReciprocalEnergyAndVirial(const gmx_pme_t* pme, CodePath mode, PmeSolveAlgorithm method)
{
    PmeOutput output;
    switch (mode)
    {
        case CodePath::CPU:
            switch (method)
            {
                case PmeSolveAlgorithm::Coulomb:
                    get_pme_ener_vir_q(pme->solve_work, pme->nthread, &output);
                    break;

                case PmeSolveAlgorithm::LennardJones:
                    get_pme_ener_vir_lj(pme->solve_work, pme->nthread, &output);
                    break;

                default: GMX_THROW(InternalError("Test not implemented for this mode"));
            }
            break;
        case CodePath::GPU:
            switch (method)
            {
                case PmeSolveAlgorithm::Coulomb: pme_gpu_getEnergyAndVirial(*pme, &output); break;

                default: GMX_THROW(InternalError("Test not implemented for this mode"));
            }
            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
    return output;
}

} // namespace test
} // namespace gmx
