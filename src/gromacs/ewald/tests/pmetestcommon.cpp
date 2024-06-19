/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements common routines for PME tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "pmetestcommon.h"

#include <cstring>

#include <algorithm>
#include <iterator>
#include <utility>

#include <gtest/gtest.h>

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme_coordinate_receiver_gpu.h"
#include "gromacs/ewald/pme_gather.h"
#include "gromacs/ewald/pme_gpu_calculate_splines.h"
#include "gromacs/ewald/pme_gpu_constants.h"
#include "gromacs/ewald/pme_gpu_internal.h"
#include "gromacs/ewald/pme_gpu_staging.h"
#include "gromacs/ewald/pme_gpu_types_host.h"
#include "gromacs/ewald/pme_grid.h"
#include "gromacs/ewald/pme_internal.h"
#include "gromacs/ewald/pme_redistribute.h"
#include "gromacs/ewald/pme_solve.h"
#include "gromacs/ewald/pme_spread.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/boxmatrix.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/locality.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/test_hardware_environment.h"
#include "testutils/testasserts.h"
#include "testutils/testinit.h"

class DeviceContext;
class DeviceStream;
class GpuEventSynchronizer;
struct t_inputrec;

namespace gmx
{
namespace test
{

//! A couple of valid inputs for boxes.
const std::map<std::string, Matrix3x3> c_inputBoxes = {
    { "rect", { { 8.0F, 0.0F, 0.0F, 0.0F, 3.4F, 0.0F, 0.0F, 0.0F, 2.0F } } },
    { "tric", { { 7.0F, 0.0F, 0.0F, 0.0F, 4.1F, 0.0F, 3.5F, 2.0F, 12.2F } } },
};

//! Valid PME orders for testing
std::vector<int> c_inputPmeOrders{ 3, 4, 5 };

MessageStringCollector getSkipMessagesIfNecessary(const t_inputrec& inputRec, const CodePath codePath)
{
    // Note that we can't call GTEST_SKIP() from within this method,
    // because it only returns from the current function. So we
    // collect all the reasons why the test cannot run, return them
    // and skip in a higher stack frame.

    MessageStringCollector messages;
    messages.startContext("Test is being skipped because:");

    if (codePath == CodePath::CPU)
    {
        // Everything is implemented, no reason to skip
        return messages;
    }

    std::string errorMessage;
    messages.appendIf(!pme_gpu_supports_build(&errorMessage), errorMessage);
    messages.appendIf(!pme_gpu_supports_input(inputRec, &errorMessage), errorMessage);
    return messages;
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
PmeSafePointer pmeInitWrapper(const t_inputrec*    inputRec,
                              const CodePath       mode,
                              const DeviceContext* deviceContext,
                              const DeviceStream*  deviceStream,
                              const PmeGpuProgram* pmeGpuProgram,
                              const Matrix3x3&     box,
                              const real           ewaldCoeff_q,
                              const real           ewaldCoeff_lj)
{
    const MDLogger dummyLogger;
    const matrix   dummyBox = { { 0 } };
    const auto     runMode  = (mode == CodePath::CPU) ? PmeRunMode::CPU : PmeRunMode::Mixed;
    t_commrec      dummyCommrec;
    NumPmeDomains  numPmeDomains = { 1, 1 };
    // TODO: Need to use proper value when GPU PME decomposition code path is tested
    const real     haloExtentForAtomDisplacement = 1.0;
    gmx_pme_t*     pmeDataRaw                    = gmx_pme_init(&dummyCommrec,
                                         numPmeDomains,
                                         inputRec,
                                         dummyBox,
                                         haloExtentForAtomDisplacement,
                                         false,
                                         false,
                                         true,
                                         ewaldCoeff_q,
                                         ewaldCoeff_lj,
                                         1,
                                         runMode,
                                         nullptr,
                                         deviceContext,
                                         deviceStream,
                                         pmeGpuProgram,
                                         dummyLogger,
                                         nullptr);
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

PmeSafePointer pmeInitEmpty(const t_inputrec* inputRec)
{
    const Matrix3x3 defaultBox = { { 1.0F, 0.0F, 0.0F, 0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 1.0F } };
    return pmeInitWrapper(inputRec, CodePath::CPU, nullptr, nullptr, nullptr, defaultBox, 0.0F, 0.0F);
}

//! Make a GPU state-propagator manager
std::unique_ptr<StatePropagatorDataGpu> makeStatePropagatorDataGpu(const gmx_pme_t&     pme,
                                                                   const DeviceContext* deviceContext,
                                                                   const DeviceStream* deviceStream)
{
    // TODO: Pin the host buffer and use async memory copies
    // TODO: Special constructor for PME-only rank / PME-tests is used here. There should be a mechanism to
    //       restrict one from using other constructor here.
    return std::make_unique<StatePropagatorDataGpu>(
            deviceStream, *deviceContext, GpuApiCallBehavior::Sync, pme_gpu_get_block_size(&pme), nullptr);
}

//! PME initialization with atom data
void pmeInitAtoms(gmx_pme_t*               pme,
                  StatePropagatorDataGpu*  stateGpu,
                  const CodePath           mode,
                  const CoordinatesVector& coordinates,
                  const ChargesVector&     charges)
{
    const Index atomCount = coordinates.size();
    GMX_RELEASE_ASSERT(atomCount == gmx::ssize(charges), "Mismatch in atom data");
    PmeAtomComm* atc = nullptr;

    switch (mode)
    {
        case CodePath::CPU:
            atc              = &(pme->atc[0]);
            atc->x           = coordinates;
            atc->coefficient = charges;
            gmx_pme_reinit_atoms(pme, atomCount, charges, {});
            /* With decomposition there would be more boilerplate atc code here, e.g. do_redist_pos_coeffs */
            break;

        case CodePath::GPU:
            // TODO: Avoid use of atc in the GPU code path
            atc = &(pme->atc[0]);
            // We need to set atc->n for passing the size in the tests
            atc->setNumAtoms(atomCount);
            gmx_pme_reinit_atoms(pme, atomCount, charges, {});

            stateGpu->reinit(atomCount, atomCount);
            stateGpu->copyCoordinatesToGpu(arrayRefFromArray(coordinates.data(), coordinates.size()),
                                           gmx::AtomLocality::Local);
            pme_gpu_set_kernelparam_coordinates(pme->gpu, stateGpu->getCoordinates());

            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Getting local PME real grid pointer for test I/O
static real* pmeGetRealGridInternal(const gmx_pme_t* pme)
{
    const size_t gridIndex = 0;
    return pme->gridsCoulomb[gridIndex].fftgrid;
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
            gmx_parallel_3dfft_real_limits(
                    pme->gridsCoulomb[gridIndex].pfft_setup.get(), gridSize, gridOffsetUnused, paddedGridSize);
            break;

        case CodePath::GPU:
            pme_gpu_get_real_grid_sizes(pme->gpu, &gridSize, &paddedGridSize);
            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Getting local PME complex grid pointer for test I/O
static t_complex* pmeGetComplexGridInternal(gmx_pme_t* pme)
{
    const size_t gridIndex = 0;
    return pme->gridsCoulomb[gridIndex].cfftgrid;
}

//! Getting local PME complex grid dimensions
static void pmeGetComplexGridSizesInternal(const gmx_pme_t* pme,
                                           IVec& gridSize,       //NOLINT(google-runtime-references)
                                           IVec& paddedGridSize) //NOLINT(google-runtime-references)
{
    const size_t gridIndex = 0;
    IVec         gridOffsetUnused, complexOrderUnused;
    gmx_parallel_3dfft_complex_limits(pme->gridsCoulomb[gridIndex].pfft_setup.get(),
                                      complexOrderUnused,
                                      gridSize,
                                      gridOffsetUnused,
                                      paddedGridSize); // TODO: what about YZX ordering?
}

//! Getting the PME grid memory buffer and its sizes - template definition
template<typename ValueType>
static void pmeGetGridAndSizesInternal(gmx_pme_t* /*unused*/,
                                       CodePath /*unused*/,
                                       ValueType*& /*unused*/, //NOLINT(google-runtime-references)
                                       IVec& /*unused*/,       //NOLINT(google-runtime-references)
                                       IVec& /*unused*/) = delete; //NOLINT(google-runtime-references)

//! Getting the PME real grid memory buffer and its sizes
template<>
void pmeGetGridAndSizesInternal<real>(gmx_pme_t* pme, CodePath mode, real*& grid, IVec& gridSize, IVec& paddedGridSize)
{
    grid = pmeGetRealGridInternal(pme);
    pmeGetRealGridSizesInternal(pme, mode, gridSize, paddedGridSize);
}

//! Getting the PME complex grid memory buffer and its sizes
template<>
void pmeGetGridAndSizesInternal<t_complex>(gmx_pme_t* pme,
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
    PmeAtomComm*    atc                          = &(pme->atc[0]);
    const size_t    gridIndex                    = 0;
    PmeAndFftGrids& grids                        = pme->gridsCoulomb[gridIndex];
    const bool      computeSplinesForZeroCharges = true;

    switch (mode)
    {
        case CodePath::CPU:
            spread_on_grid(pme, atc, &grids, computeSplines, spreadCharges, computeSplinesForZeroCharges);
            if (spreadCharges && !pme->bUseThreads)
            {
                wrap_periodic_pmegrid(pme, grids.pmeGrids.grid.grid);
                copy_pmegrid_to_fftgrid(pme, &grids);
            }
            break;

/* The compiler will complain about passing fftgrid (converting double ** to float **) if using
 * double precision. GPUs are not used with double precision anyhow. */
#if !GMX_DOUBLE
        case CodePath::GPU:
        {
            const real lambdaQ = 1.0;
            // no synchronization needed as x is transferred in the PME stream
            GpuEventSynchronizer* xReadyOnDevice = nullptr;

            bool                           useGpuDirectComm         = false;
            gmx::PmeCoordinateReceiverGpu* pmeCoordinateReceiverGpu = nullptr;

            pme_gpu_spread(pme->gpu,
                           xReadyOnDevice,
                           pme->gridsCoulomb,
                           computeSplines,
                           spreadCharges,
                           lambdaQ,
                           useGpuDirectComm,
                           pmeCoordinateReceiverGpu,
                           false,
                           nullptr);
        }
        break;
#endif

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
void pmePerformSolve(gmx_pme_t*        pme,
                     CodePath          mode,
                     PmeSolveAlgorithm method,
                     real              cellVolume,
                     GridOrdering      gridOrdering,
                     bool              computeEnergyAndVirial)
{
    t_complex*   h_grid              = pmeGetComplexGridInternal(pme);
    const bool   useLorentzBerthelot = false;
    const size_t threadIndex         = 0;
    const size_t gridIndex           = 0;
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
                    pme->pmeSolve->solveCoulombYZX(
                            *pme, h_grid, cellVolume, computeEnergyAndVirial, threadIndex);
                    break;

                case PmeSolveAlgorithm::LennardJones:
                    // Swap the complex grid to be passed to solve with h_grid
                    std::swap(pme->gridsLJ[0].cfftgrid, h_grid);
                    pme->pmeSolve->solveLJYZX(
                            *pme, pme->gridsLJ, useLorentzBerthelot, cellVolume, computeEnergyAndVirial, threadIndex);
                    // Swap back
                    std::swap(pme->gridsLJ[0].cfftgrid, h_grid);
                    break;

                default: GMX_THROW(InternalError("Test not implemented for this mode"));
            }
            break;

        case CodePath::GPU:
            switch (method)
            {
                case PmeSolveAlgorithm::Coulomb:
                    pme_gpu_solve(pme->gpu, gridIndex, h_grid, gridOrdering, computeEnergyAndVirial);
                    break;

                default: GMX_THROW(InternalError("Test not implemented for this mode"));
            }
            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! PME force gathering
void pmePerformGather(gmx_pme_t* pme, CodePath mode, ForcesVector& forces)
{
    PmeAtomComm* atc       = &(pme->atc[0]);
    const Index  atomCount = atc->numAtoms();
    GMX_RELEASE_ASSERT(forces.ssize() == atomCount, "Invalid force buffer size");
    const real      scale       = 1.0;
    const size_t    threadIndex = 0;
    const size_t    gridIndex   = 0;
    PmeAndFftGrids& grids       = pme->gridsCoulomb[gridIndex];
    ArrayRef<real>  pmegrid     = grids.pmeGrids.grid.grid;

    switch (mode)
    {
        case CodePath::CPU:
            atc->f = forces;
            if (atc->nthread == 1)
            {
                // something which is normally done in serial spline computation (make_thread_local_ind())
                atc->spline[threadIndex].n = atomCount;
            }
            copy_fftgrid_to_pmegrid(pme, &grids, pme->nthread, threadIndex);
            unwrap_periodic_pmegrid(pme, pmegrid);
            gather_f_bsplines(pme, pmegrid, true, atc, &atc->spline[threadIndex], scale);
            break;

/* The compiler will complain about passing fftgrid (converting double ** to float **) if using
 * double precision. GPUs are not used with double precision anyhow. */
#if !GMX_DOUBLE
        case CodePath::GPU:
        {
            // Variable initialization needs a non-switch scope
            const bool computeEnergyAndVirial = false;
            const real lambdaQ                = 1.0;
            PmeOutput  output = pme_gpu_getOutput(pme, computeEnergyAndVirial, lambdaQ);
            GMX_ASSERT(forces.size() == output.forces_.size(),
                       "Size of force buffers did not match");
            pme_gpu_gather(pme->gpu, pme->gridsCoulomb, lambdaQ, nullptr, computeEnergyAndVirial);
            std::copy(std::begin(output.forces_), std::end(output.forces_), std::begin(forces));
        }
        break;
#endif

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

//! A binary enum for spline data layout transformation
enum class PmeLayoutTransform
{
    GpuToHost,
    HostToGpu
};

/*! \brief Gets a unique index to an element in a spline parameter buffer.
 *
 * These theta/dtheta buffers are laid out for GPU spread/gather
 * kernels. The index is wrt the execution block, in range(0,
 * atomsPerBlock * order * DIM).
 *
 * This is a wrapper, only used in unit tests.
 * \param[in] order            PME order
 * \param[in] splineIndex      Spline contribution index (from 0 to \p order - 1)
 * \param[in] dimIndex         Dimension index (from 0 to 2)
 * \param[in] atomIndex        Atom index wrt the block.
 * \param[in] atomsPerWarp     Number of atoms processed by a warp.
 *
 * \returns Index into theta or dtheta array using GPU layout.
 */
static int getSplineParamFullIndex(int order, int splineIndex, int dimIndex, int atomIndex, int atomsPerWarp)
{
    if (order != c_pmeGpuOrder)
    {
        throw order;
    }
    constexpr int fixedOrder = c_pmeGpuOrder;
    GMX_UNUSED_VALUE(fixedOrder);

    const int atomWarpIndex = atomIndex % atomsPerWarp;
    const int warpIndex     = atomIndex / atomsPerWarp;
    int       indexBase, result;
    switch (atomsPerWarp)
    {
        case 1:
            indexBase = getSplineParamIndexBase<fixedOrder, 1>(warpIndex, atomWarpIndex);
            result    = getSplineParamIndex<fixedOrder, 1>(indexBase, dimIndex, splineIndex);
            break;

        case 2:
            indexBase = getSplineParamIndexBase<fixedOrder, 2>(warpIndex, atomWarpIndex);
            result    = getSplineParamIndex<fixedOrder, 2>(indexBase, dimIndex, splineIndex);
            break;

        case 4:
            indexBase = getSplineParamIndexBase<fixedOrder, 4>(warpIndex, atomWarpIndex);
            result    = getSplineParamIndex<fixedOrder, 4>(indexBase, dimIndex, splineIndex);
            break;

        case 8:
            indexBase = getSplineParamIndexBase<fixedOrder, 8>(warpIndex, atomWarpIndex);
            result    = getSplineParamIndex<fixedOrder, 8>(indexBase, dimIndex, splineIndex);
            break;

        default:
            GMX_THROW(NotImplementedError(
                    formatString("Test function call not unrolled for atomsPerWarp = %d in "
                                 "getSplineParamFullIndex",
                                 atomsPerWarp)));
    }
    return result;
}

/*! \brief Rearranges the atom spline data between the GPU and host layouts.
 * Only used for test purposes so far, likely to be horribly slow.
 *
 * \param[in]  pmeGpu     The PME GPU structure.
 * \param[out] atc        The PME CPU atom data structure (with a single-threaded layout).
 * \param[in]  type       The spline data type (values or derivatives).
 * \param[in]  dimIndex   Dimension index.
 * \param[in]  transform  Layout transform type
 */
static void pme_gpu_transform_spline_atom_data(PmeGpu*            pmeGpu,
                                               const PmeAtomComm* atc,
                                               PmeSplineDataType  type,
                                               int                dimIndex,
                                               PmeLayoutTransform transform)
{
    // The GPU atom spline data is laid out in a different way currently than the CPU one.
    // This function converts the data from GPU to CPU layout (in the host memory).
    // It is only intended for testing purposes so far.
    // Ideally we should use similar layouts on CPU and GPU if we care about mixed modes and their
    // performance (e.g. spreading on GPU, gathering on CPU).
    GMX_RELEASE_ASSERT(atc->nthread == 1, "Only the serial PME data layout is supported");
    const uintmax_t threadIndex  = 0;
    const auto      atomCount    = atc->numAtoms();
    const auto      atomsPerWarp = pme_gpu_get_atoms_per_warp(pmeGpu);
    GMX_RELEASE_ASSERT(atomsPerWarp > 0, "Can not get GPU warp size");
    const auto pmeOrder = pmeGpu->common->pme_order;
    GMX_ASSERT(pmeOrder == c_pmeGpuOrder, "Only PME order 4 is implemented");

    real*  cpuSplineBuffer;
    float* h_splineBuffer;
    switch (type)
    {
        case PmeSplineDataType::Values:
            cpuSplineBuffer = atc->spline[threadIndex].theta.coefficients[dimIndex];
            h_splineBuffer  = pmeGpu->staging.h_theta.data();
            break;

        case PmeSplineDataType::Derivatives:
            cpuSplineBuffer = atc->spline[threadIndex].dtheta.coefficients[dimIndex];
            h_splineBuffer  = pmeGpu->staging.h_dtheta.data();
            break;

        default: GMX_THROW(InternalError("Unknown spline data type"));
    }

    for (auto atomIndex = 0; atomIndex < atomCount; atomIndex++)
    {
        for (auto orderIndex = 0; orderIndex < pmeOrder; orderIndex++)
        {
            const auto gpuValueIndex =
                    getSplineParamFullIndex(pmeOrder, orderIndex, dimIndex, atomIndex, atomsPerWarp);
            const auto cpuValueIndex = atomIndex * pmeOrder + orderIndex;
            GMX_ASSERT(cpuValueIndex < atomCount * pmeOrder,
                       "Atom spline data index out of bounds (while transforming GPU data layout "
                       "for host)");
            switch (transform)
            {
                case PmeLayoutTransform::GpuToHost:
                    cpuSplineBuffer[cpuValueIndex] = h_splineBuffer[gpuValueIndex];
                    break;

                case PmeLayoutTransform::HostToGpu:
                    h_splineBuffer[gpuValueIndex] = cpuSplineBuffer[cpuValueIndex];
                    break;

                default: GMX_THROW(InternalError("Unknown layout transform"));
            }
        }
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
    const Index        atomCount = atc->numAtoms();
    const Index        pmeOrder  = pme->pme_order;
    const Index        dimSize   = pmeOrder * atomCount;
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
    const Index  atomCount = atc->numAtoms();
    GMX_RELEASE_ASSERT(atomCount == gmx::ssize(gridLineIndices),
                       "Mismatch in gridline indices size");

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
            memcpy(pme_gpu_staging(pme->gpu).h_gridlineIndices.data(),
                   gridLineIndices.data(),
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
static void pmeSetGridInternal(gmx_pme_t*                              pme,
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
            std::memset(grid, 0, paddedGridSize[XX] * paddedGridSize[YY] * paddedGridSize[ZZ] * sizeof(ValueType));
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
void pmeSetRealGrid(gmx_pme_t* pme, CodePath mode, const SparseRealGridValuesInput& gridValues)
{
    pmeSetGridInternal<real>(pme, mode, GridOrdering::XYZ, gridValues);
}

//! Setting complex grid to be used in solve
void pmeSetComplexGrid(gmx_pme_t*                          pme,
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
    const PmeAtomComm* atc       = pme->atc.data();
    const size_t       atomCount = atc->numAtoms();

    GridLineIndicesVector gridLineIndices;
    switch (mode)
    {
        case CodePath::GPU:
        {
            auto* gridlineIndicesAsIVec =
                    reinterpret_cast<IVec*>(pme_gpu_staging(pme->gpu).h_gridlineIndices.data());
            ArrayRef<IVec> gridlineIndicesArrayRef = arrayRefFromArray(gridlineIndicesAsIVec, atomCount);
            gridLineIndices = { gridlineIndicesArrayRef.begin(), gridlineIndicesArrayRef.end() };
        }
        break;

        case CodePath::CPU: gridLineIndices = { atc->idx.begin(), atc->idx.end() }; break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
    return gridLineIndices;
}

//! Getting real or complex grid - only non zero values
template<typename ValueType>
static SparseGridValuesOutput<ValueType> pmeGetGridInternal(gmx_pme_t* pme, CodePath mode, GridOrdering gridOrdering)
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
SparseRealGridValuesOutput pmeGetRealGrid(gmx_pme_t* pme, CodePath mode)
{
    return pmeGetGridInternal<real>(pme, mode, GridOrdering::XYZ);
}

//! Getting the complex grid output of pmePerformSolve()
SparseComplexGridValuesOutput pmeGetComplexGrid(gmx_pme_t* pme, CodePath mode, GridOrdering gridOrdering)
{
    return pmeGetGridInternal<t_complex>(pme, mode, gridOrdering);
}

//! Getting the reciprocal energy and virial
PmeOutput pmeGetReciprocalEnergyAndVirial(const gmx_pme_t* pme, CodePath mode, PmeSolveAlgorithm method)
{
    PmeOutput  output;
    const real lambdaQ = 1.0;
    switch (mode)
    {
        case CodePath::CPU:
            switch (method)
            {
                case PmeSolveAlgorithm::Coulomb:
                    pme->pmeSolve->getCoulombEnergyAndVirial(&output);
                    break;

                case PmeSolveAlgorithm::LennardJones:
                    pme->pmeSolve->getLJEnergyAndVirial(&output);
                    break;

                default: GMX_THROW(InternalError("Test not implemented for this mode"));
            }
            break;
        case CodePath::GPU:
            switch (method)
            {
                case PmeSolveAlgorithm::Coulomb:
                    pme_gpu_getEnergyAndVirial(*pme, lambdaQ, &output);
                    break;

                default: GMX_THROW(InternalError("Test not implemented for this mode"));
            }
            break;

        default: GMX_THROW(InternalError("Test not implemented for this mode"));
    }
    return output;
}

std::string makeRefDataFileName()
{
    // By default, the reference data filename is set via a call to
    // gmx::TestFileManager::getTestSpecificFileName() that queries
    // GoogleTest and gets a string that includes the return value for
    // nameOfTest(). The logic here must match that of the call to
    // ::testing::RegisterTest, so that it works as intended. In
    // particular, the name must include a "WorksOn" substring that
    // precedes the name of the hardware context, so that this can be
    // removed.
    //
    // Get the info about the test
    const ::testing::TestInfo* testInfo = ::testing::UnitTest::GetInstance()->current_test_info();

    // Get the test name and prepare to remove the part describing the
    // hardware context.
    std::string testName(testInfo->name());
    auto        worksOnPos = testName.find("WorksOn");
    GMX_RELEASE_ASSERT(worksOnPos != testName.size(),
                       "Test name must include the 'WorksOn' fragment");

    // Build the complete refdata filename like
    // getTestSpecificFilename() would do it for a non-dynamical
    // parameterized test.
    std::string refDataFileName = formatString("%s_%sWorksWith_%s.xml",
                                               testInfo->test_suite_name(),
                                               testName.substr(0, worksOnPos).c_str(),
                                               testInfo->value_param());
    // Use the check that the name isn't too long
    checkTestNameLength(refDataFileName);
    return refDataFileName;
}

std::string makeHardwareContextName(const int hardwareContextIndex)
{
    std::optional<int> gpuId = getPmeTestHardwareContexts()[hardwareContextIndex].gpuId();
    std::string        description;
    if (gpuId.has_value())
    {
        description = "GPU" + std::to_string(gpuId.value());
    }
    else
    {
        description = "CPU";
    }
    return description;
}

PmeTestHardwareContext::PmeTestHardwareContext() : codePath_(CodePath::CPU) {}

PmeTestHardwareContext::PmeTestHardwareContext(TestDevice* testDevice) :
    codePath_(CodePath::GPU), testDevice_(testDevice)
{
    testDevice_->deviceContext().activate();
    pmeGpuProgram_ = buildPmeGpuProgram(testDevice_->deviceContext());
}

//! Returns a human-readable context description line
std::string PmeTestHardwareContext::description() const
{
    switch (codePath_)
    {
        case CodePath::CPU: return "CPU";
        case CodePath::GPU: return "GPU (" + testDevice_->description() + ")";
        default: return "Unknown code path.";
    }
}

std::optional<int> PmeTestHardwareContext::gpuId() const
{
    switch (codePath_)
    {
        case CodePath::CPU: return std::nullopt;
        case CodePath::GPU: return testDevice_->id();
        default: return std::nullopt;
    }
}

void PmeTestHardwareContext::activate() const
{
    if (codePath_ == CodePath::GPU)
    {
        testDevice_->deviceContext().activate();
    }
}

ArrayRef<const PmeTestHardwareContext> getPmeTestHardwareContexts()
{
    // The test hardware contexts used in PME tests
    static std::vector<PmeTestHardwareContext> s_pmeTestHardwareContexts;

    // Lazily make s_pmeTestHardwareContexts to avoid static
    // initialization fiasco.
    if (s_pmeTestHardwareContexts.empty())
    {
        // Add CPU
        s_pmeTestHardwareContexts.emplace_back();
        // Add GPU devices
        const auto& testDeviceList = getTestHardwareEnvironment()->getTestDeviceList();
        for (const auto& testDevice : testDeviceList)
        {
            s_pmeTestHardwareContexts.emplace_back(testDevice.get());
        }
    }
    return s_pmeTestHardwareContexts;
}

void registerTestsDynamically()
{
    auto       contexts = getPmeTestHardwareContexts();
    Range<int> contextIndexRange(0, contexts.size());
    registerDynamicalPmeSplineSpreadTests(contextIndexRange);
    registerDynamicalPmeSolveTests(contextIndexRange);
    registerDynamicalPmeGatherTests(contextIndexRange);
}

} // namespace test
} // namespace gmx
