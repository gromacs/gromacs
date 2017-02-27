/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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

#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-gather.h"
#include "gromacs/ewald/pme-grid.h"
#include "gromacs/ewald/pme-internal.h"
#include "gromacs/ewald/pme-spread.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace test
{

//! PME initialization - internal
static PmeSafePointer pmeInitInternal(const t_inputrec *inputRec, size_t atomCount)
{
    gmx_pme_t *pmeDataRaw = nullptr;
    gmx_pme_init(&pmeDataRaw, nullptr, 1, 1, inputRec,
                 atomCount, false, false, true, 0.0, 0.0, 1);
    PmeSafePointer pme(pmeDataRaw); // taking ownership
    return pme;
}

//! Simple PME initialization based on input, no atom data
PmeSafePointer pmeInitEmpty(const t_inputrec *inputRec)
{
    return pmeInitInternal(inputRec, 0);
    // hiding the fact that PME actually needs to know the number of atoms in advance
}

//! PME initialization with atom data and system box
PmeSafePointer pmeInitWithAtoms(const t_inputrec        *inputRec,
                                const CoordinatesVector &coordinates,
                                const ChargesVector     &charges,
                                const Matrix3x3          box
                                )
{
    const size_t    atomCount = coordinates.size();
    GMX_RELEASE_ASSERT(atomCount == charges.size(), "Mismatch in atom data");
    PmeSafePointer  pmeSafe = pmeInitInternal(inputRec, atomCount);
    pme_atomcomm_t *atc     = &(pmeSafe->atc[0]);
    atc->x           = const_cast<rvec *>(as_rvec_array(coordinates.data()));
    atc->coefficient = const_cast<real *>(charges.data());
    /* With decomposition there would be more boilerplate atc code here, e.g. do_redist_pos_coeffs */

    // TODO get rid of this with proper matrix type
    matrix boxTemp;
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            boxTemp[i][j] = box[i * DIM + j];
        }
    }
    const char *boxError = check_box(-1, boxTemp);
    GMX_RELEASE_ASSERT(boxError == nullptr, boxError);
    invertBoxMatrix(boxTemp, pmeSafe->recipbox);

    return pmeSafe;
}

//! Getting local PME real grid pointer for test I/O
static real *pmeGetRealGridInternal(const gmx_pme_t *pme)
{
    const size_t gridIndex = 0;
    return pme->fftgrid[gridIndex];
}

//! Getting local PME real grid dimensions
static void pmeGetRealGridSizesInternal(const gmx_pme_t      *pme,
                                        IVec                 &gridSize,
                                        IVec                 &paddedGridSize)
{
    const size_t gridIndex = 0;
    IVec         gridOffsetUnused;
    gmx_parallel_3dfft_real_limits(pme->pfft_setup[gridIndex], gridSize, gridOffsetUnused, paddedGridSize);
}

//! PME spline calculation and charge spreading
void pmePerformSplineAndSpread(gmx_pme_t *pme, CodePath mode, // TODO const qualifiers
                               bool computeSplines, bool spreadCharges)
{
    GMX_RELEASE_ASSERT(pme != nullptr, "PME data is not initialized");
    pme_atomcomm_t *atc                          = &(pme->atc[0]);
    const size_t    gridIndex                    = 0;
    const bool      computeSplinesForZeroCharges = true;
    real           *fftgrid                      = spreadCharges ? pme->fftgrid[gridIndex] : nullptr;

    switch (mode)
    {
        case CodePath::CPU:
            spread_on_grid(pme, atc, &pme->pmegrid[gridIndex], computeSplines, spreadCharges,
                           fftgrid, computeSplinesForZeroCharges, gridIndex);
            if (spreadCharges && !pme->bUseThreads)
            {
                wrap_periodic_pmegrid(pme, pme->pmegrid[gridIndex].grid.grid);
                copy_pmegrid_to_fftgrid(pme, pme->pmegrid[gridIndex].grid.grid, fftgrid, gridIndex);
            }
            break;

        default:
            GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Getting the internal spline data buffer pointer
static real *pmeGetSplineDataInternal(const gmx_pme_t *pme, PmeSplineDataType type, int dimIndex)
{
    GMX_RELEASE_ASSERT((0 <= dimIndex) && (dimIndex < DIM), "Invalid dimension index");
    const pme_atomcomm_t *atc          = &(pme->atc[0]);
    const size_t          threadIndex  = 0;
    real                 *splineBuffer = nullptr;
    switch (type)
    {
        case PmeSplineDataType::Values:
            splineBuffer = atc->spline[threadIndex].theta[dimIndex];
            break;

        case PmeSplineDataType::Derivatives:
            splineBuffer = atc->spline[threadIndex].dtheta[dimIndex];
            break;

        default:
            GMX_THROW(InternalError("Unknown spline data type"));
    }
    return splineBuffer;
}

//! PME force gathering
void pmePerformGather(gmx_pme_t *pme, CodePath mode,
                      PmeGatherInputHandling inputTreatment, ForcesVector &forces)
{
    pme_atomcomm_t *atc                     = &(pme->atc[0]);
    const size_t    atomCount               = atc->n;
    GMX_RELEASE_ASSERT(forces.size() == atomCount, "Bad force buffer size");
    const bool      forceReductionWithInput = (inputTreatment == PmeGatherInputHandling::ReduceWith);
    const real      scale                   = 1.0;
    const size_t    threadIndex             = 0;
    const size_t    gridIndex               = 0;
    real           *grid                    = pme->pmegrid[gridIndex].grid.grid;
    switch (mode)
    {
        case CodePath::CPU:
            atc->f = as_rvec_array(forces.begin());
            if (atc->nthread == 1)
            {
                // something which is normally done in serial spline computation (make_thread_local_ind())
                atc->spline[threadIndex].n = atomCount;
            }
            copy_fftgrid_to_pmegrid(pme, pme->fftgrid[gridIndex], grid, gridIndex, pme->nthread, threadIndex);
            unwrap_periodic_pmegrid(pme, grid);
            gather_f_bsplines(pme, grid, !forceReductionWithInput, atc, &atc->spline[threadIndex], scale);
            break;

        default:
            GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Setting atom spline values/derivatives to be used in spread/gather
void pmeSetSplineData(const gmx_pme_t *pme, CodePath mode,
                      const SplineParamsDimVector &splineValues, PmeSplineDataType type, int dimIndex)
{
    const pme_atomcomm_t *atc         = &(pme->atc[0]);
    const size_t          atomCount   = atc->n;
    const size_t          pmeOrder    = pme->pme_order;
    const size_t          dimSize     = pmeOrder * atomCount;
    GMX_RELEASE_ASSERT(dimSize == splineValues.size(), "Mismatch in spline data");
    real                 *splineBuffer = pmeGetSplineDataInternal(pme, type, dimIndex);

    switch (mode)
    {
        case CodePath::CPU:
            std::copy(splineValues.begin(), splineValues.end(), splineBuffer);
            break;

        default:
            GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Setting gridline indices to be used in spread/gather
void pmeSetGridLineIndices(const gmx_pme_t *pme, CodePath mode,
                           const GridLineIndicesVector &gridLineIndices)
{
    const pme_atomcomm_t       *atc         = &(pme->atc[0]);
    const size_t                atomCount   = atc->n;
    GMX_RELEASE_ASSERT(atomCount == gridLineIndices.size(), "Mismatch in gridline indices size");

    IVec paddedGridSizeUnused, gridSize;
    pmeGetRealGridSizesInternal(pme, gridSize, paddedGridSizeUnused);

    for (const auto &index : gridLineIndices)
    {
        for (int i = 0; i < DIM; i++)
        {
            GMX_RELEASE_ASSERT((0 <= index[i]) && (index[i] < gridSize[i]), "Invalid gridline index");
        }
    }

    switch (mode)
    {
        case CodePath::CPU:
            // incompatible IVec and ivec assignment?
            //std::copy(gridLineIndices.begin(), gridLineIndices.end(), atc->idx);
            memcpy(atc->idx, gridLineIndices.data(), atomCount * sizeof(gridLineIndices[0]));
            break;

        default:
            GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Setting real grid to be used in gather
void pmeSetRealGrid(const gmx_pme_t                 *pme,
                    CodePath                         mode,
                    const SparseRealGridValuesInput &gridValues)
{
    real *grid = pmeGetRealGridInternal(pme);
    IVec  gridSize, paddedGridSize;
    pmeGetRealGridSizesInternal(pme, gridSize, paddedGridSize);
    switch (mode)
    {
        case CodePath::CPU:
            std::memset(grid, 0, paddedGridSize[XX] * paddedGridSize[YY] * paddedGridSize[ZZ] * sizeof(real));
            for (const auto &gridValue : gridValues)
            {
                for (int i = 0; i < DIM; i++)
                {
                    GMX_RELEASE_ASSERT((0 <= gridValue.first[i]) && (gridValue.first[i] < gridSize[i]), "Invalid grid value index");
                }
                const size_t gridValueIndex = (gridValue.first[XX] * paddedGridSize[YY] + gridValue.first[YY]) * paddedGridSize[ZZ] + gridValue.first[ZZ];
                grid[gridValueIndex] = gridValue.second;
            }
            break;

        default:
            GMX_THROW(InternalError("Test not implemented for this mode"));
    }
}

//! Getting the single dimension's spline values or derivatives
SplineParamsDimVector pmeGetSplineData(const gmx_pme_t *pme, CodePath mode,
                                       PmeSplineDataType type, int dimIndex)
{
    GMX_RELEASE_ASSERT(pme != nullptr, "PME data is not initialized");
    const pme_atomcomm_t    *atc         = &(pme->atc[0]);
    const size_t             atomCount   = atc->n;
    const size_t             pmeOrder    = pme->pme_order;
    const size_t             dimSize     = pmeOrder * atomCount;

    real                    *sourceBuffer = pmeGetSplineDataInternal(pme, type, dimIndex);
    SplineParamsDimVector    result;
    switch (mode)
    {
        case CodePath::CPU:
            result = SplineParamsDimVector::fromArray(sourceBuffer, dimSize);
            break;

        default:
            GMX_THROW(InternalError("Test not implemented for this mode"));
    }
    return result;
}

//! Getting the gridline indices
GridLineIndicesVector pmeGetGridlineIndices(const gmx_pme_t *pme, CodePath mode)
{
    GMX_RELEASE_ASSERT(pme != nullptr, "PME data is not initialized");
    const pme_atomcomm_t *atc         = &(pme->atc[0]);
    const size_t          atomCount   = atc->n;

    GridLineIndicesVector gridLineIndices;
    switch (mode)
    {
        case CodePath::CPU:
            gridLineIndices = GridLineIndicesVector::fromArray(reinterpret_cast<IVec *>(atc->idx), atomCount);
            break;

        default:
            GMX_THROW(InternalError("Test not implemented for this mode"));
    }
    return gridLineIndices;
}

//! Getting the real grid (spreading output of PmePerformSplineAndSpread())
SparseRealGridValuesOutput pmeGetRealGrid(const gmx_pme_t *pme, CodePath mode)
{
    GMX_RELEASE_ASSERT(pme != nullptr, "PME data is not initialized");

    SparseRealGridValuesOutput gridValues;
    IVec                       gridSize, paddedGridSize;

    real                      *grid = pmeGetRealGridInternal(pme);
    pmeGetRealGridSizesInternal(pme, gridSize, paddedGridSize);

    switch (mode)
    {
        case CodePath::CPU:
            gridValues.clear();

            for (int ix = 0; ix < gridSize[XX]; ix++)
            {
                for (int iy = 0; iy < gridSize[YY]; iy++)
                {
                    for (int iz = 0; iz < gridSize[ZZ]; iz++)
                    {
                        const size_t gridValueIndex = (ix * paddedGridSize[YY] + iy) * paddedGridSize[ZZ] + iz;
                        const real   value          = grid[gridValueIndex];
                        if (value != 0.0)
                        {
                            auto key = formatString("Cell %d %d %d", ix, iy, iz);
                            gridValues[key] = value;
                        }
                    }
                }
            }
            break;

        default:
            GMX_THROW(InternalError("Test not implemented for this mode"));
    }
    return gridValues;
}

}
}
