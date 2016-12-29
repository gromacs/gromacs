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
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

//! Getting local PME real grid pointer for test I/O
real *PmeGetRealGrid(const PmeSafePointer &pmeSafe)
{
    const size_t gridIndex = 0;
    return pmeSafe->fftgrid[gridIndex];
}

//! Getting local PME real grid dimensions
void PmeGetRealGridSizes(const PmeSafePointer &pmeSafe,
                         IVec                 &gridSize,
                         IVec                 &paddedGridSize)
{
    const size_t gridIndex = 0;
    IVec         gridOffsetUnused;
    gmx_parallel_3dfft_real_limits(pmeSafe->pfft_setup[gridIndex], gridSize, gridOffsetUnused, paddedGridSize);
}

//! PME initialization - internal
PmeSafePointer PmeInitInternal(const t_inputrec *inputRec, size_t atomCount)
{
    gmx_pme_t *pmeDataRaw = NULL;
    gmx_pme_init(&pmeDataRaw, NULL, 1, 1, inputRec,
                 atomCount, FALSE, FALSE, TRUE, 0.0, 0.0, 1);
    PmeSafePointer pme(pmeDataRaw); // taking ownership
    return pme;
}

//! Simple PME initialization based on input, no atom data
PmeSafePointer PmeInitEmpty(const t_inputrec *inputRec)
{
    return PmeInitInternal(inputRec, 0);
    // hiding the fact that PME actually needs to know the number of atoms in advance
}

//! PME initialization with atom data and system box
PmeSafePointer PmeInitWithAtoms(const t_inputrec        *inputRec,
                                const CoordinatesVector &coordinates,
                                const ChargesVector     &charges,
                                const Matrix3x3          box
                                )
{
    const size_t    atomCount = coordinates.size();
    GMX_RELEASE_ASSERT(atomCount == charges.size(), "Mismatch in atom data");
    PmeSafePointer  pmeSafe = PmeInitInternal(inputRec, atomCount);
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
    gmx::invertBoxMatrix(boxTemp, pmeSafe->recipbox);

    return pmeSafe;
}

//! PME spline calculation and charge spreading
void PmePerformSplineAndSpread(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                               bool computeSplines, bool spreadCharges)
{
    gmx_pme_t      *pmeUnsafe                    = pmeSafe.get();
    pme_atomcomm_t *atc                          = &(pmeSafe->atc[0]);
    const size_t    gridIndex                    = 0;
    const bool      computeSplinesForZeroCharges = true;
    real           *fftgrid                      = spreadCharges ? pmeSafe->fftgrid[gridIndex] : nullptr;

    switch (mode)
    {
        case PmeCodePath::CPU:
            spread_on_grid(pmeUnsafe, atc, &pmeSafe->pmegrid[gridIndex], computeSplines, spreadCharges,
                           fftgrid, computeSplinesForZeroCharges, gridIndex);
            if (spreadCharges && !pmeSafe->bUseThreads)
            {
                wrap_periodic_pmegrid(pmeUnsafe, pmeSafe->pmegrid[gridIndex].grid.grid);
                copy_pmegrid_to_fftgrid(pmeUnsafe, pmeSafe->pmegrid[gridIndex].grid.grid, fftgrid, gridIndex);
            }
            break;

        default:
            GMX_THROW(gmx::InternalError("Test not implemented for this mode"));
    }
}

//! PME force gathering
void PmePerformGather(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                      PmeGatherInputHandling inputTreatment, ForcesVector &forces)
{
    gmx_pme_t      *pmeUnsafe               = pmeSafe.get();
    pme_atomcomm_t *atc                     = &(pmeSafe->atc[0]);
    const size_t    atomCount               = atc->n;
    GMX_RELEASE_ASSERT(forces.size() == atomCount, "Bad force buffer size");
    const bool      forceReductionWithInput = (inputTreatment == PmeGatherInputHandling::ReduceWith);
    const real      scale                   = 1.0;
    const size_t    threadIndex             = 0;
    const size_t    gridIndex               = 0;
    real           *grid                    = pmeSafe->pmegrid[gridIndex].grid.grid;
    switch (mode)
    {
        case PmeCodePath::CPU:
            atc->f = as_rvec_array(forces.data());
            if (atc->nthread == 1)
            {
                // something which is normally done in serial spline computation (make_thread_local_ind())
                atc->spline[threadIndex].n = atomCount;
            }
            copy_fftgrid_to_pmegrid(pmeUnsafe, pmeSafe->fftgrid[gridIndex], grid, gridIndex, pmeSafe->nthread, threadIndex);
            unwrap_periodic_pmegrid(pmeUnsafe, grid);
            gather_f_bsplines(pmeUnsafe, grid, !forceReductionWithInput,
                              atc, &atc->spline[threadIndex], scale);
            break;

        default:
            GMX_THROW(gmx::InternalError("Test not implemented for this mode"));
    }
}

//! Setting atom spline values/derivatives to be used in spread/gather
void PmeSetSplineData(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                      const SplineParamsVector &splineValues, PmeSplineDataType type)
{
    const pme_atomcomm_t *atc         = &(pmeSafe->atc[0]);
    const size_t          atomCount   = atc->n;
    const size_t          pmeOrder    = pmeSafe->pme_order;
    const size_t          dimSize     = pmeOrder * atomCount;
    GMX_RELEASE_ASSERT(DIM * dimSize == splineValues.size(), "Mismatch in spline data");
    const size_t          threadIndex   = 0;
    real                **targetBuffers = nullptr;
    switch (type)
    {
        case PmeSplineDataType::Values:
            targetBuffers = atc->spline[threadIndex].theta;
            break;

        case PmeSplineDataType::Derivatives:
            targetBuffers = atc->spline[threadIndex].dtheta;
            break;

        default:
            GMX_THROW(InternalError("Unknown spline data type"));
    }

    switch (mode)
    {
        case PmeCodePath::CPU:
            // spline values - XX...XXYY...YYZZ...ZZ
            for (int i = 0; i < DIM; i++)
            {
                std::copy(splineValues.begin() + i * dimSize,
                          splineValues.begin() + (i + 1) * dimSize,
                          targetBuffers[i]);
            }
            break;

        default:
            GMX_THROW(gmx::InternalError("Test not implemented for this mode"));
    }
}

//! Setting gridline indices to be used in spread/gather
void PmeSetGridLineIndices(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                           const GridLineIndicesVector &gridLineIndices)
{
    pme_atomcomm_t       *atc         = &(pmeSafe->atc[0]);
    const size_t          atomCount   = atc->n;
    GMX_RELEASE_ASSERT(atomCount == gridLineIndices.size(), "Mismatch in gridline indices size");

    IVec paddedGridSizeUnused, gridSize;
    PmeGetRealGridSizes(pmeSafe, gridSize, paddedGridSizeUnused);
    for (const auto &index : gridLineIndices)
    {
        for (int i = 0; i < DIM; i++)
        {
            GMX_RELEASE_ASSERT((0 <= index[i]) && (index[i] < gridSize[i]), "Invalid gridline index");
        }
    }

    switch (mode)
    {
        case PmeCodePath::CPU:
            // incompatible IVec and ivec assignment?
            //std::copy(gridLineIndices.begin(), gridLineIndices.end(), atc->idx);
            memcpy(atc->idx, gridLineIndices.data(), atomCount * sizeof(gridLineIndices[0]));
            break;

        default:
            GMX_THROW(gmx::InternalError("Test not implemented for this mode"));
    }
}

//! Setting real grid to be used in gather
void PmeSetRealGrid(const PmeSafePointer       &pmeSafe,
                    PmeCodePath                 mode,
                    const SparseRealGridValues &gridValues)
{
    real *grid = PmeGetRealGrid(pmeSafe);
    IVec  gridSize, paddedGridSize;
    PmeGetRealGridSizes(pmeSafe, gridSize, paddedGridSize);
    switch (mode)
    {
        case PmeCodePath::CPU:
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
            GMX_THROW(gmx::InternalError("Test not implemented for this mode"));
    }
}

//! Fetching the spline computation outputs of PmePerformSplineAndSpread()
void PmeFetchOutputsSpline(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                           SplineParamsVector &splineValues,
                           SplineParamsVector &splineDerivatives,
                           GridLineIndicesVector &gridLineIndices)

{
    const pme_atomcomm_t *atc         = &(pmeSafe->atc[0]);
    const size_t          atomCount   = atc->n;
    const size_t          pmeOrder    = pmeSafe->pme_order;
    const size_t          threadIndex = 0; // relying on running single threaded on CPU
    switch (mode)
    {
        case PmeCodePath::CPU:
            gridLineIndices.assign(atc->idx, atc->idx + atomCount);
            // spline values - XX...XXYY...YYZZ...ZZ
            splineValues.clear();
            splineDerivatives.clear();
            for (int i = 0; i < DIM; i++)
            {
                splineValues.insert(splineValues.end(), atc->spline[threadIndex].theta[i],
                                    atc->spline[threadIndex].theta[i] + atomCount * pmeOrder);
                splineDerivatives.insert(splineDerivatives.end(), atc->spline[threadIndex].dtheta[i],
                                         atc->spline[threadIndex].dtheta[i] + atomCount * pmeOrder);
            }
            break;

        default:
            GMX_THROW(gmx::InternalError("Test not implemented for this mode"));
    }
}

//! Fetching the spreading output of PmePerformSplineAndSpread()
void PmeFetchOutputsSpread(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                           SparseRealGridValues &gridValues)
{
    real *grid = PmeGetRealGrid(pmeSafe);
    IVec  gridSize, paddedGridSize;
    PmeGetRealGridSizes(pmeSafe, gridSize, paddedGridSize);

    switch (mode)
    {
        case PmeCodePath::CPU:
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
                            IVec key = {ix, iy, iz};
                            gridValues[key] = value;
                        }
                    }
                }
            }
            break;

        default:
            GMX_THROW(gmx::InternalError("Test not implemented for this mode"));
    }
}

}
