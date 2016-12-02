/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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

#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-grid.h"
#include "gromacs/ewald/pme-internal.h"
#include "gromacs/ewald/pme-spread.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{
//! PME initialization - internal
PmeSafePointer PmeInitInternal(const t_inputrec *inputRec, size_t atomCount)
{
    gmx_pme_t *pmeDataRaw = NULL;
    gmx_pme_init(&pmeDataRaw, NULL, 1, 1, inputRec,
                 atomCount, FALSE, FALSE, TRUE, 0.0, 0.0, 1);
    PmeSafePointer pme(pmeDataRaw, gmx_pme_destroy); // taking ownership
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
            // TODO: there could be 2 GPU modes - separate spline and spread, and fused.
            // In the fused case everything faster in a single kernel.
            // Everything still ends up on the host in the end - have to make sure of that!
    }
}

//! Fetching the spline computation outputs of PmePerformSplineAndSpread()
void PmeFetchOutputsSpline(const PmeSafePointer &pmeSafe, PmeCodePath mode,
                           CoordinatesVector &fractCoordinates,
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
            fractCoordinates.assign(atc->fractx, atc->fractx + atomCount);
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
                           SparseGridValues &gridValues)
{
    const size_t          gridIndex = 0;
    IVec                  gridSize, gridOffsetUnused, paddedGridSize;

    switch (mode)
    {
        case PmeCodePath::CPU:
            gridValues.clear();

            gmx_parallel_3dfft_real_limits(pmeSafe->pfft_setup[gridIndex], gridSize, gridOffsetUnused, paddedGridSize);
            for (int ix = 0; ix < gridSize[XX]; ix++)
            {
                for (int iy = 0; iy < gridSize[YY]; iy++)
                {
                    for (int iz = 0; iz < gridSize[ZZ]; iz++)
                    {
                        const size_t gridValueIndex = (ix * paddedGridSize[YY] + iy) * paddedGridSize[ZZ] + iz;
                        const real   value          = pmeSafe->fftgrid[gridIndex][gridValueIndex];
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
