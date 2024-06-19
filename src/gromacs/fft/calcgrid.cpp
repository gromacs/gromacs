/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include "gmxpre.h"

#include "calcgrid.h"

#include <cmath>

#include <algorithm>
#include <filesystem>

#include "gromacs/utility/fatalerror.h"

/* The grid sizes below are based on timing of a 3D cubic grid in fftw
 * compiled with SSE using 4 threads in fft5d.c.
 * A grid size is removed when a larger grid is faster.
 */

/* Small grid size array */
constexpr int g_initNR            = 15;
constexpr int grid_init[g_initNR] = { 6, 8, 10, 12, 14, 16, 20, 24, 25, 28, 32, 36, 40, 42, 44 };

/* For larger grid sizes, a prefactor with any power of 2 can be added.
 * Only sizes divisible by 4 should be used, 90 is allowed, 140 not.
 */
constexpr int g_baseNR            = 14;
constexpr int grid_base[g_baseNR] = { 45, 48, 50, 52, 54, 56, 60, 64, 70, 72, 75, 80, 81, 84 };

real calcFftGrid(FILE* fp, const matrix box, real gridSpacing, int minGridPointsPerDim, int* nx, int* ny, int* nz)
{
    int  d, n[DIM];
    int  i;
    rvec box_size;
    int  nmin, fac2, attempt;
    rvec spacing;
    real max_spacing;

    if ((*nx <= 0 || *ny <= 0 || *nz <= 0) && gridSpacing <= 0)
    {
        gmx_fatal(FARGS, "invalid fourier grid spacing: %g", gridSpacing);
    }

    static_assert(grid_base[g_baseNR - 1] % 4 == 0,
                  "the last entry in grid_base is not a multiple of 4");

    /* New grid calculation setup:
     *
     * To maintain similar accuracy for triclinic PME grids as for rectangular
     * ones, the max grid spacing should set along the box vectors rather than
     * cartesian X/Y/Z directions. This will lead to slightly larger grids, but
     * it is much better than having to go to pme_order=6.
     *
     * Thus, instead of just extracting the diagonal elements to box_size[d], we
     * now calculate the cartesian length of the vectors.
     *
     * /Erik Lindahl, 20060402.
     */
    for (d = 0; d < DIM; d++)
    {
        box_size[d] = 0;
        for (i = 0; i < DIM; i++)
        {
            box_size[d] += box[d][i] * box[d][i];
        }
        box_size[d] = std::sqrt(box_size[d]);
    }

    n[XX] = *nx;
    n[YY] = *ny;
    n[ZZ] = *nz;

    if ((*nx <= 0) || (*ny <= 0) || (*nz <= 0))
    {
        if (nullptr != fp)
        {
            fprintf(fp,
                    "Calculating fourier grid dimensions for%s%s%s\n",
                    *nx > 0 ? "" : " X",
                    *ny > 0 ? "" : " Y",
                    *nz > 0 ? "" : " Z");
        }
    }

    max_spacing = 0;
    for (d = 0; d < DIM; d++)
    {
        if (n[d] <= 0)
        {
            nmin = static_cast<int>(box_size[d] / gridSpacing + 0.999);
            nmin = std::max(nmin, minGridPointsPerDim);

            i = g_initNR - 1;
            if (grid_init[i] >= nmin)
            {
                /* Take the smallest possible grid in the list */
                while (i > 0 && grid_init[i - 1] >= nmin)
                {
                    i--;
                }
                n[d] = grid_init[i];
            }
            else
            {
                /* Determine how many pre-factors of 2 we need */
                fac2 = 1;
                i    = g_baseNR - 1;
                while (fac2 * grid_base[i] < nmin)
                {
                    fac2 *= 2;
                }
                /* Find the smallest grid that is >= nmin */
                do
                {
                    attempt = fac2 * grid_base[i];
                    /* We demand a factor of 4, avoid 140, allow 90 */
                    if (((attempt % 4 == 0 && attempt != 140) || attempt == 90) && attempt >= nmin)
                    {
                        n[d] = attempt;
                    }
                    i--;
                } while (i > 0);
            }
        }

        spacing[d]  = box_size[d] / n[d];
        max_spacing = std::max(max_spacing, spacing[d]);
    }
    *nx = n[XX];
    *ny = n[YY];
    *nz = n[ZZ];
    if (nullptr != fp)
    {
        fprintf(fp,
                "Using a fourier grid of %dx%dx%d, spacing %.3f %.3f %.3f\n",
                *nx,
                *ny,
                *nz,
                spacing[XX],
                spacing[YY],
                spacing[ZZ]);
    }

    return max_spacing;
}
