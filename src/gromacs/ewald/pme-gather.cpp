/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "gromacs/math/vec.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/smalloc.h"

#include "pme-internal.h"
#include "pme-simd.h"
#include "pme-spline-work.h"

#define DO_FSPLINE(order)                      \
    for (ithx = 0; (ithx < order); ithx++)              \
    {                                              \
        index_x = (i0+ithx)*pny*pnz;               \
        tx      = thx[ithx];                       \
        dx      = dthx[ithx];                      \
                                               \
        for (ithy = 0; (ithy < order); ithy++)          \
        {                                          \
            index_xy = index_x+(j0+ithy)*pnz;      \
            ty       = thy[ithy];                  \
            dy       = dthy[ithy];                 \
            fxy1     = fz1 = 0;                    \
                                               \
            for (ithz = 0; (ithz < order); ithz++)      \
            {                                      \
                gval  = grid[index_xy+(k0+ithz)];  \
                fxy1 += thz[ithz]*gval;            \
                fz1  += dthz[ithz]*gval;           \
            }                                      \
            fx += dx*ty*fxy1;                      \
            fy += tx*dy*fxy1;                      \
            fz += tx*ty*fz1;                       \
        }                                          \
    }


void gather_f_bsplines(struct gmx_pme_t *pme, real *grid,
                       gmx_bool bClearF, pme_atomcomm_t *atc,
                       splinedata_t *spline,
                       real scale)
{
    /* sum forces for local particles */
    int    nn, n, ithx, ithy, ithz, i0, j0, k0;
    int    index_x, index_xy;
    int    nx, ny, nz, pny, pnz;
    int   *idxptr;
    real   tx, ty, dx, dy, coefficient;
    real   fx, fy, fz, gval;
    real   fxy1, fz1;
    real  *thx, *thy, *thz, *dthx, *dthy, *dthz;
    int    norder;
    real   rxx, ryx, ryy, rzx, rzy, rzz;
    int    order;

#ifdef PME_SIMD4_SPREAD_GATHER
    // cppcheck-suppress unreadVariable cppcheck seems not to analyze code from pme-simd4.h
    struct pme_spline_work *work = pme->spline_work;
#ifndef PME_SIMD4_UNALIGNED
    real                    thz_buffer[GMX_SIMD4_WIDTH*3],  *thz_aligned;
    real                    dthz_buffer[GMX_SIMD4_WIDTH*3], *dthz_aligned;

    thz_aligned  = gmx_simd4_align_r(thz_buffer);
    dthz_aligned = gmx_simd4_align_r(dthz_buffer);
#endif
#endif

    order = pme->pme_order;
    nx    = pme->nkx;
    ny    = pme->nky;
    nz    = pme->nkz;
    pny   = pme->pmegrid_ny;
    pnz   = pme->pmegrid_nz;

    rxx   = pme->recipbox[XX][XX];
    ryx   = pme->recipbox[YY][XX];
    ryy   = pme->recipbox[YY][YY];
    rzx   = pme->recipbox[ZZ][XX];
    rzy   = pme->recipbox[ZZ][YY];
    rzz   = pme->recipbox[ZZ][ZZ];

    for (nn = 0; nn < spline->n; nn++)
    {
        n           = spline->ind[nn];
        coefficient = scale*atc->coefficient[n];

        if (bClearF)
        {
            atc->f[n][XX] = 0;
            atc->f[n][YY] = 0;
            atc->f[n][ZZ] = 0;
        }
        if (coefficient != 0)
        {
            fx     = 0;
            fy     = 0;
            fz     = 0;
            idxptr = atc->idx[n];
            norder = nn*order;

            i0   = idxptr[XX];
            j0   = idxptr[YY];
            k0   = idxptr[ZZ];

            /* Pointer arithmetic alert, next six statements */
            thx  = spline->theta[XX] + norder;
            thy  = spline->theta[YY] + norder;
            thz  = spline->theta[ZZ] + norder;
            dthx = spline->dtheta[XX] + norder;
            dthy = spline->dtheta[YY] + norder;
            dthz = spline->dtheta[ZZ] + norder;

            switch (order)
            {
                case 4:
#ifdef PME_SIMD4_SPREAD_GATHER
#ifdef PME_SIMD4_UNALIGNED
#define PME_GATHER_F_SIMD4_ORDER4
#else
#define PME_GATHER_F_SIMD4_ALIGNED
#define PME_ORDER 4
#endif
#include "pme-simd4.h"
#else
                    DO_FSPLINE(4);
#endif
                    break;
                case 5:
#ifdef PME_SIMD4_SPREAD_GATHER
#define PME_GATHER_F_SIMD4_ALIGNED
#define PME_ORDER 5
#include "pme-simd4.h"
#else
                    DO_FSPLINE(5);
#endif
                    break;
                default:
                    DO_FSPLINE(order);
                    break;
            }

            atc->f[n][XX] += -coefficient*( fx*nx*rxx );
            atc->f[n][YY] += -coefficient*( fx*nx*ryx + fy*ny*ryy );
            atc->f[n][ZZ] += -coefficient*( fx*nx*rzx + fy*ny*rzy + fz*nz*rzz );
        }
    }
    /* Since the energy and not forces are interpolated
     * the net force might not be exactly zero.
     * This can be solved by also interpolating F, but
     * that comes at a cost.
     * A better hack is to remove the net force every
     * step, but that must be done at a higher level
     * since this routine doesn't see all atoms if running
     * in parallel. Don't know how important it is?  EL 990726
     */
}


real gather_energy_bsplines(struct gmx_pme_t *pme, real *grid,
                            pme_atomcomm_t *atc)
{
    splinedata_t *spline;
    int     n, ithx, ithy, ithz, i0, j0, k0;
    int     index_x, index_xy;
    int *   idxptr;
    real    energy, pot, tx, ty, coefficient, gval;
    real    *thx, *thy, *thz;
    int     norder;
    int     order;

    spline = &atc->spline[0];

    order = pme->pme_order;

    energy = 0;
    for (n = 0; (n < atc->n); n++)
    {
        coefficient      = atc->coefficient[n];

        if (coefficient != 0)
        {
            idxptr = atc->idx[n];
            norder = n*order;

            i0   = idxptr[XX];
            j0   = idxptr[YY];
            k0   = idxptr[ZZ];

            /* Pointer arithmetic alert, next three statements */
            thx  = spline->theta[XX] + norder;
            thy  = spline->theta[YY] + norder;
            thz  = spline->theta[ZZ] + norder;

            pot = 0;
            for (ithx = 0; (ithx < order); ithx++)
            {
                index_x = (i0+ithx)*pme->pmegrid_ny*pme->pmegrid_nz;
                tx      = thx[ithx];

                for (ithy = 0; (ithy < order); ithy++)
                {
                    index_xy = index_x+(j0+ithy)*pme->pmegrid_nz;
                    ty       = thy[ithy];

                    for (ithz = 0; (ithz < order); ithz++)
                    {
                        gval  = grid[index_xy+(k0+ithz)];
                        pot  += tx*ty*thz[ithz]*gval;
                    }

                }
            }

            energy += pot*coefficient;
        }
    }

    return energy;
}
