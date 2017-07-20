/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/smalloc.h"

#include "pme-internal.h"
#include "pme-simd.h"
#include "pme-spline-work.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

template <typename INT>
static inline RVec
do_fspline(INT order, const struct gmx_pme_t *pme, const real *grid,
           const pme_atomcomm_t *atc, const splinedata_t *spline, int nn)
{
    const int  pny   = pme->pmegrid_ny;
    const int  pnz   = pme->pmegrid_nz;

    const int  n      = spline->ind[nn];
    const int* idxptr = atc->idx[n];
    const int  norder = nn*order;

    const int  i0   = idxptr[XX];
    const int  j0   = idxptr[YY];
    const int  k0   = idxptr[ZZ];

    /* Pointer arithmetic alert, next six statements */
    const real* thx  = spline->theta[XX] + norder;
    const real* thy  = spline->theta[YY] + norder;
    const real* thz  = spline->theta[ZZ] + norder;
    const real* dthx = spline->dtheta[XX] + norder;
    const real* dthy = spline->dtheta[YY] + norder;
    const real* dthz = spline->dtheta[ZZ] + norder;

    RVec        f(0, 0, 0);

    for (int ithx = 0; (ithx < order); ithx++)
    {
        const int  index_x = (i0+ithx)*pny*pnz;
        const real tx      = thx[ithx];
        const real dx      = dthx[ithx];

        for (int ithy = 0; (ithy < order); ithy++)
        {
            const int  index_xy = index_x+(j0+ithy)*pnz;
            const real ty       = thy[ithy];
            const real dy       = dthy[ithy];
            real       fxy1     = 0, fz1 = 0;

            for (int ithz = 0; (ithz < order); ithz++)
            {
                const real gval  = grid[index_xy+(k0+ithz)];
                fxy1 += thz[ithz]*gval;
                fz1  += dthz[ithz]*gval;
            }
            f[XX] += dx*ty*fxy1;
            f[YY] += tx*dy*fxy1;
            f[ZZ] += tx*ty*fz1;
        }
    }

    return f;
}

#ifdef PME_SIMD4_UNALIGNED //TODO: Consider always have at least a dummy implementation of Simd (enough for first phase of two-phase lookup) and then use enable_if instead of #ifdef
/* Gather for one charge with pme_order=4 with unaligned SIMD4 load+store.
 * This code does not assume any memory alignment for the grid.
 */
static inline RVec
do_fspline(std::integral_constant<int, 4>, const struct gmx_pme_t *pme, const real *grid,
           const pme_atomcomm_t *atc, const splinedata_t *spline, int nn)
{
    Simd4Real  fx_S, fy_S, fz_S;

    Simd4Real  tx_S, ty_S, tz_S;
    Simd4Real  dx_S, dy_S, dz_S;

    Simd4Real  gval_S;

    Simd4Real  fxy1_S;
    Simd4Real  fz1_S;

    const int  pny   = pme->pmegrid_ny;
    const int  pnz   = pme->pmegrid_nz;

    const int  n      = spline->ind[nn];
    const int* idxptr = atc->idx[n];
    const int  norder = nn*4;

    const int  i0   = idxptr[XX];
    const int  j0   = idxptr[YY];
    const int  k0   = idxptr[ZZ];

    /* Pointer arithmetic alert, next six statements */
    const real* thx  = spline->theta[XX] + norder;
    const real* thy  = spline->theta[YY] + norder;
    const real* thz  = spline->theta[ZZ] + norder;
    const real* dthx = spline->dtheta[XX] + norder;
    const real* dthy = spline->dtheta[YY] + norder;
    const real* dthz = spline->dtheta[ZZ] + norder;

    fx_S = setZero();
    fy_S = setZero();
    fz_S = setZero();

    /* With order 4 the z-spline is actually aligned */
    tz_S  = load4(thz);
    dz_S  = load4(dthz);

    for (int ithx = 0; (ithx < 4); ithx++)
    {
        const int index_x  = (i0+ithx)*pny*pnz;
        tx_S     = Simd4Real(thx[ithx]);
        dx_S     = Simd4Real(dthx[ithx]);

        for (int ithy = 0; (ithy < 4); ithy++)
        {
            const int index_xy = index_x+(j0+ithy)*pnz;
            ty_S     = Simd4Real(thy[ithy]);
            dy_S     = Simd4Real(dthy[ithy]);

            gval_S = load4U(grid+index_xy+k0);

            fxy1_S = tz_S * gval_S;
            fz1_S  = dz_S * gval_S;

            fx_S = fma(dx_S * ty_S, fxy1_S, fx_S);
            fy_S = fma(tx_S * dy_S, fxy1_S, fy_S);
            fz_S = fma(tx_S * ty_S, fz1_S, fz_S);
        }
    }

    return {
               reduce(fx_S), reduce(fy_S), reduce(fz_S)
    };
}
#endif

#ifdef PME_SIMD4_SPREAD_GATHER
/* This code assumes that the grid is allocated 4-real aligned
 * and that pnz is a multiple of 4.
 * This code supports pme_order <= 5.
 */
template <int ORDER>
static inline typename std::enable_if<ORDER == 4 || ORDER == 5, RVec>::type
do_fspline(std::integral_constant<int, ORDER> order, const struct gmx_pme_t *pme, const real *grid,
           const pme_atomcomm_t *atc, const splinedata_t *spline, int nn)
{
    int              offset;

    Simd4Real        fx_S, fy_S, fz_S;

    Simd4Real        tx_S, ty_S, tz_S0, tz_S1;
    Simd4Real        dx_S, dy_S, dz_S0, dz_S1;

    Simd4Real        gval_S0;
    Simd4Real        gval_S1;

    Simd4Real        fxy1_S0;
    Simd4Real        fz1_S0;
    Simd4Real        fxy1_S1;
    Simd4Real        fz1_S1;
    Simd4Real        fxy1_S;
    Simd4Real        fz1_S;

    const int        pny   = pme->pmegrid_ny;
    const int        pnz   = pme->pmegrid_nz;

    const int        n      = spline->ind[nn];
    const int      * idxptr = atc->idx[n];
    const int        norder = nn*order;

    const int        i0   = idxptr[XX];
    const int        j0   = idxptr[YY];
    const int        k0   = idxptr[ZZ];

    /* Pointer arithmetic alert, next six statements */
    const real            * thx  = spline->theta[XX] + norder;
    const real            * thy  = spline->theta[YY] + norder;
    const real            * thz  = spline->theta[ZZ] + norder;
    const real            * dthx = spline->dtheta[XX] + norder;
    const real            * dthy = spline->dtheta[YY] + norder;
    const real            * dthz = spline->dtheta[ZZ] + norder;

    struct pme_spline_work *work = pme->spline_work;

    offset = k0 & 3;

    fx_S = setZero();
    fy_S = setZero();
    fz_S = setZero();

#ifdef PME_SIMD4_UNALIGNED
    tz_S0 = load4U(thz-offset);
    tz_S1 = load4U(thz-offset+4);
    dz_S0 = load4U(dthz-offset);
    dz_S1 = load4U(dthz-offset+4);
#else
    GMX_ALIGNED(real, GMX_SIMD4_WIDTH)  thz_aligned[GMX_SIMD4_WIDTH*2];
    GMX_ALIGNED(real, GMX_SIMD4_WIDTH)  dthz_aligned[GMX_SIMD4_WIDTH*2];
    {
        int i;
        /* Copy (d)thz to an aligned buffer (unused buffer parts are masked) */
        for (i = 0; i < order; i++)
        {
            thz_aligned[offset+i]  = thz[i];
            dthz_aligned[offset+i] = dthz[i];
        }
        tz_S0 = load4(thz_aligned);
        tz_S1 = load4(thz_aligned+4);
        dz_S0 = load4(dthz_aligned);
        dz_S1 = load4(dthz_aligned+4);
    }
#endif
    tz_S0 = selectByMask(tz_S0, work->mask_S0[offset]);
    dz_S0 = selectByMask(dz_S0, work->mask_S0[offset]);
    tz_S1 = selectByMask(tz_S1, work->mask_S1[offset]);
    dz_S1 = selectByMask(dz_S1, work->mask_S1[offset]);

    for (int ithx = 0; (ithx < order); ithx++)
    {
        const int index_x  = (i0+ithx)*pny*pnz;
        tx_S     = Simd4Real(thx[ithx]);
        dx_S     = Simd4Real(dthx[ithx]);

        for (int ithy = 0; (ithy < order); ithy++)
        {
            const int index_xy = index_x+(j0+ithy)*pnz;
            ty_S     = Simd4Real(thy[ithy]);
            dy_S     = Simd4Real(dthy[ithy]);

            gval_S0 = load4(grid+index_xy+k0-offset);
            gval_S1 = load4(grid+index_xy+k0-offset+4);

            fxy1_S0 = tz_S0 * gval_S0;
            fz1_S0  = dz_S0 * gval_S0;
            fxy1_S1 = tz_S1 * gval_S1;
            fz1_S1  = dz_S1 * gval_S1;

            fxy1_S = fxy1_S0 + fxy1_S1;
            fz1_S  = fz1_S0 + fz1_S1;

            fx_S = fma(dx_S * ty_S, fxy1_S, fx_S);
            fy_S = fma(tx_S * dy_S, fxy1_S, fy_S);
            fz_S = fma(tx_S * ty_S, fz1_S, fz_S);
        }
    }

    return {
               reduce(fx_S), reduce(fy_S), reduce(fz_S)
    };
}
#endif


void gather_f_bsplines(const struct gmx_pme_t *pme, const real *grid,
                       gmx_bool bClearF, const pme_atomcomm_t *atc,
                       const splinedata_t *spline,
                       real scale)
{
    const int  nx    = pme->nkx;
    const int  ny    = pme->nky;
    const int  nz    = pme->nkz;

    const real rxx   = pme->recipbox[XX][XX];
    const real ryx   = pme->recipbox[YY][XX];
    const real ryy   = pme->recipbox[YY][YY];
    const real rzx   = pme->recipbox[ZZ][XX];
    const real rzy   = pme->recipbox[ZZ][YY];
    const real rzz   = pme->recipbox[ZZ][ZZ];

    const int  order  = pme->pme_order;

    /* sum forces for local particles */
    for (int nn = 0; nn < spline->n; nn++)
    {
        const int  n           = spline->ind[nn];
        const real coefficient = scale*atc->coefficient[n];

        if (bClearF)
        {
            atc->f[n][XX] = 0;
            atc->f[n][YY] = 0;
            atc->f[n][ZZ] = 0;
        }
        if (coefficient != 0)
        {
            RVec f;

            switch (order)
            {
                case 4:
                    f = do_fspline(std::integral_constant<int, 4>(), //TODO: When C++14 is allowed use generic lambda to avoid passing arguments 3x (and potentially replacing switch statement with a dispatch function)
                                   pme, grid, atc, spline, nn);
                    break;
                case 5:
                    f = do_fspline(std::integral_constant<int, 5>(),
                                   pme, grid, atc, spline, nn);
                    break;
                default:
                    f = do_fspline(order,
                                   pme, grid, atc, spline, nn);
                    break;
            }

            atc->f[n][XX] += -coefficient*( f[XX]*nx*rxx );
            atc->f[n][YY] += -coefficient*( f[XX]*nx*ryx + f[YY]*ny*ryy );
            atc->f[n][ZZ] += -coefficient*( f[XX]*nx*rzx + f[YY]*ny*rzy + f[ZZ]*nz*rzz );
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
    int           n, ithx, ithy, ithz, i0, j0, k0;
    int           index_x, index_xy;
    int          *idxptr;
    real          energy, pot, tx, ty, coefficient, gval;
    real         *thx, *thy, *thz;
    int           norder;
    int           order;

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
