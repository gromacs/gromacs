/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "pme-gather.h"

#include "config.h"

#include "gromacs/ewald/pme-internal.h"
#include "gromacs/ewald/pme-simd.h"
#include "gromacs/ewald/pme-spline-work.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"

static inline void
do_fspline(int order,
           int pny, int pnz, int i0, int j0, int k0,
           const real *thx, const real *thy, const real *thz,
           const real *dthx, const real *dthy, const real *dthz,
           const real *grid,
           real &fx, real &fy, real &fz)
{
    for (int ithx = 0; (ithx < order); ithx++)
    {
        int  index_x = (i0+ithx)*pny*pnz;
        real tx      = thx[ithx];
        real dx      = dthx[ithx];

        for (int ithy = 0; (ithy < order); ithy++)
        {
            int  index_xy = index_x+(j0+ithy)*pnz;
            real ty       = thy[ithy];
            real dy       = dthy[ithy];
            real fxy1     = 0, fz1 = 0;

            for (int ithz = 0; (ithz < order); ithz++)
            {
                real gval  = grid[index_xy+(k0+ithz)];
                fxy1 += thz[ithz]*gval;
                fz1  += dthz[ithz]*gval;
            }
            fx += dx*ty*fxy1;
            fy += tx*dy*fxy1;
            fz += tx*ty*fz1;
        }
    }
}

template <int order>
static inline void
do_fspline_order(int pny, int pnz, int i0, int j0, int k0,
                 const real *thx, const real *thy, const real *thz,
                 const real *dthx, const real *dthy, const real *dthz,
                 const real *grid,
                 const struct pme_spline_work gmx_unused *work,
                 real gmx_unused *thz_aligned, real gmx_unused *dthz_aligned,
                 real &fx, real &fy, real &fz)
{
    do_fspline(order, pny, pnz, i0, j0, k0,
               thx, thy, thz, dthx, dthy, dthz,
               grid, fx, fy, fz);
}

#ifdef PME_SIMD4_SPREAD_GATHER
/* This code assumes that the grid is allocated 4-real aligned
 * and that pnz is a multiple of 4.
 * This code supports pme_order <= 5.
 */
static inline void
do_fspline_simd(int order, int pny, int pnz, int i0, int j0, int k0,
                const real *thx, const real *thy, const real *thz,
                const real *dthx, const real *dthy, const real *dthz,
                const real *grid, const struct pme_spline_work *work,
                real gmx_unused *thz_aligned, real gmx_unused *dthz_aligned,
                real &fx, real &fy, real &fz)
{
    int              offset;

    gmx_simd4_real_t fx_S, fy_S, fz_S;

    gmx_simd4_real_t tx_S, ty_S, tz_S0, tz_S1;
    gmx_simd4_real_t dx_S, dy_S, dz_S0, dz_S1;

    gmx_simd4_real_t gval_S0;
    gmx_simd4_real_t gval_S1;

    gmx_simd4_real_t fxy1_S0;
    gmx_simd4_real_t fz1_S0;
    gmx_simd4_real_t fxy1_S1;
    gmx_simd4_real_t fz1_S1;
    gmx_simd4_real_t fxy1_S;
    gmx_simd4_real_t fz1_S;

    offset = k0 & 3;

    fx_S = gmx_simd4_setzero_r();
    fy_S = gmx_simd4_setzero_r();
    fz_S = gmx_simd4_setzero_r();

#ifdef PME_SIMD4_UNALIGNED
    tz_S0 = gmx_simd4_loadu_r(thz-offset);
    tz_S1 = gmx_simd4_loadu_r(thz-offset+4);
    dz_S0 = gmx_simd4_loadu_r(dthz-offset);
    dz_S1 = gmx_simd4_loadu_r(dthz-offset+4);
#else
    /* Copy (d)thz to an aligned buffer (unused buffer parts are masked) */
    for (int i = 0; i < order; i++)
    {
        thz_aligned[offset+i]  = thz[i];
        dthz_aligned[offset+i] = dthz[i];
    }
    tz_S0 = gmx_simd4_load_r(thz_aligned);
    tz_S1 = gmx_simd4_load_r(thz_aligned+4);
    dz_S0 = gmx_simd4_load_r(dthz_aligned);
    dz_S1 = gmx_simd4_load_r(dthz_aligned+4);
#endif
    tz_S0 = gmx_simd4_blendzero_r(tz_S0, work->mask_S0[offset]);
    dz_S0 = gmx_simd4_blendzero_r(dz_S0, work->mask_S0[offset]);
    tz_S1 = gmx_simd4_blendzero_r(tz_S1, work->mask_S1[offset]);
    dz_S1 = gmx_simd4_blendzero_r(dz_S1, work->mask_S1[offset]);

    for (int ithx = 0; (ithx < order); ithx++)
    {
        int index_x  = (i0+ithx)*pny*pnz;
        tx_S     = gmx_simd4_set1_r(thx[ithx]);
        dx_S     = gmx_simd4_set1_r(dthx[ithx]);

        for (int ithy = 0; (ithy < order); ithy++)
        {
            int index_xy = index_x+(j0+ithy)*pnz;
            ty_S     = gmx_simd4_set1_r(thy[ithy]);
            dy_S     = gmx_simd4_set1_r(dthy[ithy]);

            gval_S0 = gmx_simd4_load_r(grid+index_xy+k0-offset);
            gval_S1 = gmx_simd4_load_r(grid+index_xy+k0-offset+4);

            fxy1_S0 = gmx_simd4_mul_r(tz_S0, gval_S0);
            fz1_S0  = gmx_simd4_mul_r(dz_S0, gval_S0);
            fxy1_S1 = gmx_simd4_mul_r(tz_S1, gval_S1);
            fz1_S1  = gmx_simd4_mul_r(dz_S1, gval_S1);

            fxy1_S = gmx_simd4_add_r(fxy1_S0, fxy1_S1);
            fz1_S  = gmx_simd4_add_r(fz1_S0, fz1_S1);

            fx_S = gmx_simd4_fmadd_r(gmx_simd4_mul_r(dx_S, ty_S), fxy1_S, fx_S);
            fy_S = gmx_simd4_fmadd_r(gmx_simd4_mul_r(tx_S, dy_S), fxy1_S, fy_S);
            fz_S = gmx_simd4_fmadd_r(gmx_simd4_mul_r(tx_S, ty_S), fz1_S, fz_S);
        }
    }

    fx += gmx_simd4_reduce_r(fx_S);
    fy += gmx_simd4_reduce_r(fy_S);
    fz += gmx_simd4_reduce_r(fz_S);
}
#ifdef PME_SIMD4_UNALIGNED
/* Gather for one charge with pme_order=4 with unaligned SIMD4 load+store.
 * This code does not assume any memory alignment for the grid.
 */
static inline void
do_fspline_simd_unaligned4(int pny, int pnz, int i0, int j0, int k0,
                           const real *thx, const real *thy, const real *thz,
                           const real *dthx, const real *dthy, const real *dthz,
                           const real *grid,
                           real &fx, real &fy, real &fz)
{
    gmx_simd4_real_t fx_S, fy_S, fz_S;

    gmx_simd4_real_t tx_S, ty_S, tz_S;
    gmx_simd4_real_t dx_S, dy_S, dz_S;

    gmx_simd4_real_t gval_S;

    gmx_simd4_real_t fxy1_S;
    gmx_simd4_real_t fz1_S;

    fx_S = gmx_simd4_setzero_r();
    fy_S = gmx_simd4_setzero_r();
    fz_S = gmx_simd4_setzero_r();

    /* With order 4 the z-spline is actually aligned */
    tz_S  = gmx_simd4_load_r(thz);
    dz_S  = gmx_simd4_load_r(dthz);

    for (int ithx = 0; (ithx < 4); ithx++)
    {
        int index_x  = (i0+ithx)*pny*pnz;
        tx_S     = gmx_simd4_set1_r(thx[ithx]);
        dx_S     = gmx_simd4_set1_r(dthx[ithx]);

        for (int ithy = 0; (ithy < 4); ithy++)
        {
            int index_xy = index_x+(j0+ithy)*pnz;
            ty_S     = gmx_simd4_set1_r(thy[ithy]);
            dy_S     = gmx_simd4_set1_r(dthy[ithy]);

            gval_S = gmx_simd4_loadu_r(grid+index_xy+k0);

            fxy1_S = gmx_simd4_mul_r(tz_S, gval_S);
            fz1_S  = gmx_simd4_mul_r(dz_S, gval_S);

            fx_S = gmx_simd4_fmadd_r(gmx_simd4_mul_r(dx_S, ty_S), fxy1_S, fx_S);
            fy_S = gmx_simd4_fmadd_r(gmx_simd4_mul_r(tx_S, dy_S), fxy1_S, fy_S);
            fz_S = gmx_simd4_fmadd_r(gmx_simd4_mul_r(tx_S, ty_S), fz1_S, fz_S);
        }
    }

    fx += gmx_simd4_reduce_r(fx_S);
    fy += gmx_simd4_reduce_r(fy_S);
    fz += gmx_simd4_reduce_r(fz_S);
}
#endif
template <>
inline void
do_fspline_order<4>(int pny, int pnz, int i0, int j0, int k0,
                    const real *thx, const real *thy, const real *thz,
                    const real *dthx, const real *dthy, const real *dthz,
                    const real *grid,
                    const struct pme_spline_work gmx_unused *work,
                    real gmx_unused *thz_aligned, real gmx_unused *dthz_aligned,
                    real &fx, real &fy, real &fz)
{
#ifdef PME_SIMD4_UNALIGNED
    do_fspline_simd_unaligned4(pny, pnz, i0, j0, k0,
                               thx, thy, thz, dthx, dthy, dthz,
                               grid, fx, fy, fz);
#else
    do_fspline_simd(4, pny, pnz, i0, j0, k0,
                    thx, thy, thz, dthx, dthy, dthz,
                    grid, work, thz_aligned, dthz_aligned, fx, fy, fz);
#endif
}

template <>
inline void
do_fspline_order<5>(int pny, int pnz, int i0, int j0, int k0,
                    const real *thx, const real *thy, const real *thz,
                    const real *dthx, const real *dthy, const real *dthz,
                    const real *grid,
                    const struct pme_spline_work gmx_unused *work,
                    real gmx_unused *thz_aligned, real gmx_unused *dthz_aligned,
                    real &fx, real &fy, real &fz)
{
    do_fspline_simd(5, pny, pnz, i0, j0, k0,
                    thx, thy, thz, dthx, dthy, dthz,
                    grid, work, thz_aligned, dthz_aligned, fx, fy, fz);
}
#endif

void gather_f_bsplines(struct gmx_pme *pme, real *grid,
                       gmx_bool bClearF, pme_atomcomm_t *atc,
                       splinedata_t *spline,
                       real scale)
{
    /* sum forces for local particles */
    int                     nn, n, i0, j0, k0;
    int                     nx, ny, nz, pny, pnz;
    int                 *   idxptr;
    real                    coefficient;
    real                    fx, fy, fz;
    real                   *thx, *thy, *thz, *dthx, *dthy, *dthz;
    int                     norder;
    real                    rxx, ryx, ryy, rzx, rzy, rzz;
    int                     order;

    struct pme_spline_work *work = pme->spline_work;
    real                    thz_buffer[GMX_SIMD4_WIDTH*3],  *thz_aligned;
    real                    dthz_buffer[GMX_SIMD4_WIDTH*3], *dthz_aligned;

#ifdef PME_SIMD4_SPREAD_GATHER
    thz_aligned  = gmx_simd4_align_r(thz_buffer);
    dthz_aligned = gmx_simd4_align_r(dthz_buffer);
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
                    do_fspline_order<4>(pny, pnz, i0, j0, k0,
                                        thx, thy, thz, dthx, dthy, dthz,
                                        grid, work, thz_aligned, dthz_aligned,
                                        fx, fy, fz);
                    break;
                case 5:
                    do_fspline_order<5>(pny, pnz, i0, j0, k0,
                                        thx, thy, thz, dthx, dthy, dthz,
                                        grid, work, thz_aligned, dthz_aligned,
                                        fx, fy, fz);
                    break;
                default:
                    do_fspline(order, pny, pnz, i0, j0, k0,
                               thx, thy, thz, dthx, dthy, dthz,
                               grid, fx, fy, fz);
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


real gather_energy_bsplines(struct gmx_pme *pme, real *grid,
                            pme_atomcomm_t *atc)
{
    splinedata_t *spline;
    int           n, ithx, ithy, ithz, i0, j0, k0;
    int           index_x, index_xy;
    int       *   idxptr;
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
