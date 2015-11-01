/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

/* This include file has code between ifdef's to make sure
 * that this performance sensitive code is inlined
 * and to remove conditionals and variable loop bounds at compile time.
 */

#ifdef PME_SPREAD_SIMD4_ORDER4
/* Spread one charge with pme_order=4 with unaligned SIMD4 load+store.
 * This code does not assume any memory alignment for the grid.
 */
{
    gmx::Simd4Real ty_S0, ty_S1, ty_S2, ty_S3;
    gmx::Simd4Real tz_S;
    gmx::Simd4Real vx_S;
    gmx::Simd4Real vx_tz_S;
    gmx::Simd4Real sum_S0, sum_S1, sum_S2, sum_S3;
    gmx::Simd4Real gri_S0, gri_S1, gri_S2, gri_S3;

    ty_S0 = gmx::simd4Set1(thy[0]);
    ty_S1 = gmx::simd4Set1(thy[1]);
    ty_S2 = gmx::simd4Set1(thy[2]);
    ty_S3 = gmx::simd4Set1(thy[3]);

    /* With order 4 the z-spline is actually aligned */
    tz_S  = gmx::simd4Load(thz);

    for (ithx = 0; (ithx < 4); ithx++)
    {
        index_x = (i0+ithx)*pny*pnz;
        valx    = coefficient*thx[ithx];

        vx_S   = gmx::simd4Set1(valx);

        vx_tz_S = gmx::simd4Mul(vx_S, tz_S);

        gri_S0 = gmx::simd4LoadU(grid+index_x+(j0+0)*pnz+k0);
        gri_S1 = gmx::simd4LoadU(grid+index_x+(j0+1)*pnz+k0);
        gri_S2 = gmx::simd4LoadU(grid+index_x+(j0+2)*pnz+k0);
        gri_S3 = gmx::simd4LoadU(grid+index_x+(j0+3)*pnz+k0);

        sum_S0 = gmx::simd4Fmadd(vx_tz_S, ty_S0, gri_S0);
        sum_S1 = gmx::simd4Fmadd(vx_tz_S, ty_S1, gri_S1);
        sum_S2 = gmx::simd4Fmadd(vx_tz_S, ty_S2, gri_S2);
        sum_S3 = gmx::simd4Fmadd(vx_tz_S, ty_S3, gri_S3);

        gmx::simd4StoreU(grid+index_x+(j0+0)*pnz+k0, sum_S0);
        gmx::simd4StoreU(grid+index_x+(j0+1)*pnz+k0, sum_S1);
        gmx::simd4StoreU(grid+index_x+(j0+2)*pnz+k0, sum_S2);
        gmx::simd4StoreU(grid+index_x+(j0+3)*pnz+k0, sum_S3);
    }
}
#undef PME_SPREAD_SIMD4_ORDER4
#endif


#ifdef PME_GATHER_F_SIMD4_ORDER4
/* Gather for one charge with pme_order=4 with unaligned SIMD4 load+store.
 * This code does not assume any memory alignment for the grid.
 */
{
    gmx::Simd4Real fx_S, fy_S, fz_S;

    gmx::Simd4Real tx_S, ty_S, tz_S;
    gmx::Simd4Real dx_S, dy_S, dz_S;

    gmx::Simd4Real gval_S;

    gmx::Simd4Real fxy1_S;
    gmx::Simd4Real fz1_S;

    fx_S = gmx::simd4SetZero();
    fy_S = gmx::simd4SetZero();
    fz_S = gmx::simd4SetZero();

    /* With order 4 the z-spline is actually aligned */
    tz_S  = gmx::simd4Load(thz);
    dz_S  = gmx::simd4Load(dthz);

    for (ithx = 0; (ithx < 4); ithx++)
    {
        index_x  = (i0+ithx)*pny*pnz;
        tx_S     = gmx::simd4Set1(thx[ithx]);
        dx_S     = gmx::simd4Set1(dthx[ithx]);

        for (ithy = 0; (ithy < 4); ithy++)
        {
            index_xy = index_x+(j0+ithy)*pnz;
            ty_S     = gmx::simd4Set1(thy[ithy]);
            dy_S     = gmx::simd4Set1(dthy[ithy]);

            gval_S = gmx::simd4LoadU(grid+index_xy+k0);

            fxy1_S = gmx::simd4Mul(tz_S, gval_S);
            fz1_S  = gmx::simd4Mul(dz_S, gval_S);

            fx_S = gmx::simd4Fmadd(gmx::simd4Mul(dx_S, ty_S), fxy1_S, fx_S);
            fy_S = gmx::simd4Fmadd(gmx::simd4Mul(tx_S, dy_S), fxy1_S, fy_S);
            fz_S = gmx::simd4Fmadd(gmx::simd4Mul(tx_S, ty_S), fz1_S, fz_S);
        }
    }

    fx += gmx::simd4Reduce(fx_S);
    fy += gmx::simd4Reduce(fy_S);
    fz += gmx::simd4Reduce(fz_S);
}
#undef PME_GATHER_F_SIMD4_ORDER4
#endif


#ifdef PME_SPREAD_SIMD4_ALIGNED
/* This code assumes that the grid is allocated 4-real aligned
 * and that pnz is a multiple of 4.
 * This code supports pme_order <= 5.
 */
{
    int              offset;
    int              index;
    gmx::Simd4Real   ty_S0, ty_S1, ty_S2, ty_S3, ty_S4;
    gmx::Simd4Real   tz_S0;
    gmx::Simd4Real   tz_S1;
    gmx::Simd4Real   vx_S;
    gmx::Simd4Real   vx_tz_S0;
    gmx::Simd4Real   vx_tz_S1;
    gmx::Simd4Real   sum_S00, sum_S01, sum_S02, sum_S03, sum_S04;
    gmx::Simd4Real   sum_S10, sum_S11, sum_S12, sum_S13, sum_S14;
    gmx::Simd4Real   gri_S00, gri_S01, gri_S02, gri_S03, gri_S04;
    gmx::Simd4Real   gri_S10, gri_S11, gri_S12, gri_S13, gri_S14;

    offset = k0 & 3;

    ty_S0 = gmx::simd4Set1(thy[0]);
    ty_S1 = gmx::simd4Set1(thy[1]);
    ty_S2 = gmx::simd4Set1(thy[2]);
    ty_S3 = gmx::simd4Set1(thy[3]);
#if PME_ORDER == 5
    ty_S4 = gmx::simd4Set1(thy[4]);
#endif

#ifdef PME_SIMD4_UNALIGNED
    tz_S0 = gmx::simd4LoadU(thz-offset);
    tz_S1 = gmx::simd4LoadU(thz-offset+4);
#else
    {
        int i;
        /* Copy thz to an aligned buffer (unused buffer parts are masked) */
        for (i = 0; i < PME_ORDER; i++)
        {
            thz_aligned[offset+i] = thz[i];
        }
        tz_S0 = gmx::simd4Load(thz_aligned);
        tz_S1 = gmx::simd4Load(thz_aligned+4);
    }
#endif
    tz_S0 = gmx::simd4Mask(tz_S0, work->mask_S0[offset]);
    tz_S1 = gmx::simd4Mask(tz_S1, work->mask_S1[offset]);

    for (ithx = 0; (ithx < PME_ORDER); ithx++)
    {
        index = (i0+ithx)*pny*pnz + j0*pnz + k0 - offset;
        valx  = coefficient*thx[ithx];

        vx_S   = gmx::simd4Set1(valx);

        vx_tz_S0 = gmx::simd4Mul(vx_S, tz_S0);
        vx_tz_S1 = gmx::simd4Mul(vx_S, tz_S1);

        gri_S00 = gmx::simd4Load(grid+index+0*pnz);
        gri_S01 = gmx::simd4Load(grid+index+1*pnz);
        gri_S02 = gmx::simd4Load(grid+index+2*pnz);
        gri_S03 = gmx::simd4Load(grid+index+3*pnz);
#if PME_ORDER == 5
        gri_S04 = gmx::simd4Load(grid+index+4*pnz);
#endif
        gri_S10 = gmx::simd4Load(grid+index+0*pnz+4);
        gri_S11 = gmx::simd4Load(grid+index+1*pnz+4);
        gri_S12 = gmx::simd4Load(grid+index+2*pnz+4);
        gri_S13 = gmx::simd4Load(grid+index+3*pnz+4);
#if PME_ORDER == 5
        gri_S14 = gmx::simd4Load(grid+index+4*pnz+4);
#endif

        sum_S00 = gmx::simd4Fmadd(vx_tz_S0, ty_S0, gri_S00);
        sum_S01 = gmx::simd4Fmadd(vx_tz_S0, ty_S1, gri_S01);
        sum_S02 = gmx::simd4Fmadd(vx_tz_S0, ty_S2, gri_S02);
        sum_S03 = gmx::simd4Fmadd(vx_tz_S0, ty_S3, gri_S03);
#if PME_ORDER == 5
        sum_S04 = gmx::simd4Fmadd(vx_tz_S0, ty_S4, gri_S04);
#endif
        sum_S10 = gmx::simd4Fmadd(vx_tz_S1, ty_S0, gri_S10);
        sum_S11 = gmx::simd4Fmadd(vx_tz_S1, ty_S1, gri_S11);
        sum_S12 = gmx::simd4Fmadd(vx_tz_S1, ty_S2, gri_S12);
        sum_S13 = gmx::simd4Fmadd(vx_tz_S1, ty_S3, gri_S13);
#if PME_ORDER == 5
        sum_S14 = gmx::simd4Fmadd(vx_tz_S1, ty_S4, gri_S14);
#endif

        gmx::simd4Store(grid+index+0*pnz, sum_S00);
        gmx::simd4Store(grid+index+1*pnz, sum_S01);
        gmx::simd4Store(grid+index+2*pnz, sum_S02);
        gmx::simd4Store(grid+index+3*pnz, sum_S03);
#if PME_ORDER == 5
        gmx::simd4Store(grid+index+4*pnz, sum_S04);
#endif
        gmx::simd4Store(grid+index+0*pnz+4, sum_S10);
        gmx::simd4Store(grid+index+1*pnz+4, sum_S11);
        gmx::simd4Store(grid+index+2*pnz+4, sum_S12);
        gmx::simd4Store(grid+index+3*pnz+4, sum_S13);
#if PME_ORDER == 5
        gmx::simd4Store(grid+index+4*pnz+4, sum_S14);
#endif
    }
}
#undef PME_ORDER
#undef PME_SPREAD_SIMD4_ALIGNED
#endif


#ifdef PME_GATHER_F_SIMD4_ALIGNED
/* This code assumes that the grid is allocated 4-real aligned
 * and that pnz is a multiple of 4.
 * This code supports pme_order <= 5.
 */
{
    int              offset;

    gmx::Simd4Real   fx_S, fy_S, fz_S;

    gmx::Simd4Real   tx_S, ty_S, tz_S0, tz_S1;
    gmx::Simd4Real   dx_S, dy_S, dz_S0, dz_S1;

    gmx::Simd4Real   gval_S0;
    gmx::Simd4Real   gval_S1;

    gmx::Simd4Real   fxy1_S0;
    gmx::Simd4Real   fz1_S0;
    gmx::Simd4Real   fxy1_S1;
    gmx::Simd4Real   fz1_S1;
    gmx::Simd4Real   fxy1_S;
    gmx::Simd4Real   fz1_S;

    offset = k0 & 3;

    fx_S = gmx::simd4SetZero();
    fy_S = gmx::simd4SetZero();
    fz_S = gmx::simd4SetZero();

#ifdef PME_SIMD4_UNALIGNED
    tz_S0 = gmx::simd4LoadU(thz-offset);
    tz_S1 = gmx::simd4LoadU(thz-offset+4);
    dz_S0 = gmx::simd4LoadU(dthz-offset);
    dz_S1 = gmx::simd4LoadU(dthz-offset+4);
#else
    {
        int i;
        /* Copy (d)thz to an aligned buffer (unused buffer parts are masked) */
        for (i = 0; i < PME_ORDER; i++)
        {
            thz_aligned[offset+i]  = thz[i];
            dthz_aligned[offset+i] = dthz[i];
        }
        tz_S0 = gmx::simd4Load(thz_aligned);
        tz_S1 = gmx::simd4Load(thz_aligned+4);
        dz_S0 = gmx::simd4Load(dthz_aligned);
        dz_S1 = gmx::simd4Load(dthz_aligned+4);
    }
#endif
    tz_S0 = gmx::simd4Mask(tz_S0, work->mask_S0[offset]);
    dz_S0 = gmx::simd4Mask(dz_S0, work->mask_S0[offset]);
    tz_S1 = gmx::simd4Mask(tz_S1, work->mask_S1[offset]);
    dz_S1 = gmx::simd4Mask(dz_S1, work->mask_S1[offset]);

    for (ithx = 0; (ithx < PME_ORDER); ithx++)
    {
        index_x  = (i0+ithx)*pny*pnz;
        tx_S     = gmx::simd4Set1(thx[ithx]);
        dx_S     = gmx::simd4Set1(dthx[ithx]);

        for (ithy = 0; (ithy < PME_ORDER); ithy++)
        {
            index_xy = index_x+(j0+ithy)*pnz;
            ty_S     = gmx::simd4Set1(thy[ithy]);
            dy_S     = gmx::simd4Set1(dthy[ithy]);

            gval_S0 = gmx::simd4Load(grid+index_xy+k0-offset);
            gval_S1 = gmx::simd4Load(grid+index_xy+k0-offset+4);

            fxy1_S0 = gmx::simd4Mul(tz_S0, gval_S0);
            fz1_S0  = gmx::simd4Mul(dz_S0, gval_S0);
            fxy1_S1 = gmx::simd4Mul(tz_S1, gval_S1);
            fz1_S1  = gmx::simd4Mul(dz_S1, gval_S1);

            fxy1_S = gmx::simd4Add(fxy1_S0, fxy1_S1);
            fz1_S  = gmx::simd4Add(fz1_S0, fz1_S1);

            fx_S = gmx::simd4Fmadd(gmx::simd4Mul(dx_S, ty_S), fxy1_S, fx_S);
            fy_S = gmx::simd4Fmadd(gmx::simd4Mul(tx_S, dy_S), fxy1_S, fy_S);
            fz_S = gmx::simd4Fmadd(gmx::simd4Mul(tx_S, ty_S), fz1_S, fz_S);
        }
    }

    fx += gmx::simd4Reduce(fx_S);
    fy += gmx::simd4Reduce(fy_S);
    fz += gmx::simd4Reduce(fz_S);
}
#undef PME_ORDER
#undef PME_GATHER_F_SIMD4_ALIGNED
#endif
