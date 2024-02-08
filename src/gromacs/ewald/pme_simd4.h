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

/* This include file has code between ifdef's to make sure
 * that this performance sensitive code is inlined
 * and to remove conditionals and variable loop bounds at compile time.
 */

#ifdef PME_SPREAD_SIMD4_ORDER4
/* Spread one charge with pme_order=4 with unaligned SIMD4 load+store.
 * This code does not assume any memory alignment for the grid.
 */
{
    using namespace gmx;
    Simd4Real ty_S0(thy[0]);
    Simd4Real ty_S1(thy[1]);
    Simd4Real ty_S2(thy[2]);
    Simd4Real ty_S3(thy[3]);
    Simd4Real tz_S;
    Simd4Real vx_S;
    Simd4Real vx_tz_S;
    Simd4Real sum_S0, sum_S1, sum_S2, sum_S3;
    Simd4Real gri_S0, gri_S1, gri_S2, gri_S3;

    /* With order 4 the z-spline is actually aligned */
    tz_S = load4(thz);

    for (ithx = 0; (ithx < 4); ithx++)
    {
        index_x = (i0 + ithx) * pny * pnz;
        valx    = coefficient * thx[ithx];

        vx_S = Simd4Real(valx);

        vx_tz_S = vx_S * tz_S;

        gri_S0 = load4U(grid + index_x + (j0 + 0) * pnz + k0);
        gri_S1 = load4U(grid + index_x + (j0 + 1) * pnz + k0);
        gri_S2 = load4U(grid + index_x + (j0 + 2) * pnz + k0);
        gri_S3 = load4U(grid + index_x + (j0 + 3) * pnz + k0);

        sum_S0 = fma(vx_tz_S, ty_S0, gri_S0);
        sum_S1 = fma(vx_tz_S, ty_S1, gri_S1);
        sum_S2 = fma(vx_tz_S, ty_S2, gri_S2);
        sum_S3 = fma(vx_tz_S, ty_S3, gri_S3);

        store4U(grid + index_x + (j0 + 0) * pnz + k0, sum_S0);
        store4U(grid + index_x + (j0 + 1) * pnz + k0, sum_S1);
        store4U(grid + index_x + (j0 + 2) * pnz + k0, sum_S2);
        store4U(grid + index_x + (j0 + 3) * pnz + k0, sum_S3);
    }
}
#    undef PME_SPREAD_SIMD4_ORDER4
#endif


#ifdef PME_SPREAD_SIMD4_ALIGNED
/* This code assumes that the grid is allocated 4-real aligned
 * and that pnz is a multiple of 4.
 * This code supports pme_order <= 5.
 */
{
    using namespace gmx;
    int       offset;
    int       index;
    Simd4Real ty_S0(thy[0]);
    Simd4Real ty_S1(thy[1]);
    Simd4Real ty_S2(thy[2]);
    Simd4Real ty_S3(thy[3]);
    Simd4Real tz_S0;
    Simd4Real tz_S1;
    Simd4Real vx_S;
    Simd4Real vx_tz_S0;
    Simd4Real vx_tz_S1;
    Simd4Real sum_S00, sum_S01, sum_S02, sum_S03;
    Simd4Real sum_S10, sum_S11, sum_S12, sum_S13;
    Simd4Real gri_S00, gri_S01, gri_S02, gri_S03;
    Simd4Real gri_S10, gri_S11, gri_S12, gri_S13;
#    if PME_ORDER == 5
    Simd4Real ty_S4(thy[4]);
    Simd4Real sum_S04;
    Simd4Real sum_S14;
    Simd4Real gri_S04;
    Simd4Real gri_S14;
#    endif

    offset = k0 & 3;

#    ifdef PME_SIMD4_UNALIGNED
    tz_S0 = load4U(thz - offset);
    tz_S1 = load4U(thz - offset + 4);
#    else
    {
        int i;
        /* Copy thz to an aligned buffer (unused buffer parts are masked) */
        for (i = 0; i < PME_ORDER; i++)
        {
            thz_aligned[offset + i] = thz[i];
        }
        tz_S0 = load4(thz_aligned);
        tz_S1 = load4(thz_aligned + 4);
    }
#    endif
    tz_S0 = selectByMask(tz_S0, work.mask_S0[offset]);
    tz_S1 = selectByMask(tz_S1, work.mask_S1[offset]);

    for (ithx = 0; (ithx < PME_ORDER); ithx++)
    {
        index = (i0 + ithx) * pny * pnz + j0 * pnz + k0 - offset;
        valx  = coefficient * thx[ithx];

        vx_S = Simd4Real(valx);

        vx_tz_S0 = vx_S * tz_S0;
        vx_tz_S1 = vx_S * tz_S1;

        gri_S00 = load4(grid + index + 0 * pnz);
        gri_S01 = load4(grid + index + 1 * pnz);
        gri_S02 = load4(grid + index + 2 * pnz);
        gri_S03 = load4(grid + index + 3 * pnz);
#    if PME_ORDER == 5
        gri_S04 = load4(grid + index + 4 * pnz);
#    endif
        gri_S10 = load4(grid + index + 0 * pnz + 4);
        gri_S11 = load4(grid + index + 1 * pnz + 4);
        gri_S12 = load4(grid + index + 2 * pnz + 4);
        gri_S13 = load4(grid + index + 3 * pnz + 4);
#    if PME_ORDER == 5
        gri_S14 = load4(grid + index + 4 * pnz + 4);
#    endif

        sum_S00 = fma(vx_tz_S0, ty_S0, gri_S00);
        sum_S01 = fma(vx_tz_S0, ty_S1, gri_S01);
        sum_S02 = fma(vx_tz_S0, ty_S2, gri_S02);
        sum_S03 = fma(vx_tz_S0, ty_S3, gri_S03);
#    if PME_ORDER == 5
        sum_S04 = fma(vx_tz_S0, ty_S4, gri_S04);
#    endif
        sum_S10 = fma(vx_tz_S1, ty_S0, gri_S10);
        sum_S11 = fma(vx_tz_S1, ty_S1, gri_S11);
        sum_S12 = fma(vx_tz_S1, ty_S2, gri_S12);
        sum_S13 = fma(vx_tz_S1, ty_S3, gri_S13);
#    if PME_ORDER == 5
        sum_S14 = fma(vx_tz_S1, ty_S4, gri_S14);
#    endif

        store4(grid + index + 0 * pnz, sum_S00);
        store4(grid + index + 1 * pnz, sum_S01);
        store4(grid + index + 2 * pnz, sum_S02);
        store4(grid + index + 3 * pnz, sum_S03);
#    if PME_ORDER == 5
        store4(grid + index + 4 * pnz, sum_S04);
#    endif
        store4(grid + index + 0 * pnz + 4, sum_S10);
        store4(grid + index + 1 * pnz + 4, sum_S11);
        store4(grid + index + 2 * pnz + 4, sum_S12);
        store4(grid + index + 3 * pnz + 4, sum_S13);
#    if PME_ORDER == 5
        store4(grid + index + 4 * pnz + 4, sum_S14);
#    endif
    }
}
#    undef PME_ORDER
#    undef PME_SPREAD_SIMD4_ALIGNED
#endif
