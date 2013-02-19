/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifdef PME_SPREAD_SSE_ORDER4
/* This code does not assume any memory alignment.
 * This code only works for pme_order = 4.
 */
{
    __m128 ty_SSE0, ty_SSE1, ty_SSE2, ty_SSE3;
    __m128 tz_SSE;
    __m128 vx_SSE;
    __m128 vx_tz_SSE;
    __m128 sum_SSE0, sum_SSE1, sum_SSE2, sum_SSE3;
    __m128 gri_SSE0, gri_SSE1, gri_SSE2, gri_SSE3;

    ty_SSE0 = _mm_load1_ps(&thy[0]);
    ty_SSE1 = _mm_load1_ps(&thy[1]);
    ty_SSE2 = _mm_load1_ps(&thy[2]);
    ty_SSE3 = _mm_load1_ps(&thy[3]);

    tz_SSE  = _mm_loadu_ps(thz);

    for (ithx = 0; (ithx < 4); ithx++)
    {
        index_x = (i0+ithx)*pny*pnz;
        valx    = qn*thx[ithx];

        vx_SSE   = _mm_load1_ps(&valx);

        vx_tz_SSE = _mm_mul_ps(vx_SSE, tz_SSE);

        gri_SSE0 = _mm_loadu_ps(grid+index_x+(j0+0)*pnz+k0);
        gri_SSE1 = _mm_loadu_ps(grid+index_x+(j0+1)*pnz+k0);
        gri_SSE2 = _mm_loadu_ps(grid+index_x+(j0+2)*pnz+k0);
        gri_SSE3 = _mm_loadu_ps(grid+index_x+(j0+3)*pnz+k0);

        sum_SSE0 = _mm_add_ps(gri_SSE0, _mm_mul_ps(vx_tz_SSE, ty_SSE0));
        sum_SSE1 = _mm_add_ps(gri_SSE1, _mm_mul_ps(vx_tz_SSE, ty_SSE1));
        sum_SSE2 = _mm_add_ps(gri_SSE2, _mm_mul_ps(vx_tz_SSE, ty_SSE2));
        sum_SSE3 = _mm_add_ps(gri_SSE3, _mm_mul_ps(vx_tz_SSE, ty_SSE3));

        _mm_storeu_ps(grid+index_x+(j0+0)*pnz+k0, sum_SSE0);
        _mm_storeu_ps(grid+index_x+(j0+1)*pnz+k0, sum_SSE1);
        _mm_storeu_ps(grid+index_x+(j0+2)*pnz+k0, sum_SSE2);
        _mm_storeu_ps(grid+index_x+(j0+3)*pnz+k0, sum_SSE3);
    }
}
#undef PME_SPREAD_SSE_ORDER4
#endif


#ifdef PME_GATHER_F_SSE_ORDER4
/* This code does not assume any memory alignment.
 * This code only works for pme_order = 4.
 */
{
    float  fx_tmp[4], fy_tmp[4], fz_tmp[4];

    __m128 fx_SSE, fy_SSE, fz_SSE;

    __m128 tx_SSE, ty_SSE, tz_SSE;
    __m128 dx_SSE, dy_SSE, dz_SSE;

    __m128 gval_SSE;

    __m128 fxy1_SSE;
    __m128 fz1_SSE;

    fx_SSE = _mm_setzero_ps();
    fy_SSE = _mm_setzero_ps();
    fz_SSE = _mm_setzero_ps();

    tz_SSE  = _mm_loadu_ps(thz);
    dz_SSE  = _mm_loadu_ps(dthz);

    for (ithx = 0; (ithx < 4); ithx++)
    {
        index_x  = (i0+ithx)*pny*pnz;
        tx_SSE   = _mm_load1_ps(thx+ithx);
        dx_SSE   = _mm_load1_ps(dthx+ithx);

        for (ithy = 0; (ithy < 4); ithy++)
        {
            index_xy = index_x+(j0+ithy)*pnz;
            ty_SSE   = _mm_load1_ps(thy+ithy);
            dy_SSE   = _mm_load1_ps(dthy+ithy);

            gval_SSE = _mm_loadu_ps(grid+index_xy+k0);

            fxy1_SSE = _mm_mul_ps(tz_SSE, gval_SSE);
            fz1_SSE  = _mm_mul_ps(dz_SSE, gval_SSE);

            fx_SSE = _mm_add_ps(fx_SSE, _mm_mul_ps(_mm_mul_ps(dx_SSE, ty_SSE), fxy1_SSE));
            fy_SSE = _mm_add_ps(fy_SSE, _mm_mul_ps(_mm_mul_ps(tx_SSE, dy_SSE), fxy1_SSE));
            fz_SSE = _mm_add_ps(fz_SSE, _mm_mul_ps(_mm_mul_ps(tx_SSE, ty_SSE), fz1_SSE));
        }
    }

    _mm_storeu_ps(fx_tmp, fx_SSE);
    _mm_storeu_ps(fy_tmp, fy_SSE);
    _mm_storeu_ps(fz_tmp, fz_SSE);

    fx += fx_tmp[0]+fx_tmp[1]+fx_tmp[2]+fx_tmp[3];
    fy += fy_tmp[0]+fy_tmp[1]+fy_tmp[2]+fy_tmp[3];
    fz += fz_tmp[0]+fz_tmp[1]+fz_tmp[2]+fz_tmp[3];
}
#undef PME_GATHER_F_SSE_ORDER4
#endif


#ifdef PME_SPREAD_SSE_ALIGNED
/* This code assumes that the grid is allocated 16 bit aligned
 * and that pnz is a multiple of 4.
 * This code supports pme_order <= 5.
 */
{
    int    offset;
    int    index;
    __m128 ty_SSE0, ty_SSE1, ty_SSE2, ty_SSE3, ty_SSE4;
    __m128 tz_SSE0;
    __m128 tz_SSE1;
    __m128 vx_SSE;
    __m128 vx_tz_SSE0;
    __m128 vx_tz_SSE1;
    __m128 sum_SSE00, sum_SSE01, sum_SSE02, sum_SSE03, sum_SSE04;
    __m128 sum_SSE10, sum_SSE11, sum_SSE12, sum_SSE13, sum_SSE14;
    __m128 gri_SSE00, gri_SSE01, gri_SSE02, gri_SSE03, gri_SSE04;
    __m128 gri_SSE10, gri_SSE11, gri_SSE12, gri_SSE13, gri_SSE14;

    offset = k0 & 3;

    ty_SSE0 = _mm_load1_ps(&thy[0]);
    ty_SSE1 = _mm_load1_ps(&thy[1]);
    ty_SSE2 = _mm_load1_ps(&thy[2]);
    ty_SSE3 = _mm_load1_ps(&thy[3]);
#if PME_ORDER == 5
    ty_SSE4 = _mm_load1_ps(&thy[4]);
#endif

    tz_SSE0 = _mm_loadu_ps(thz-offset);
    tz_SSE1 = _mm_loadu_ps(thz-offset+4);
    tz_SSE0 = _mm_and_ps(tz_SSE0, work->mask_SSE0[offset]);
    tz_SSE1 = _mm_and_ps(tz_SSE1, work->mask_SSE1[offset]);

    for (ithx = 0; (ithx < PME_ORDER); ithx++)
    {
        index = (i0+ithx)*pny*pnz + j0*pnz + k0 - offset;
        valx  = qn*thx[ithx];

        vx_SSE   = _mm_load1_ps(&valx);

        vx_tz_SSE0 = _mm_mul_ps(vx_SSE, tz_SSE0);
        vx_tz_SSE1 = _mm_mul_ps(vx_SSE, tz_SSE1);

        gri_SSE00 = _mm_load_ps(grid+index+0*pnz);
        gri_SSE01 = _mm_load_ps(grid+index+1*pnz);
        gri_SSE02 = _mm_load_ps(grid+index+2*pnz);
        gri_SSE03 = _mm_load_ps(grid+index+3*pnz);
#if PME_ORDER == 5
        gri_SSE04 = _mm_load_ps(grid+index+4*pnz);
#endif
        gri_SSE10 = _mm_load_ps(grid+index+0*pnz+4);
        gri_SSE11 = _mm_load_ps(grid+index+1*pnz+4);
        gri_SSE12 = _mm_load_ps(grid+index+2*pnz+4);
        gri_SSE13 = _mm_load_ps(grid+index+3*pnz+4);
#if PME_ORDER == 5
        gri_SSE14 = _mm_load_ps(grid+index+4*pnz+4);
#endif

        sum_SSE00 = _mm_add_ps(gri_SSE00, _mm_mul_ps(vx_tz_SSE0, ty_SSE0));
        sum_SSE01 = _mm_add_ps(gri_SSE01, _mm_mul_ps(vx_tz_SSE0, ty_SSE1));
        sum_SSE02 = _mm_add_ps(gri_SSE02, _mm_mul_ps(vx_tz_SSE0, ty_SSE2));
        sum_SSE03 = _mm_add_ps(gri_SSE03, _mm_mul_ps(vx_tz_SSE0, ty_SSE3));
#if PME_ORDER == 5
        sum_SSE04 = _mm_add_ps(gri_SSE04, _mm_mul_ps(vx_tz_SSE0, ty_SSE4));
#endif
        sum_SSE10 = _mm_add_ps(gri_SSE10, _mm_mul_ps(vx_tz_SSE1, ty_SSE0));
        sum_SSE11 = _mm_add_ps(gri_SSE11, _mm_mul_ps(vx_tz_SSE1, ty_SSE1));
        sum_SSE12 = _mm_add_ps(gri_SSE12, _mm_mul_ps(vx_tz_SSE1, ty_SSE2));
        sum_SSE13 = _mm_add_ps(gri_SSE13, _mm_mul_ps(vx_tz_SSE1, ty_SSE3));
#if PME_ORDER == 5
        sum_SSE14 = _mm_add_ps(gri_SSE14, _mm_mul_ps(vx_tz_SSE1, ty_SSE4));
#endif

        _mm_store_ps(grid+index+0*pnz, sum_SSE00);
        _mm_store_ps(grid+index+1*pnz, sum_SSE01);
        _mm_store_ps(grid+index+2*pnz, sum_SSE02);
        _mm_store_ps(grid+index+3*pnz, sum_SSE03);
#if PME_ORDER == 5
        _mm_store_ps(grid+index+4*pnz, sum_SSE04);
#endif
        _mm_store_ps(grid+index+0*pnz+4, sum_SSE10);
        _mm_store_ps(grid+index+1*pnz+4, sum_SSE11);
        _mm_store_ps(grid+index+2*pnz+4, sum_SSE12);
        _mm_store_ps(grid+index+3*pnz+4, sum_SSE13);
#if PME_ORDER == 5
        _mm_store_ps(grid+index+4*pnz+4, sum_SSE14);
#endif
    }
}
#undef PME_ORDER
#undef PME_SPREAD_SSE_ALIGNED
#endif


#ifdef PME_GATHER_F_SSE_ALIGNED
/* This code assumes that the grid is allocated 16 bit aligned
 * and that pnz is a multiple of 4.
 * This code supports pme_order <= 5.
 */
{
    int    offset;

    float  fx_tmp[4], fy_tmp[4], fz_tmp[4];

    __m128 fx_SSE, fy_SSE, fz_SSE;

    __m128 tx_SSE, ty_SSE, tz_SSE0, tz_SSE1;
    __m128 dx_SSE, dy_SSE, dz_SSE0, dz_SSE1;

    __m128 gval_SSE0;
    __m128 gval_SSE1;

    __m128 fxy1_SSE0;
    __m128 fz1_SSE0;
    __m128 fxy1_SSE1;
    __m128 fz1_SSE1;
    __m128 fxy1_SSE;
    __m128 fz1_SSE;

    offset = k0 & 3;

    fx_SSE = _mm_setzero_ps();
    fy_SSE = _mm_setzero_ps();
    fz_SSE = _mm_setzero_ps();

    tz_SSE0 = _mm_loadu_ps(thz-offset);
    dz_SSE0 = _mm_loadu_ps(dthz-offset);
    tz_SSE1 = _mm_loadu_ps(thz-offset+4);
    dz_SSE1 = _mm_loadu_ps(dthz-offset+4);
    tz_SSE0 = _mm_and_ps(tz_SSE0, work->mask_SSE0[offset]);
    dz_SSE0 = _mm_and_ps(dz_SSE0, work->mask_SSE0[offset]);
    tz_SSE1 = _mm_and_ps(tz_SSE1, work->mask_SSE1[offset]);
    dz_SSE1 = _mm_and_ps(dz_SSE1, work->mask_SSE1[offset]);

    for (ithx = 0; (ithx < PME_ORDER); ithx++)
    {
        index_x  = (i0+ithx)*pny*pnz;
        tx_SSE   = _mm_load1_ps(thx+ithx);
        dx_SSE   = _mm_load1_ps(dthx+ithx);

        for (ithy = 0; (ithy < PME_ORDER); ithy++)
        {
            index_xy = index_x+(j0+ithy)*pnz;
            ty_SSE   = _mm_load1_ps(thy+ithy);
            dy_SSE   = _mm_load1_ps(dthy+ithy);

            gval_SSE0 = _mm_load_ps(grid+index_xy+k0-offset);
            gval_SSE1 = _mm_load_ps(grid+index_xy+k0-offset+4);

            fxy1_SSE0 = _mm_mul_ps(tz_SSE0, gval_SSE0);
            fz1_SSE0  = _mm_mul_ps(dz_SSE0, gval_SSE0);
            fxy1_SSE1 = _mm_mul_ps(tz_SSE1, gval_SSE1);
            fz1_SSE1  = _mm_mul_ps(dz_SSE1, gval_SSE1);

            fxy1_SSE = _mm_add_ps(fxy1_SSE0, fxy1_SSE1);
            fz1_SSE  = _mm_add_ps(fz1_SSE0, fz1_SSE1);

            fx_SSE = _mm_add_ps(fx_SSE, _mm_mul_ps(_mm_mul_ps(dx_SSE, ty_SSE), fxy1_SSE));
            fy_SSE = _mm_add_ps(fy_SSE, _mm_mul_ps(_mm_mul_ps(tx_SSE, dy_SSE), fxy1_SSE));
            fz_SSE = _mm_add_ps(fz_SSE, _mm_mul_ps(_mm_mul_ps(tx_SSE, ty_SSE), fz1_SSE));
        }
    }

    _mm_store_ps(fx_tmp, fx_SSE);
    _mm_store_ps(fy_tmp, fy_SSE);
    _mm_store_ps(fz_tmp, fz_SSE);

    fx += fx_tmp[0]+fx_tmp[1]+fx_tmp[2]+fx_tmp[3];
    fy += fy_tmp[0]+fy_tmp[1]+fy_tmp[2]+fy_tmp[3];
    fz += fz_tmp[0]+fz_tmp[1]+fz_tmp[2]+fz_tmp[3];
}
#undef PME_ORDER
#undef PME_GATHER_F_SSE_ALIGNED
#endif
