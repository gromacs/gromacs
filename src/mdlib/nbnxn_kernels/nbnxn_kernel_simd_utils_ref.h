/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS Development Team
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
#ifndef _nbnxn_kernel_simd_utils_ref_h_
#define _nbnxn_kernel_simd_utils_ref_h_

/* This files contains all functions/macros for the SIMD kernels
 * which have explicit dependencies on the j-cluster size and/or SIMD-width.
 * The functionality which depends on the j-cluster size is:
 *   LJ-parameter lookup
 *   force table lookup
 *   energy group pair energy storage
 */


#if GMX_SIMD_WIDTH_HERE > 4
/* The 4xn kernel operates on 4-wide i-force registers */

/* float/double SIMD register type */
typedef struct {
    real r[4];
} gmx_mm_pr4;

static gmx_inline gmx_mm_pr4
gmx_load_pr4(const real *r)
{
    gmx_mm_pr4 a;
    int        i;

    for (i = 0; i < 4; i++)
    {
        a.r[i] = r[i];
    }

    return a;
}

static gmx_inline void
gmx_store_pr4(real *dest, gmx_mm_pr4 src)
{
    gmx_mm_pr4 a;
    int        i;

    for (i = 0; i < 4; i++)
    {
        dest[i] = src.r[i];
    }
}

static gmx_inline gmx_mm_pr4
gmx_add_pr4(gmx_mm_pr4 a, gmx_mm_pr4 b)
{
    gmx_mm_pr4 c;
    int        i;

    for (i = 0; i < 4; i++)
    {
        c.r[i] = a.r[i] + b.r[i];
    }

    return c;
}

#endif


#ifndef TAB_FDV0
static gmx_inline void
load_table_f(const real *tab_coul_F, gmx_epi32 ti_S, int *ti,
             gmx_mm_pr *ctab0_S, gmx_mm_pr *ctab1_S)
{
    int i;

    for (i = 0; i < GMX_SIMD_WIDTH_HERE; i++)
    {
        ctab0_S->r[i] = tab_coul_F[ti_S.r[i]];
        ctab1_S->r[i] = tab_coul_F[ti_S.r[i]+1];
    }

    *ctab1_S  = gmx_sub_pr(*ctab1_S, *ctab0_S);
}

static gmx_inline void
load_table_f_v(const real *tab_coul_F, const real *tab_coul_V,
               gmx_epi32 ti_S, int *ti,
               gmx_mm_pr *ctab0_S, gmx_mm_pr *ctab1_S, gmx_mm_pr *ctabv_S)
{
    int i;

    load_table_f(tab_coul_F, ti_S, ti, ctab0_S, ctab1_S);

    for (i = 0; i < GMX_SIMD_WIDTH_HERE; i++)
    {
        ctabv_S->r[i] = tab_coul_V[ti_S.r[i]];
    }
}
#endif

#ifdef TAB_FDV0
static gmx_inline void
load_table_f(const real *tab_coul_FDV0, gmx_epi32 ti_S, int *ti,
             gmx_mm_pr *ctab0_S, gmx_mm_pr *ctab1_S)
{
    int i;

    for (i = 0; i < GMX_SIMD_WIDTH_HERE; i++)
    {
        ctab0_S->r[i] = tab_coul_FDV0[ti_S.r[i]*4];
        ctab1_S->r[i] = tab_coul_FDV0[ti_S.r[i]*4+1];
    }
}

static gmx_inline void
load_table_f_v(const real *tab_coul_FDV0,
               gmx_epi32 ti_S, int *ti,
               gmx_mm_pr *ctab0_S, gmx_mm_pr *ctab1_S, gmx_mm_pr *ctabv_S)
{
    int i;

    load_table_f(tab_coul_FDV0, ti_S, ti, ctab0_S, ctab1_S);

    for (i = 0; i < GMX_SIMD_WIDTH_HERE; i++)
    {
        ctabv_S->r[i] = tab_coul_FDV0[ti_S.r[i]*4+2];
    }
}
#endif

/* Sum the elements within each input register and store the sums in out.
 * Note that 4/8-way SIMD requires gmx_mm_transpose_sum4_pr instead.
 */
#if GMX_SIMD_WIDTH_HERE == 2
static gmx_inline gmx_mm_pr
gmx_mm_transpose_sum2_pr(gmx_mm_pr in0, gmx_mm_pr in1)
{
    gmx_mm_pr sum;

    sum.r[0] = in0.r[0] + in0.r[1];
    sum.r[1] = in1.r[0] + in1.r[1];

    return sum;
}
#endif

#if GMX_SIMD_WIDTH_HERE >= 4
#if GMX_SIMD_WIDTH_HERE == 4
static gmx_inline gmx_mm_pr
#else
static gmx_inline gmx_mm_pr4
#endif
gmx_mm_transpose_sum4_pr(gmx_mm_pr in0, gmx_mm_pr in1,
                         gmx_mm_pr in2, gmx_mm_pr in3)
{
#if GMX_SIMD_WIDTH_HERE == 4
    gmx_mm_pr  sum;
#else
    gmx_mm_pr4 sum;
#endif
    int        i;

    sum.r[0] = 0;
    sum.r[1] = 0;
    sum.r[2] = 0;
    sum.r[3] = 0;

    for (i = 0; i < GMX_SIMD_WIDTH_HERE; i++)
    {
        sum.r[0] += in0.r[i];
        sum.r[1] += in1.r[i];
        sum.r[2] += in2.r[i];
        sum.r[3] += in3.r[i];
    }

    return sum;
}
#endif

#ifdef GMX_DOUBLE
/* In double precision it can be faster to first calculate single precision
 * square roots for two double precision registers at once and then use
 * double precision Newton-Raphson iteration to reach full double precision.
 * For this reference code we just use a plain-C sqrt.
 */
static gmx_inline void
gmx_mm_invsqrt2_pd(gmx_mm_pr in0, gmx_mm_pr in1,
                   gmx_mm_pr *out0, gmx_mm_pr *out1)
{
    out0 = gmx_invsqrt_pr(in0);
    out1 = gmx_invsqrt_pr(in1);
}
#endif

static gmx_inline void
load_lj_pair_params(const real *nbfp, const int *type, int aj,
                    gmx_mm_pr *c6_S, gmx_mm_pr *c12_S)
{
    int i;

    for (i = 0; i < GMX_SIMD_WIDTH_HERE; i++)
    {
        c6_S->r[i]  = nbfp[type[aj+i]*NBFP_STRIDE];
        c12_S->r[i] = nbfp[type[aj+i]*NBFP_STRIDE+1];
    }
}

#endif /* _nbnxn_kernel_simd_utils_ref_h_ */
