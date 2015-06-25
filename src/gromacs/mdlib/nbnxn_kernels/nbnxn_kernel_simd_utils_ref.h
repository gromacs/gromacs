/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
#ifndef _nbnxn_kernel_simd_utils_ref_h_
#define _nbnxn_kernel_simd_utils_ref_h_

#include "gromacs/mdlib/nbnxn_simd.h"
#include "gromacs/simd/simd_math.h"

typedef gmx_simd_int32_t        gmx_simd_ref_exclfilter;
typedef gmx_simd_ref_exclfilter gmx_exclfilter;
static const int filter_stride = GMX_SIMD_INT32_WIDTH/GMX_SIMD_REAL_WIDTH;

/* Set the stride for the lookup of the two LJ parameters from their
   (padded) array. Only strides of 2 and 4 are currently supported. */
#if defined GMX_NBNXN_SIMD_2XNN
static const int nbfp_stride = 4;
#elif defined GMX_DOUBLE
static const int nbfp_stride = 2;
#else
static const int nbfp_stride = 4;
#endif

#if GMX_SIMD_REAL_WIDTH > 4
/* The 4xn kernel operates on 4-wide i-force registers */

/* float/double SIMD register type */
typedef struct {
    real r[4];
} gmx_simd4_real_t;

static gmx_inline gmx_simd4_real_t
gmx_simd4_load_r(const real *r)
{
    gmx_simd4_real_t a;
    int              i;

    for (i = 0; i < 4; i++)
    {
        a.r[i] = r[i];
    }

    return a;
}

static gmx_inline void
gmx_simd4_store_r(real *dest, gmx_simd4_real_t src)
{
    gmx_simd4_real_t a;
    int              i;

    for (i = 0; i < 4; i++)
    {
        dest[i] = src.r[i];
    }
}

static gmx_inline gmx_simd4_real_t
gmx_simd4_add_r(gmx_simd4_real_t a, gmx_simd4_real_t b)
{
    gmx_simd4_real_t c;
    int              i;

    for (i = 0; i < 4; i++)
    {
        c.r[i] = a.r[i] + b.r[i];
    }

    return c;
}

static gmx_inline real
gmx_simd4_reduce_r(gmx_simd4_real_t a)
{
    return a.r[0] + a.r[1] + a.r[2] + a.r[3];
}

#endif


#ifdef GMX_NBNXN_SIMD_2XNN

/* Half-width operations are required for the 2xnn kernels */

/* Half-width SIMD real type */
/* float/double SIMD register type */
typedef struct {
    real r[GMX_SIMD_REAL_WIDTH/2];
} gmx_mm_hpr;

/* Half-width SIMD operations */

/* Load reals at half-width aligned pointer b into half-width SIMD register a */
static gmx_inline void
gmx_load_hpr(gmx_mm_hpr *a, const real *b)
{
    int i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        a->r[i] = b[i];
    }
}

/* Set all entries in half-width SIMD register *a to b */
static gmx_inline void
gmx_set1_hpr(gmx_mm_hpr *a, real b)
{
    int i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        a->r[i] = b;
    }
}

/* Load one real at b and one real at b+1 into halves of a, respectively */
static gmx_inline void
gmx_load1p1_pr(gmx_simd_real_t *a, const real *b)
{
    int i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        a->r[                        i] = b[0];
        a->r[GMX_SIMD_REAL_WIDTH/2 + i] = b[1];
    }
}

/* Load reals at half-width aligned pointer b into two halves of a */
static gmx_inline void
gmx_loaddh_pr(gmx_simd_real_t *a, const real *b)
{
    int i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        a->r[i]                         = b[i];
        a->r[GMX_SIMD_REAL_WIDTH/2 + i] = b[i];
    }
}

/* Store half-width SIMD register b into half width aligned memory a */
static gmx_inline void
gmx_store_hpr(real *a, gmx_mm_hpr b)
{
    int i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        a[i] = b.r[i];
    }
}

static gmx_inline gmx_mm_hpr
gmx_add_hpr(gmx_mm_hpr a, gmx_mm_hpr b)
{
    gmx_mm_hpr c;
    int        i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        c.r[i] = a.r[i] + b.r[i];
    }

    return c;
}

static gmx_inline gmx_mm_hpr
gmx_sub_hpr(gmx_mm_hpr a, gmx_mm_hpr b)
{
    gmx_mm_hpr c;
    int        i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        c.r[i] = a.r[i] - b.r[i];
    }

    return c;
}

/* Sum over 4 half SIMD registers */
static gmx_inline gmx_mm_hpr
gmx_sum4_hpr(gmx_simd_real_t a, gmx_simd_real_t b)
{
    gmx_mm_hpr c;
    int        i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        c.r[i] =
            a.r[i] +
            a.r[GMX_SIMD_REAL_WIDTH/2+i] +
            b.r[i] +
            b.r[GMX_SIMD_REAL_WIDTH/2+i];
    }

    return c;
}

#ifdef GMX_NBNXN_SIMD_2XNN
/* Sum the elements of halfs of each input register and store sums in out */
static gmx_inline gmx_simd4_real_t
gmx_mm_transpose_sum4h_pr(gmx_simd_real_t a, gmx_simd_real_t b)
{
    gmx_simd4_real_t sum;
    int              i;

    sum.r[0] = 0;
    sum.r[1] = 0;
    sum.r[2] = 0;
    sum.r[3] = 0;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        sum.r[0] += a.r[i];
        sum.r[1] += a.r[GMX_SIMD_REAL_WIDTH/2+i];
        sum.r[2] += b.r[i];
        sum.r[3] += b.r[GMX_SIMD_REAL_WIDTH/2+i];
    }

    return sum;
}
#endif

static gmx_inline void
gmx_pr_to_2hpr(gmx_simd_real_t a, gmx_mm_hpr *b, gmx_mm_hpr *c)
{
    int i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        b->r[i] = a.r[i];
        c->r[i] = a.r[GMX_SIMD_REAL_WIDTH/2 + i];
    }
}
static gmx_inline void
gmx_2hpr_to_pr(gmx_mm_hpr a, gmx_mm_hpr b, gmx_simd_real_t *c)
{
    int i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        c->r[i]                         = a.r[i];
        c->r[GMX_SIMD_REAL_WIDTH/2 + i] = b.r[i];
    }
}

#endif /* GMX_NBNXN_SIMD_2XNN */


#ifndef TAB_FDV0
static gmx_inline void
load_table_f(const real *tab_coul_F, gmx_simd_int32_t ti_S,
             int gmx_unused *ti,
             gmx_simd_real_t *ctab0_S, gmx_simd_real_t *ctab1_S)
{
    int i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        ctab0_S->r[i] = tab_coul_F[ti_S.i[i]];
        ctab1_S->r[i] = tab_coul_F[ti_S.i[i]+1];
    }

    *ctab1_S  = gmx_simd_sub_r(*ctab1_S, *ctab0_S);
}

static gmx_inline void
load_table_f_v(const real *tab_coul_F, const real *tab_coul_V,
               gmx_simd_int32_t ti_S, int *ti,
               gmx_simd_real_t *ctab0_S, gmx_simd_real_t *ctab1_S,
               gmx_simd_real_t *ctabv_S)
{
    int i;

    load_table_f(tab_coul_F, ti_S, ti, ctab0_S, ctab1_S);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        ctabv_S->r[i] = tab_coul_V[ti_S.i[i]];
    }
}
#endif

#ifdef TAB_FDV0
static gmx_inline void
load_table_f(const real *tab_coul_FDV0, gmx_simd_int32_t ti_S, int *ti,
             gmx_simd_real_t *ctab0_S, gmx_simd_real_t *ctab1_S)
{
    int i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        ctab0_S->r[i] = tab_coul_FDV0[ti_S.i[i]*4];
        ctab1_S->r[i] = tab_coul_FDV0[ti_S.i[i]*4+1];
    }
}

static gmx_inline void
load_table_f_v(const real *tab_coul_FDV0,
               gmx_simd_int32_t ti_S, int *ti,
               gmx_simd_real_t *ctab0_S, gmx_simd_real_t *ctab1_S,
               gmx_simd_real_t *ctabv_S)
{
    int i;

    load_table_f(tab_coul_FDV0, ti_S, ti, ctab0_S, ctab1_S);

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        ctabv_S->r[i] = tab_coul_FDV0[ti_S.i[i]*4+2];
    }
}
#endif

/* Sum the elements within each input register and store the sums in out.
 * Note that 4/8-way SIMD requires gmx_mm_transpose_sum4_pr instead.
 */
#if GMX_SIMD_REAL_WIDTH == 2
static gmx_inline gmx_simd_real_t
gmx_mm_transpose_sum2_pr(gmx_simd_real_t in0, gmx_simd_real_t in1)
{
    gmx_simd_real_t sum;

    sum.r[0] = in0.r[0] + in0.r[1];
    sum.r[1] = in1.r[0] + in1.r[1];

    return sum;
}
#endif

#if GMX_SIMD_REAL_WIDTH >= 4
#if GMX_SIMD_REAL_WIDTH == 4
static gmx_inline gmx_simd_real_t
#else
static gmx_inline gmx_simd4_real_t
#endif
gmx_mm_transpose_sum4_pr(gmx_simd_real_t in0, gmx_simd_real_t in1,
                         gmx_simd_real_t in2, gmx_simd_real_t in3)
{
#if GMX_SIMD_REAL_WIDTH == 4
    gmx_simd_real_t  sum;
#else
    gmx_simd4_real_t sum;
#endif
    int              i;

    sum.r[0] = 0;
    sum.r[1] = 0;
    sum.r[2] = 0;
    sum.r[3] = 0;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
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
gmx_mm_invsqrt2_pd(gmx_simd_real_t in0, gmx_simd_real_t in1,
                   gmx_simd_real_t *out0, gmx_simd_real_t *out1)
{
    *out0 = gmx_simd_invsqrt_r(in0);
    *out1 = gmx_simd_invsqrt_r(in1);
}
#endif

static gmx_inline void
load_lj_pair_params(const real *nbfp, const int *type, int aj,
                    gmx_simd_real_t *c6_S, gmx_simd_real_t *c12_S)
{
    int i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        c6_S->r[i]  = nbfp[type[aj+i]*nbfp_stride];
        c12_S->r[i] = nbfp[type[aj+i]*nbfp_stride+1];
    }
}

#ifdef GMX_NBNXN_SIMD_2XNN
static gmx_inline void
load_lj_pair_params2(const real *nbfp0, const real *nbfp1,
                     const int *type, int aj,
                     gmx_simd_real_t *c6_S, gmx_simd_real_t *c12_S)
{
    int i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH/2; i++)
    {
        c6_S->r[i]                          = nbfp0[type[aj+i]*nbfp_stride];
        c6_S->r[GMX_SIMD_REAL_WIDTH/2 + i]  = nbfp1[type[aj+i]*nbfp_stride];
        c12_S->r[i]                         = nbfp0[type[aj+i]*nbfp_stride+1];
        c12_S->r[GMX_SIMD_REAL_WIDTH/2 + i] = nbfp1[type[aj+i]*nbfp_stride+1];
    }
}
#endif

/* Code for handling loading exclusions and converting them into
   interactions. The x86 code might use either integer- or real-type
   SIMD, but the reference code does not need to know. */

#define gmx_load1_exclfilter(e)       gmx_simd_ref_load1_exclfilter(e)
#define gmx_load_exclusion_filter(e)  gmx_simd_ref_load_exclusion_filter((int *) e)
#define gmx_checkbitmask_pb(m0, m1)   gmx_simd_ref_checkbitmask_pb(m0, m1)

static gmx_inline gmx_simd_ref_exclfilter
gmx_simd_ref_load1_exclfilter(int src)
{
    gmx_simd_ref_exclfilter a;
    int                     i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        a.i[i] = src;
    }

    return a;
}

static gmx_inline gmx_simd_ref_exclfilter
gmx_simd_ref_load_exclusion_filter(const int *src)
{
    gmx_simd_ref_exclfilter a;
    int                     i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        a.i[i] = src[i];
    }

    return a;
}

/* For topology exclusion-pair checking we need: ((a & b) ? True :
 * False). The x86 implementations use hardware-suitable integer-
 * and/or real-valued SIMD operations and a bit-wise "and" to do
 * this. The reference implementation normally uses logical operations
 * for logic, but in this case the i- and j-atom exclusion masks
 * computed during searching expect to be combined with bit-wise
 * "and".
 *
 * If the same bit is set in both input masks, return TRUE, else
 * FALSE. This function is only called with a single bit set in b.
 */
static gmx_inline gmx_simd_bool_t
gmx_simd_ref_checkbitmask_pb(gmx_simd_ref_exclfilter a, gmx_simd_ref_exclfilter b)
{
    gmx_simd_bool_t c;
    int             i;

    for (i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        c.b[i] = ((a.i[i] & b.i[i]) != 0);
    }

    return c;
}

#endif /* _nbnxn_kernel_simd_utils_ref_h_ */
