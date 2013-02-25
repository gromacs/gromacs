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
#ifndef _nbnxn_kernel_sse_utils_h_
#define _nbnxn_kernel_sse_utils_h_

/* This files contains all functions/macros for the SIMD kernels
 * which have explicit dependencies on the j-cluster size and/or SIMD-width.
 * The functionality which depends on the j-cluster size is:
 *   LJ-parameter lookup
 *   force table lookup
 *   energy group pair energy storage
 */


#ifdef GMX_SIMD_PLAIN_C

static gmx_inline void
load_table_f(const real *tab_coul_F, gmx_epi32 ti_S, int *ti,
             gmx_mm_pr *ctab0_S, gmx_mm_pr *ctab1_S)
{
    ctab0_S->x = tab_coul_F[ti_S.x];
    ctab1_S->x = tab_coul_F[ti_S.x+1];
    ctab0_S->y = tab_coul_F[ti_S.y];
    ctab1_S->y = tab_coul_F[ti_S.y+1];

    *ctab1_S  = gmx_sub_pr(*ctab1_S, *ctab0_S);
}

static gmx_inline void
load_table_f_v(const real *tab_coul_F, const real *tab_coul_V,
               gmx_epi32 ti_S, int *ti,
               gmx_mm_pr *ctab0_S, gmx_mm_pr *ctab1_S, gmx_mm_pr *ctabv_S)
{
    load_table_f(tab_coul_F, ti_S, ti, ctab0_S, ctab1_S);

    ctabv_S->x = tab_coul_V[ti_S.x];
    ctabv_S->y = tab_coul_V[ti_S.y];
}

/* Sum the elements within each input register and store the sums in out.
 * Note that 4/8-way SIMD requires gmx_mm_transpose_sum4_pr instead.
 */
static gmx_inline gmx_mm_pr
gmx_mm_transpose_sum2_pr(gmx_mm_pr in0, gmx_mm_pr in1)
{
    real tmp;

    tmp   = in0.y;
    in0.y = in1.x;
    in1.x = tmp;

    return gmx_add_pr(in0, in1);
}

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
    out0->x = 1.0/sqrt(in0.x);
    out0->y = 1.0/sqrt(in0.y);
    out1->x = 1.0/sqrt(in1.x);
    out1->y = 1.0/sqrt(in1.y);
}
#endif

static gmx_inline void
load_lj_pair_params(const real *nbfp, const int *type, int aj,
                    gmx_mm_pr *c6_S, gmx_mm_pr *c12_S)
{
    c6_S->x  = nbfp[type[aj+0]*NBFP_STRIDE];
    c6_S->y  = nbfp[type[aj+1]*NBFP_STRIDE];
    c12_S->x = nbfp[type[aj+0]*NBFP_STRIDE+1];
    c12_S->y = nbfp[type[aj+1]*NBFP_STRIDE+1];
}

#endif /* GMX_SIMD_PLAIN_C */


/* Include SIMD architecture specific versions of the 4/5 functions above */
#ifdef GMX_X86_SSE2
/* Include x86 SSE2 compatible SIMD functions */
#if defined GMX_X86_AVX_256 && !defined GMX_USE_HALF_WIDTH_SIMD_HERE
#ifdef GMX_DOUBLE
#include "nbnxn_kernel_simd_utils_x86_256d.h"
#else
#include "nbnxn_kernel_simd_utils_x86_256s.h"
#endif
#else
#ifdef GMX_DOUBLE
#include "nbnxn_kernel_simd_utils_x86_128d.h"
#else
#include "nbnxn_kernel_simd_utils_x86_128s.h"
#endif
#endif
#endif


/* Add energy register to possibly multiple terms in the energy array */
static inline void add_ener_grp(gmx_mm_pr e_S, real *v, const int *offset_jj)
{
    int jj;

    /* We need to balance the number of store operations with
     * the rapidly increases number of combinations of energy groups.
     * We add to a temporary buffer for 1 i-group vs 2 j-groups.
     */
    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        gmx_mm_pr v_S;

        v_S = gmx_load_pr(v+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE);
        gmx_store_pr(v+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE, gmx_add_pr(v_S, e_S));
    }
}

/* Note that for 2xnn kernels add_ener_grp_halves also needs to be implemented.
 * For 256bit AVX this is in nbnxn_kernel_simd_utils_x86_256s.h
 */

#endif /* _nbnxn_kernel_sse_utils_h_ */
