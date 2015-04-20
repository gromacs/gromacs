/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
#ifndef _nbnxn_kernel_simd_utils_x86_128s_h_
#define _nbnxn_kernel_simd_utils_x86_128s_h_

#include "config.h"

#include "gromacs/legacyheaders/types/simple.h"

/* This files contains all functions/macros for the SIMD kernels
 * which have explicit dependencies on the j-cluster size and/or SIMD-width.
 * The functionality which depends on the j-cluster size is:
 *   LJ-parameter lookup
 *   force table lookup
 *   energy group pair energy storage
 */

typedef gmx_simd_int32_t gmx_exclfilter;
static const int filter_stride = GMX_SIMD_INT32_WIDTH/GMX_SIMD_REAL_WIDTH;

/* Collect element 0 and 1 of the 4 inputs to out0 and out1, respectively */
static gmx_inline void gmx_simdcall
gmx_shuffle_4_ps_fil01_to_2_ps(__m128 in0, __m128 in1, __m128 in2, __m128 in3,
                               __m128 *out0, __m128 *out1)
{
    __m128 _c01, _c23;

    _c01  = _mm_movelh_ps(in0, in1);
    _c23  = _mm_movelh_ps(in2, in3);
    *out0 = _mm_shuffle_ps(_c01, _c23, _MM_SHUFFLE(2, 0, 2, 0));
    *out1 = _mm_shuffle_ps(_c01, _c23, _MM_SHUFFLE(3, 1, 3, 1));
}

/* Collect element 2 of the 4 inputs to out */
static gmx_inline __m128 gmx_simdcall
gmx_shuffle_4_ps_fil2_to_1_ps(__m128 in0, __m128 in1, __m128 in2, __m128 in3)
{
    __m128 _c01, _c23;

    _c01 = _mm_shuffle_ps(in0, in1, _MM_SHUFFLE(3, 2, 3, 2));
    _c23 = _mm_shuffle_ps(in2, in3, _MM_SHUFFLE(3, 2, 3, 2));

    return _mm_shuffle_ps(_c01, _c23, _MM_SHUFFLE(2, 0, 2, 0));
}

/* Sum the elements within each input register and store the sums in out */
static gmx_inline __m128 gmx_simdcall
gmx_mm_transpose_sum4_pr(__m128 in0, __m128 in1,
                         __m128 in2, __m128 in3)
{
    _MM_TRANSPOSE4_PS(in0, in1, in2, in3);
    in0  = _mm_add_ps(in0, in1);
    in2  = _mm_add_ps(in2, in3);

    return _mm_add_ps(in0, in2);
}

static gmx_inline void
load_lj_pair_params(const real *nbfp, const int *type, int aj,
                    __m128 *c6_S, __m128 *c12_S)
{
    __m128 clj_S[UNROLLJ];
    int    p;

    for (p = 0; p < UNROLLJ; p++)
    {
        /* Here we load 4 aligned floats, but we need just 2 */
        clj_S[p] = gmx_simd_load_r(nbfp+type[aj+p]*nbfp_stride);
    }
    gmx_shuffle_4_ps_fil01_to_2_ps(clj_S[0], clj_S[1], clj_S[2], clj_S[3], c6_S, c12_S);
}

/* The load_table functions below are performance critical.
 * The routines issue UNROLLI*UNROLLJ _mm_load_ps calls.
 * As these all have latencies, scheduling is crucial.
 * The Intel compilers and CPUs seem to do a good job at this.
 * But AMD CPUs perform significantly worse with gcc than with icc.
 * Performance is improved a bit by using the extract function UNROLLJ times,
 * instead of doing an _mm_store_si128 for every i-particle.
 * This is only faster when we use FDV0 formatted tables, where we also need
 * to multiple the index by 4, which can be done by a SIMD bit shift.
 * With single precision AVX, 8 extracts are much slower than 1 store.
 * Because of this, the load_table_f function always takes the ti
 * parameter, which should contain a buffer that is aligned with
 * prepare_table_load_buffer(), but it is only used with full-width
 * AVX_256. */

static gmx_inline void gmx_simdcall
load_table_f(const real *tab_coul_FDV0, gmx_simd_int32_t ti_S, int gmx_unused *ti,
             __m128 *ctab0_S, __m128 *ctab1_S)
{
    int    idx[4];
    __m128 ctab_S[4];

    /* Table has 4 entries, left-shift index by 2 */
    ti_S = _mm_slli_epi32(ti_S, 2);
    /* Without SSE4.1 the extract macro needs an immediate: unroll */
    idx[0]    = gmx_simd_extract_i(ti_S, 0);
    ctab_S[0] = _mm_load_ps(tab_coul_FDV0+idx[0]);
    idx[1]    = gmx_simd_extract_i(ti_S, 1);
    ctab_S[1] = _mm_load_ps(tab_coul_FDV0+idx[1]);
    idx[2]    = gmx_simd_extract_i(ti_S, 2);
    ctab_S[2] = _mm_load_ps(tab_coul_FDV0+idx[2]);
    idx[3]    = gmx_simd_extract_i(ti_S, 3);
    ctab_S[3] = _mm_load_ps(tab_coul_FDV0+idx[3]);

    /* Shuffle the force table entries to a convenient order */
    gmx_shuffle_4_ps_fil01_to_2_ps(ctab_S[0], ctab_S[1], ctab_S[2], ctab_S[3], ctab0_S, ctab1_S);
}

static gmx_inline void gmx_simdcall
load_table_f_v(const real *tab_coul_FDV0, gmx_simd_int32_t ti_S, int gmx_unused *ti,
               __m128 *ctab0_S, __m128 *ctab1_S, __m128 *ctabv_S)
{
    int    idx[4];
    __m128 ctab_S[4];

    /* Table has 4 entries, left-shift index by 2 */
    ti_S = _mm_slli_epi32(ti_S, 2);
    /* Without SSE4.1 the extract macro needs an immediate: unroll */
    idx[0]    = gmx_simd_extract_i(ti_S, 0);
    ctab_S[0] = _mm_load_ps(tab_coul_FDV0+idx[0]);
    idx[1]    = gmx_simd_extract_i(ti_S, 1);
    ctab_S[1] = _mm_load_ps(tab_coul_FDV0+idx[1]);
    idx[2]    = gmx_simd_extract_i(ti_S, 2);
    ctab_S[2] = _mm_load_ps(tab_coul_FDV0+idx[2]);
    idx[3]    = gmx_simd_extract_i(ti_S, 3);
    ctab_S[3] = _mm_load_ps(tab_coul_FDV0+idx[3]);

    /* Shuffle the force table entries to a convenient order */
    gmx_shuffle_4_ps_fil01_to_2_ps(ctab_S[0], ctab_S[1], ctab_S[2], ctab_S[3], ctab0_S, ctab1_S);

    *ctabv_S = gmx_shuffle_4_ps_fil2_to_1_ps(ctab_S[0], ctab_S[1], ctab_S[2], ctab_S[3]);
}

static gmx_inline gmx_exclfilter gmx_simdcall
gmx_load1_exclfilter(int e)
{
    return _mm_set1_epi32(e);
}

static gmx_inline gmx_exclfilter gmx_simdcall
gmx_load_exclusion_filter(const unsigned *i)
{
    return gmx_simd_load_i(i);
}

static gmx_inline gmx_simd_bool_t gmx_simdcall
gmx_checkbitmask_pb(gmx_exclfilter m0, gmx_exclfilter m1)
{
    return _mm_castsi128_ps(_mm_cmpeq_epi32(_mm_andnot_si128(m0, m1), _mm_setzero_si128()));
}

#endif /* _nbnxn_kernel_simd_utils_x86_s128s_h_ */
