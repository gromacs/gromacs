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
#ifndef _nbnxn_kernel_simd_utils_x86_128d_h_
#define _nbnxn_kernel_simd_utils_x86_128d_h_

#include "gromacs/legacyheaders/types/simple.h"

/* This files contains all functions/macros for the SIMD kernels
 * which have explicit dependencies on the j-cluster size and/or SIMD-width.
 * The functionality which depends on the j-cluster size is:
 *   LJ-parameter lookup
 *   force table lookup
 *   energy group pair energy storage
 */

#define gmx_mm_extract_epi32(x, imm) _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))

typedef gmx_simd_int32_t gmx_exclfilter;
/* This is set to a constant for now, since the code does not adapt automatically just
 * because we set the SIMD widths to other values.
 */
static const int filter_stride = 2;

/* Transpose 2 double precision registers */
static gmx_inline void
gmx_mm_transpose2_op_pd(__m128d in0, __m128d in1,
                        __m128d *out0, __m128d *out1)
{
    *out0 = _mm_unpacklo_pd(in0, in1);
    *out1 = _mm_unpackhi_pd(in0, in1);
}

/* Sum the elements within each input register and store the sums in out */
static gmx_inline __m128d
gmx_mm_transpose_sum2_pr(__m128d in0, __m128d in1)
{
    __m128d tr0, tr1;

    gmx_mm_transpose2_op_pd(in0, in1, &tr0, &tr1);

    return _mm_add_pd(tr0, tr1);
}

static gmx_inline __m128
gmx_mm128_invsqrt_ps_single(__m128 x)
{
    const __m128 half  = _mm_set_ps(0.5, 0.5, 0.5, 0.5);
    const __m128 three = _mm_set_ps(3.0, 3.0, 3.0, 3.0);

    __m128       lu = _mm_rsqrt_ps(x);

    return _mm_mul_ps(half, _mm_mul_ps(_mm_sub_ps(three, _mm_mul_ps(_mm_mul_ps(lu, lu), x)), lu));
}

/* Do 2 double precision invsqrt operations.
 * Doing the SIMD rsqrt and the first Newton Raphson iteration
 * in single precision gives full double precision accuracy.
 */
static gmx_inline void
gmx_mm_invsqrt2_pd(__m128d in0, __m128d in1,
                   __m128d *out0, __m128d *out1)
{
    const __m128d half  = _mm_set1_pd(0.5);
    const __m128d three = _mm_set1_pd(3.0);
    __m128        s, ir;
    __m128d       lu0, lu1;

    s     = _mm_movelh_ps(_mm_cvtpd_ps(in0), _mm_cvtpd_ps(in1));
    ir    = gmx_mm128_invsqrt_ps_single(s);
    lu0   = _mm_cvtps_pd(ir);
    lu1   = _mm_cvtps_pd(_mm_movehl_ps(ir, ir));
    *out0 = _mm_mul_pd(half, _mm_mul_pd(_mm_sub_pd(three, _mm_mul_pd(_mm_mul_pd(lu0, lu0), in0)), lu0));
    *out1 = _mm_mul_pd(half, _mm_mul_pd(_mm_sub_pd(three, _mm_mul_pd(_mm_mul_pd(lu1, lu1), in1)), lu1));
}

static gmx_inline void
load_lj_pair_params(const real *nbfp, const int *type, int aj,
                    __m128d *c6_S, __m128d *c12_S)
{
    __m128d clj_S[UNROLLJ];
    int     p;

    for (p = 0; p < UNROLLJ; p++)
    {
        clj_S[p] = _mm_load_pd(nbfp+type[aj+p]*nbfp_stride);
    }
    gmx_mm_transpose2_op_pd(clj_S[0], clj_S[1], c6_S, c12_S);
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

static gmx_inline void
load_table_f(const real *tab_coul_F, gmx_simd_int32_t ti_S, int gmx_unused *ti,
             __m128d *ctab0_S, __m128d *ctab1_S)
{
    int     idx[2];
    __m128d ctab_S[2];

    /* Without SSE4.1 the extract macro needs an immediate: unroll */
    idx[0]    = gmx_mm_extract_epi32(ti_S, 0);
    ctab_S[0] = _mm_loadu_pd(tab_coul_F+idx[0]);
    idx[1]    = gmx_mm_extract_epi32(ti_S, 1);
    ctab_S[1] = _mm_loadu_pd(tab_coul_F+idx[1]);

    /* Shuffle the force table entries to a convenient order */
    gmx_mm_transpose2_op_pd(ctab_S[0], ctab_S[1], ctab0_S, ctab1_S);
    /* The second force table entry should contain the difference */
    *ctab1_S = _mm_sub_pd(*ctab1_S, *ctab0_S);
}

static gmx_inline void
load_table_f_v(const real *tab_coul_F, const real *tab_coul_V,
               gmx_simd_int32_t ti_S, int gmx_unused *ti,
               __m128d *ctab0_S, __m128d *ctab1_S, __m128d *ctabv_S)
{
    int     idx[2];
    __m128d ctab_S[4];

    /* Without SSE4.1 the extract macro needs an immediate: unroll */
    idx[0]    = gmx_mm_extract_epi32(ti_S, 0);
    ctab_S[0] = _mm_loadu_pd(tab_coul_F+idx[0]);
    idx[1]    = gmx_mm_extract_epi32(ti_S, 1);
    ctab_S[1] = _mm_loadu_pd(tab_coul_F+idx[1]);

    /* Shuffle the force table entries to a convenient order */
    gmx_mm_transpose2_op_pd(ctab_S[0], ctab_S[1], ctab0_S, ctab1_S);
    /* The second force table entry should contain the difference */
    *ctab1_S = _mm_sub_pd(*ctab1_S, *ctab0_S);

    ctab_S[2] = _mm_loadu_pd(tab_coul_V+idx[0]);
    ctab_S[3] = _mm_loadu_pd(tab_coul_V+idx[1]);

    /* Shuffle the energy table entries to a single register */
    *ctabv_S = _mm_shuffle_pd(ctab_S[2], ctab_S[3], _MM_SHUFFLE2(0, 0));
}

static gmx_inline gmx_exclfilter
gmx_load1_exclfilter(int e)
{
    return _mm_set1_epi32(e);
}

static gmx_inline gmx_exclfilter
gmx_load_exclusion_filter(const unsigned *i)
{
    /* For now this has to be an explicit-float load since we use stride==2 */
    return gmx_simd_load_fi(i);
}

static gmx_inline gmx_simd_bool_t
gmx_checkbitmask_pb(gmx_exclfilter m0, gmx_exclfilter m1)
{
    return _mm_castsi128_pd(_mm_cmpeq_epi32(_mm_andnot_si128(m0, m1), _mm_setzero_si128()));
}

#endif /* _nbnxn_kernel_simd_utils_x86_s128d_h_ */
