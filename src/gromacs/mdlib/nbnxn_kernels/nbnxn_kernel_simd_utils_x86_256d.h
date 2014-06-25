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
#ifndef _nbnxn_kernel_simd_utils_x86_256d_h_
#define _nbnxn_kernel_simd_utils_x86_256d_h_

/* This files contains all functions/macros for the SIMD kernels
 * which have explicit dependencies on the j-cluster size and/or SIMD-width.
 * The functionality which depends on the j-cluster size is:
 *   LJ-parameter lookup
 *   force table lookup
 *   energy group pair energy storage
 */

typedef gmx_simd_real_t gmx_exclfilter;
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
static gmx_inline __m256d
gmx_mm_transpose_sum4_pr(__m256d in0, __m256d in1,
                         __m256d in2, __m256d in3)
{
    in0 = _mm256_hadd_pd(in0, in1);
    in2 = _mm256_hadd_pd(in2, in3);

    return _mm256_add_pd(_mm256_permute2f128_pd(in0, in2, 0x20), _mm256_permute2f128_pd(in0, in2, 0x31));
}

static gmx_inline __m256
gmx_mm256_invsqrt_ps_single(__m256 x)
{
    const __m256 half  = _mm256_set_ps(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5);
    const __m256 three = _mm256_set_ps(3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0);

    __m256       lu = _mm256_rsqrt_ps(x);

    return _mm256_mul_ps(half, _mm256_mul_ps(_mm256_sub_ps(three, _mm256_mul_ps(_mm256_mul_ps(lu, lu), x)), lu));
}

/* Put two 128-bit 4-float registers into one 256-bit 8-float register */
static gmx_inline __m256
gmx_2_m128_to_m256(__m128 in0, __m128 in1)
{
    return _mm256_insertf128_ps(_mm256_castps128_ps256(in0), in1, 1);
}

/* Put two 128-bit 2-double registers into one 256-bit 4-double register */
static gmx_inline __m256d
gmx_2_m128d_to_m256d(__m128d in0, __m128d in1)
{
    return _mm256_insertf128_pd(_mm256_castpd128_pd256(in0), in1, 1);
}

/* Do 2 double precision invsqrt operations.
 * Doing the SIMD rsqrt and the first Newton Raphson iteration
 * in single precision gives full double precision accuracy.
 */
static gmx_inline void
gmx_mm_invsqrt2_pd(__m256d in0, __m256d in1,
                   __m256d *out0, __m256d *out1)
{
    const __m256d half  = _mm256_set1_pd(0.5);
    const __m256d three = _mm256_set1_pd(3.0);
    __m256        s, ir;
    __m256d       lu0, lu1;

    s     =  gmx_2_m128_to_m256(_mm256_cvtpd_ps(in0), _mm256_cvtpd_ps(in1));
    ir    = gmx_mm256_invsqrt_ps_single(s);
    lu0   = _mm256_cvtps_pd(_mm256_castps256_ps128(ir));
    lu1   = _mm256_cvtps_pd(_mm256_extractf128_ps(ir, 1));
    *out0 = _mm256_mul_pd(half, _mm256_mul_pd(_mm256_sub_pd(three, _mm256_mul_pd(_mm256_mul_pd(lu0, lu0), in0)), lu0));
    *out1 = _mm256_mul_pd(half, _mm256_mul_pd(_mm256_sub_pd(three, _mm256_mul_pd(_mm256_mul_pd(lu1, lu1), in1)), lu1));
}

static gmx_inline void
load_lj_pair_params(const real *nbfp, const int *type, int aj,
                    __m256d *c6_S, __m256d *c12_S)
{
    __m128d clj_S[UNROLLJ], c6t_S[2], c12t_S[2];
    int     p;

    for (p = 0; p < UNROLLJ; p++)
    {
        clj_S[p] = _mm_load_pd(nbfp+type[aj+p]*nbfp_stride);
    }
    gmx_mm_transpose2_op_pd(clj_S[0], clj_S[1], &c6t_S[0], &c12t_S[0]);
    gmx_mm_transpose2_op_pd(clj_S[2], clj_S[3], &c6t_S[1], &c12t_S[1]);
    *c6_S  = gmx_2_m128d_to_m256d(c6t_S[0], c6t_S[1]);
    *c12_S = gmx_2_m128d_to_m256d(c12t_S[0], c12t_S[1]);
}

/* The load_table functions below are performance critical. They
 * always take the ti parameter, which should contain a buffer that
 * is aligned with prepare_table_load_buffer(), but it is only used
 * with full-width AVX_256. */

static gmx_inline void
load_table_f(const real *tab_coul_F, __m128i ti_S, int *ti,
             __m256d *ctab0_S, __m256d *ctab1_S)
{
    __m128d ctab_S[4], tr_S[4];
    int     j;

    _mm_store_si128((__m128i *)ti, ti_S);
    for (j = 0; j < 4; j++)
    {
        ctab_S[j] = _mm_loadu_pd(tab_coul_F+ti[j]);
    }
    /* Shuffle the force table entries to a convenient order */
    gmx_mm_transpose2_op_pd(ctab_S[0], ctab_S[1], &tr_S[0], &tr_S[1]);
    gmx_mm_transpose2_op_pd(ctab_S[2], ctab_S[3], &tr_S[2], &tr_S[3]);
    *ctab0_S = gmx_2_m128d_to_m256d(tr_S[0], tr_S[2]);
    *ctab1_S = gmx_2_m128d_to_m256d(tr_S[1], tr_S[3]);
    /* The second force table entry should contain the difference */
    *ctab1_S = _mm256_sub_pd(*ctab1_S, *ctab0_S);
}

static gmx_inline void
load_table_f_v(const real *tab_coul_F, const real *tab_coul_V,
               __m128i ti_S, int *ti,
               __m256d *ctab0_S, __m256d *ctab1_S, __m256d *ctabv_S)
{
    __m128d ctab_S[8], tr_S[4];
    int     j;

    _mm_store_si128((__m128i *)ti, ti_S);
    for (j = 0; j < 4; j++)
    {
        ctab_S[j] = _mm_loadu_pd(tab_coul_F+ti[j]);
    }
    /* Shuffle the force table entries to a convenient order */
    gmx_mm_transpose2_op_pd(ctab_S[0], ctab_S[1], &tr_S[0], &tr_S[1]);
    gmx_mm_transpose2_op_pd(ctab_S[2], ctab_S[3], &tr_S[2], &tr_S[3]);
    *ctab0_S = gmx_2_m128d_to_m256d(tr_S[0], tr_S[2]);
    *ctab1_S = gmx_2_m128d_to_m256d(tr_S[1], tr_S[3]);
    /* The second force table entry should contain the difference */
    *ctab1_S = _mm256_sub_pd(*ctab1_S, *ctab0_S);

    for (j = 0; j < 4; j++)
    {
        ctab_S[4+j] = _mm_loadu_pd(tab_coul_V+ti[j]);
    }
    /* Shuffle the energy table entries to a single register */
    *ctabv_S = gmx_2_m128d_to_m256d(_mm_shuffle_pd(ctab_S[4], ctab_S[5], _MM_SHUFFLE2(0, 0)), _mm_shuffle_pd(ctab_S[6], ctab_S[7], _MM_SHUFFLE2(0, 0)));
}

static gmx_inline gmx_exclfilter
gmx_load1_exclfilter(int e)
{
    return _mm256_castsi256_pd(_mm256_set1_epi32(e));
}

static gmx_inline gmx_exclfilter
gmx_load_exclusion_filter(const unsigned *i)
{
    return gmx_simd_load_r((real *) (i));
}

static gmx_inline gmx_simd_bool_t
gmx_checkbitmask_pb(gmx_exclfilter m0, gmx_exclfilter m1)
{
    /* With <= 16 bits used the cast and conversion should not be
     * required, since only mantissa bits are set and that would give
     * a non-zero float, but with the Intel compiler this does not
     * work correctly. Because AVX does not have int->double
     * conversion, we convert via float. */
    return _mm256_cmp_pd(_mm256_castps_pd(_mm256_cvtepi32_ps(_mm256_castpd_si256(_mm256_and_pd(m0, m1)))), _mm256_setzero_pd(), 0x0c);
}

#endif /* _nbnxn_kernel_simd_utils_x86_s256d_h_ */
