/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_X86_SSE4_1_H
#define GMX_SIMD_IMPL_X86_SSE4_1_H

#include "config.h"

#include <math.h>

#include <smmintrin.h>

/* x86 SSE4.1 SIMD instruction wrappers
 *
 * Please see documentation in gromacs/simd/simd.h for the available
 * defines.
 */

/* Inherit most of SSE4.1 from SSE2 */
#include "gromacs/simd/impl_x86_sse2/impl_x86_sse2.h"
/* Increment over SSE2 capabilities */
#define GMX_SIMD_X86_SSE4_1_OR_HIGHER


/* Override capability definitions from SSE2 */
#define  GMX_SIMD4_HAVE_FLOAT_DOTPRODUCT3

/* Almost all SSE4.1 instructions already exist in SSE2, but a few of them
 * can be implemented more efficiently in SSE4.1.
 */

/* SINGLE */
#undef  gmx_simd_round_f
#define gmx_simd_round_f(x)                 _mm_round_ps(x, _MM_FROUND_NINT)
#undef  gmx_simd_trunc_f
#define gmx_simd_trunc_f(x)                 _mm_round_ps(x, _MM_FROUND_TRUNC)
#undef  gmx_simd_extract_fi
#define gmx_simd_extract_fi                 _mm_extract_epi32
#undef  gmx_simd_mul_fi
#define gmx_simd_mul_fi                     _mm_mullo_epi32
#undef  gmx_simd_blendv_f
#define gmx_simd_blendv_f                   _mm_blendv_ps
#undef  gmx_simd_reduce_f
#define gmx_simd_reduce_f(a)                gmx_simd_reduce_f_sse4_1(a)
#undef  gmx_simd_blendv_fi
#define gmx_simd_blendv_fi                  _mm_blendv_epi8
#undef  gmx_simd_gather_load_bysimdint_transpose_f
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_sse4_1
#undef  gmx_simd_reduce_incr_4_return_sum_f
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_sse4_1

/* DOUBLE */
#undef  gmx_simd_round_d
#define gmx_simd_round_d(x)                 _mm_round_pd(x, _MM_FROUND_NINT)
#undef  gmx_simd_trunc_d
#define gmx_simd_trunc_d(x)                 _mm_round_pd(x, _MM_FROUND_TRUNC)
#undef  gmx_simd_extract_di
#define gmx_simd_extract_di                 _mm_extract_epi32
#undef  gmx_simd_mul_di
#define gmx_simd_mul_di                     _mm_mullo_epi32
#undef  gmx_simd_blendv_d
#define gmx_simd_blendv_d                   _mm_blendv_pd
#undef  gmx_simd_reduce_d
#define gmx_simd_reduce_d(a)                gmx_simd_reduce_d_sse4_1(a)
#undef  gmx_simd_blendv_di
#define gmx_simd_blendv_di                  _mm_blendv_epi8
#undef  gmx_simd_gather_load_bysimdint_transpose_d
#define gmx_simd_gather_load_bysimdint_transpose_d   gmx_simd_gather_load_bysimdint_transpose_d_sse4_1
#undef  gmx_simd_reduce_incr_4_return_sum_d
#define gmx_simd_reduce_incr_4_return_sum_d          gmx_simd_reduce_incr_4_return_sum_d_sse4_1

/* We only need to override the debugging versions of rsqrt/rcp. Sorry for the double-negative check. */
#ifndef NDEBUG
#    undef  gmx_simd_rcp_mask_f
#    define gmx_simd_rcp_mask_f(a, m)    _mm_and_ps(_mm_rcp_ps(_mm_blendv_ps(_mm_set1_ps(1.0f), a, m)), m)
#    undef  gmx_simd_rsqrt_mask_f
#    define gmx_simd_rsqrt_mask_f(a, m)  _mm_and_ps(_mm_rsqrt_ps(_mm_blendv_ps(_mm_set1_ps(1.0f), a, m)), m)
#    undef  gmx_simd_rcp_mask_d
#    define gmx_simd_rcp_mask_d(a, m)    _mm_and_pd(gmx_simd_rcp_d(_mm_blendv_pd(_mm_set1_pd(1.0), a, m)), m)
#    undef  gmx_simd_rsqrt_mask_d
#    define gmx_simd_rsqrt_mask_d(a, m)  _mm_and_pd(gmx_simd_rsqrt_d(_mm_blendv_pd(_mm_set1_pd(1.0), a, m)), m)
#endif

/* SIMD4 */
#undef  gmx_simd4_dotproduct3_f
#define gmx_simd4_dotproduct3_f         gmx_simd4_dotproduct3_f_sse4_1


/* SIMD reduction function */
static gmx_inline float gmx_simdcall
gmx_simd_reduce_f_sse4_1(__m128 a)
{
    float  f;
    /* Shuffle has latency 1/throughput 1, followed by add with latency 3, t-put 1.
     * This is likely faster than using _mm_hadd_ps, which has latency 5, t-put 2.
     */
    a = _mm_add_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(1, 0, 3, 2)));
    a = _mm_add_ss(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(0, 3, 2, 1)));
    /* This macro usually provides slightly better performance than _mm_cvtss_f32(). */
    _MM_EXTRACT_FLOAT(f, a, 0);
    return f;
}

/* SIMD4 Dotproduct helper function */
static gmx_inline float gmx_simdcall
gmx_simd4_dotproduct3_f_sse4_1(__m128 a, __m128 b)
{
    float f;
    _MM_EXTRACT_FLOAT(f, _mm_dp_ps(a, b, 0x71), 0);
    return f;
}

static gmx_inline double gmx_simdcall
gmx_simd_reduce_d_sse4_1(__m128d a)
{
    double  f;
    /* Shuffle has latency 1/throughput 1, followed by add with latency 3, t-put 1.
     * This is likely faster than using _mm_hadd_ps, which has latency 5, t-put 2.
     */
    a = _mm_add_sd(a, _mm_shuffle_pd(a, a, _MM_SHUFFLE2(1, 1)));
    _mm_store_sd(&f, a);
    return f;
}


/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus
template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_sse4_1(const float *       base,
                                                  gmx_simd_fint32_t   offset,
                                                  gmx_simd_float_t   &v0,
                                                  gmx_simd_float_t   &v1,
                                                  gmx_simd_float_t   &v2,
                                                  gmx_simd_float_t   &v3)
{
    /* Use optimized bit-shift multiply for the most common alignments */
    if (align == 4)
    {
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 8)
    {
        offset = _mm_slli_epi32(offset, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset = _mm_add_epi32(offset, _mm_slli_epi32(offset, 1));
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 16)
    {
        offset = _mm_slli_epi32(offset, 4);
    }

    if (align == 4 || align == 8 || align == 12 || align == 16)
    {
        v0  = _mm_load_ps(base + _mm_extract_epi32(offset, 0) );
        v1  = _mm_load_ps(base + _mm_extract_epi32(offset, 1) );
        v2  = _mm_load_ps(base + _mm_extract_epi32(offset, 2) );
        v3  = _mm_load_ps(base + _mm_extract_epi32(offset, 3) );
    }
    else
    {
        v0  = _mm_load_ps(base + align * _mm_extract_epi32(offset, 0) );
        v1  = _mm_load_ps(base + align * _mm_extract_epi32(offset, 1) );
        v2  = _mm_load_ps(base + align * _mm_extract_epi32(offset, 2) );
        v3  = _mm_load_ps(base + align * _mm_extract_epi32(offset, 3) );
    }
    _MM_TRANSPOSE4_PS(v0, v1, v2, v3);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_sse4_1(const float *       base,
                                                  gmx_simd_fint32_t   offset,
                                                  gmx_simd_float_t   &v0,
                                                  gmx_simd_float_t   &v1)
{
    __m128 t1, t2;

    /* Use optimized bit-shift multiply for the most common alignments */
    if (align == 2)
    {
        offset = _mm_slli_epi32(offset, 1);
    }
    else if (align == 4)
    {
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 6)
    {
        /* multiply by 3, then by 2 */
        offset = _mm_add_epi32(offset, _mm_slli_epi32(offset, 1));
        offset = _mm_slli_epi32(offset, 1);
    }
    else if (align == 8)
    {
        offset = _mm_slli_epi32(offset, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset = _mm_add_epi32(offset, _mm_slli_epi32(offset, 1));
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 16)
    {
        offset = _mm_slli_epi32(offset, 4);
    }

    if (align == 2 || align == 4 || align == 6 ||
        align == 8 || align == 12 || align == 16)
    {
        v0  = _mm_castpd_ps(_mm_load_sd((const double *)(base + _mm_extract_epi32(offset, 0))));
        v1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + _mm_extract_epi32(offset, 1))));
        t1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + _mm_extract_epi32(offset, 2))));
        t2  = _mm_castpd_ps(_mm_load_sd((const double *)(base + _mm_extract_epi32(offset, 3))));
    }
    else
    {
        v0  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * _mm_extract_epi32(offset, 0))));
        v1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * _mm_extract_epi32(offset, 1))));
        t1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * _mm_extract_epi32(offset, 2))));
        t2  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * _mm_extract_epi32(offset, 3))));
    }
    t1  = _mm_unpacklo_ps(v0, t1);
    t2  = _mm_unpacklo_ps(v1, t2);
    v0  = _mm_unpacklo_ps(t1, t2);
    v1  = _mm_unpackhi_ps(t1, t2);
}

static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_f_sse4_1(float *           m,
                                           gmx_simd_float_t  v0,
                                           gmx_simd_float_t  v1,
                                           gmx_simd_float_t  v2,
                                           gmx_simd_float_t  v3)
{
    /* The present implementation is quite standardized: transpose and sum.
     * However, there are a couple of alternatives worth trying and benchmarking later:
     *
     * 1) Use 3*hadd. Much fewer instructions, but high latency (5) and bad throughput (2).
     * 2) Mix hadd with shuffle/blend (blend is better t-put than shuffle).
     *
     * In a standalone test all these had similar performance, so we need to test what
     * works best in the actual kernels where it is mixed with other instructions.
     */
    _MM_TRANSPOSE4_PS(v0, v1, v2, v3);
    v0 = _mm_add_ps(v0, v1);
    v2 = _mm_add_ps(v2, v3);
    v0 = _mm_add_ps(v0, v2);
    v2 = _mm_add_ps(v0, _mm_load_ps(m));
    _mm_store_ps(m, v2);

    return gmx_simd_reduce_f_sse4_1(v0);
}
#endif


/****************************************************
 * Double precision higher-level utility functions  *
 ****************************************************/



#ifdef __cplusplus
template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_sse4_1(const double *       base,
                                                  gmx_simd_dint32_t    offset,
                                                  gmx_simd_double_t   &v0,
                                                  gmx_simd_double_t   &v1,
                                                  gmx_simd_double_t   &v2,
                                                  gmx_simd_double_t   &v3)
{
    __m128d t1, t2, t3, t4;
    /* Use optimized bit-shift multiply for the most common alignments */
    if (align == 4)
    {
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 8)
    {
        offset = _mm_slli_epi32(offset, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset = _mm_add_epi32(offset, _mm_slli_epi32(offset, 1));
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 16)
    {
        offset = _mm_slli_epi32(offset, 4);
    }

    if (align == 4 || align == 8 || align == 12 || align == 16)
    {
        t1  = _mm_load_pd(base + _mm_extract_epi32(offset, 0));
        t2  = _mm_load_pd(base + _mm_extract_epi32(offset, 1));
        t3  = _mm_load_pd(base + _mm_extract_epi32(offset, 0) + 2);
        t4  = _mm_load_pd(base + _mm_extract_epi32(offset, 1) + 2);
    }
    else
    {
        t1  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 0));
        t2  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 1));
        t3  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 0) + 2);
        t4  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 1) + 2);
    }
    v0  = _mm_unpacklo_pd(t1, t2);
    v1  = _mm_unpackhi_pd(t1, t2);
    v2  = _mm_unpacklo_pd(t3, t4);
    v3  = _mm_unpackhi_pd(t3, t4);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_sse4_1(const double *       base,
                                                  gmx_simd_dint32_t    offset,
                                                  gmx_simd_double_t   &v0,
                                                  gmx_simd_double_t   &v1)
{
    __m128d t1, t2;
    /* Use optimized bit-shift multiply for the most common alignments */
    if (align == 2)
    {
        offset = _mm_slli_epi32(offset, 1);
    }
    else if (align == 4)
    {
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 6)
    {
        /* multiply by 3, then by 2 */
        offset = _mm_add_epi32(offset, _mm_slli_epi32(offset, 1));
        offset = _mm_slli_epi32(offset, 1);
    }
    else if (align == 8)
    {
        offset = _mm_slli_epi32(offset, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset = _mm_add_epi32(offset, _mm_slli_epi32(offset, 1));
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 16)
    {
        offset = _mm_slli_epi32(offset, 4);
    }

    if (align == 2 || align == 4 || align == 6 ||
        align == 8 || align == 12 || align == 16)
    {
        t1  = _mm_load_pd(base + _mm_extract_epi32(offset, 0));
        t2  = _mm_load_pd(base + _mm_extract_epi32(offset, 1));
    }
    else
    {
        t1  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 0));
        t2  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 1));
    }
    v0  = _mm_unpacklo_pd(t1, t2);
    v1  = _mm_unpackhi_pd(t1, t2);
}

static gmx_inline double
gmx_simd_reduce_incr_4_return_sum_d_sse4_1(double *           m,
                                           gmx_simd_double_t  v0,
                                           gmx_simd_double_t  v1,
                                           gmx_simd_double_t  v2,
                                           gmx_simd_double_t  v3)
{
    __m128d t1, t2, t3, t4;

    t1 = _mm_unpacklo_pd(v0, v1);
    t2 = _mm_unpackhi_pd(v0, v1);
    t3 = _mm_unpacklo_pd(v2, v3);
    t4 = _mm_unpackhi_pd(v2, v3);

    t1 = _mm_add_pd(t1, t2);
    t3 = _mm_add_pd(t3, t4);

    t2 = _mm_add_pd(t1, _mm_load_pd(m));
    t4 = _mm_add_pd(t3, _mm_load_pd(m + 2));
    _mm_store_pd(m, t2);
    _mm_store_pd(m + 2, t4);

    t1 = _mm_add_pd(t1, t3);
    return gmx_simd_reduce_d_sse4_1(t1);
}
#endif

#endif /* GMX_SIMD_IMPL_X86_SSE4_1_H */
