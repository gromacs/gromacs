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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_256_UTIL_FLOAT_H

#include "config.h"

#include <immintrin.h>

#include "impl_x86_avx_256_common.h"

/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_f             gmx_simd_gather_load_transpose_f_avx_256
#define gmx_simd_best_pair_alignment_f               gmx_simd_best_pair_alignment_f_avx_256
#define gmx_simd_gather_loadu_transpose_f            gmx_simd_gather_loadu_transpose_f_avx_256
#define gmx_simd_transpose_scatter_storeu_f          gmx_simd_transpose_scatter_storeu_f_avx_256
#define gmx_simd_transpose_scatter_incru_f           gmx_simd_transpose_scatter_incru_f_avx_256
#define gmx_simd_transpose_scatter_decru_f           gmx_simd_transpose_scatter_decru_f_avx_256
#define gmx_simd_expand_scalars_to_triplets_f        gmx_simd_expand_scalars_to_triplets_f_avx_256
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_avx_256
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_avx_256
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_avx_256
/* Half-simd-width utilities (only needed in single, where width is 8) */
#define gmx_simd_load_dual_hsimd_f                   gmx_simd_load_dual_hsimd_f_avx_256
#define gmx_simd_loaddup_hsimd_f                     gmx_simd_loaddup_hsimd_f_avx_256
#define gmx_simd_load1_dual_hsimd_f                  gmx_simd_load1_dual_hsimd_f_avx_256
#define gmx_simd_store_dual_hsimd_f                  gmx_simd_store_dual_hsimd_f_avx_256
#define gmx_simd_decr_hsimd_f                        gmx_simd_decr_hsimd_f_avx_256
#define gmx_simd_gather_load_transpose_hsimd_f       gmx_simd_gather_load_transpose_hsimd_f_avx_256
#define gmx_simd_reduce_incr_4_return_sum_hsimd_f    gmx_simd_reduce_incr_4_return_sum_hsimd_f_avx_256

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

/*
 * This is an internal helper macro used by the three functions storing, incrementing, or
 * decrementing data. Do NOT use it outside this file.
 *
 * Input v0: [x0 x1 x2 x3 x4 x5 x6 x7]
 * Input v1: [y0 y1 y2 y3 y4 y5 y6 y7]
 * Input v2: [z0 z1 z2 z3 z4 z5 z6 z7]
 * Input v3: Unused
 *
 * Output v0: [x0 y0 z0 -  x4 y4 z4 - ]
 * Output v1: [x1 y1 z1 -  x5 y5 z5 - ]
 * Output v2: [x2 y2 z2 -  x6 y6 z6 - ]
 * Output v3: [x3 y3 z3 -  x7 y7 z7 - ]
 *
 * Here, - means undefined. Note that such values will not be zero!
 */
#define GMX_MM_TRANSPOSE_3TO4X2_FLOAT(v0, v1, v2, v3)                  \
    {                                                                   \
        __m256 gmx_mm_t1 = _mm256_unpacklo_ps(v0, v1);                  \
        __m256 gmx_mm_t2 = _mm256_unpackhi_ps(v0, v1);                  \
        v0 = _mm256_shuffle_ps(gmx_mm_t1, v2, _MM_SHUFFLE(0, 0, 1, 0)); \
        v1 = _mm256_shuffle_ps(gmx_mm_t1, v2, _MM_SHUFFLE(0, 1, 3, 2)); \
        v3 = _mm256_shuffle_ps(gmx_mm_t2, v2, _MM_SHUFFLE(0, 3, 3, 2)); \
        v2 = _mm256_shuffle_ps(gmx_mm_t2, v2, _MM_SHUFFLE(0, 2, 1, 0)); \
    }

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_transpose_f_avx_256(const float *        base,
                                         const gmx_int32_t    offset[],
                                         gmx_simd_float_t    &v0,
                                         gmx_simd_float_t    &v1,
                                         gmx_simd_float_t    &v2,
                                         gmx_simd_float_t    &v3)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;

    t1  = _mm_load_ps( base + align * offset[0] );
    t2  = _mm_load_ps( base + align * offset[1] );
    t3  = _mm_load_ps( base + align * offset[2] );
    t4  = _mm_load_ps( base + align * offset[3] );
    t5  = _mm_load_ps( base + align * offset[4] );
    t6  = _mm_load_ps( base + align * offset[5] );
    t7  = _mm_load_ps( base + align * offset[6] );
    t8  = _mm_load_ps( base + align * offset[7] );

    v0  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    v1  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    v2  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    v3  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(v0, v1);
    tB  = _mm256_unpacklo_ps(v2, v3);
    tC  = _mm256_unpackhi_ps(v0, v1);
    tD  = _mm256_unpackhi_ps(v2, v3);

    v0  = _mm256_shuffle_ps(tA, tB, _MM_SHUFFLE(1, 0, 1, 0));
    v1  = _mm256_shuffle_ps(tA, tB, _MM_SHUFFLE(3, 2, 3, 2));
    v2  = _mm256_shuffle_ps(tC, tD, _MM_SHUFFLE(1, 0, 1, 0));
    v3  = _mm256_shuffle_ps(tC, tD, _MM_SHUFFLE(3, 2, 3, 2));
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_transpose_f_avx_256(const float *        base,
                                         const gmx_int32_t    offset[],
                                         gmx_simd_float_t    &v0,
                                         gmx_simd_float_t    &v1)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;

    t1  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[0] ) );
    t2  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[1] ) );
    t3  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[2] ) );
    t4  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[3] ) );
    t5  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[4] ) );
    t6  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[5] ) );
    t7  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[6] ) );
    t8  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[7] ) );

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(tA, tC);
    tB  = _mm256_unpacklo_ps(tB, tD);
    v0  = _mm256_unpacklo_ps(tA, tB);
    v1  = _mm256_unpackhi_ps(tA, tB);
}

static const int gmx_simd_best_pair_alignment_f_avx_256 = 2;

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_loadu_transpose_f_avx_256(const float *        base,
                                          const gmx_int32_t    offset[],
                                          gmx_simd_float_t    &v0,
                                          gmx_simd_float_t    &v1,
                                          gmx_simd_float_t    &v2)
{
    __m256  t1, t2, t3, t4, t5, t6, t7, t8;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    if ( (align & 0x3) == 0)
    {
        /* With alignment 4 or better we can read a byte beyond triplets and use _mm_loadu_ps() */
        t1  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps( base + align * offset[0] )),
                                   _mm_loadu_ps( base + align * offset[4] ), 0x1);
        t2  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps(base + align * offset[1] )),
                                   _mm_loadu_ps( base + align * offset[5] ), 0x1);
        t3  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps(base + align * offset[2] )),
                                   _mm_loadu_ps( base + align * offset[6] ), 0x1);
        t4  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps(base + align * offset[3] )),
                                   _mm_loadu_ps( base + align * offset[7] ), 0x1);
    }
    else
    {
        /* Arbitrary alignment */
        t1  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[0], mask)),
                                   gmx_mm_maskload_ps( base + align * offset[4], mask), 0x1);
        t2  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[1], mask)),
                                   gmx_mm_maskload_ps( base + align * offset[5], mask), 0x1);
        t3  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[2], mask)),
                                   gmx_mm_maskload_ps( base + align * offset[6], mask), 0x1);
        t4  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[3], mask)),
                                   gmx_mm_maskload_ps( base + align * offset[7], mask), 0x1);
    }
    t5  = _mm256_unpacklo_ps(t1, t2);
    t6  = _mm256_unpacklo_ps(t3, t4);
    t7  = _mm256_unpackhi_ps(t1, t2);
    t8  = _mm256_unpackhi_ps(t3, t4);
    v0  = _mm256_shuffle_ps(t5, t6, _MM_SHUFFLE(1, 0, 1, 0));
    v1  = _mm256_shuffle_ps(t5, t6, _MM_SHUFFLE(3, 2, 3, 2));
    v2  = _mm256_shuffle_ps(t7, t8, _MM_SHUFFLE(1, 0, 1, 0));
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_storeu_f_avx_256(float *              base,
                                            const gmx_int32_t    offset[],
                                            gmx_simd_float_t     v0,
                                            gmx_simd_float_t     v1,
                                            gmx_simd_float_t     v2)
{
    __m256  v3;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    GMX_MM_TRANSPOSE_3TO4X2_FLOAT(v0, v1, v2, v3);
    gmx_mm_maskstore_ps( base + align * offset[0], mask, _mm256_castps256_ps128(v0));
    gmx_mm_maskstore_ps( base + align * offset[1], mask, _mm256_castps256_ps128(v1));
    gmx_mm_maskstore_ps( base + align * offset[2], mask, _mm256_castps256_ps128(v2));
    gmx_mm_maskstore_ps( base + align * offset[3], mask, _mm256_castps256_ps128(v3));
    gmx_mm_maskstore_ps( base + align * offset[4], mask, _mm256_extractf128_ps(v0, 0x1));
    gmx_mm_maskstore_ps( base + align * offset[5], mask, _mm256_extractf128_ps(v1, 0x1));
    gmx_mm_maskstore_ps( base + align * offset[6], mask, _mm256_extractf128_ps(v2, 0x1));
    gmx_mm_maskstore_ps( base + align * offset[7], mask, _mm256_extractf128_ps(v3, 0x1));
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_incru_f_avx_256(float *              base,
                                           const gmx_int32_t    offset[],
                                           gmx_simd_float_t     v0,
                                           gmx_simd_float_t     v1,
                                           gmx_simd_float_t     v2)
{
    __m256  v3, t0, t1, t2, t3;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    t0  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[0], mask)),
                               gmx_mm_maskload_ps( base + align * offset[4], mask), 0x1);
    t1  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[1], mask)),
                               gmx_mm_maskload_ps( base + align * offset[5], mask), 0x1);
    t2  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[2], mask)),
                               gmx_mm_maskload_ps( base + align * offset[6], mask), 0x1);
    t3  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[3], mask)),
                               gmx_mm_maskload_ps( base + align * offset[7], mask), 0x1);

    GMX_MM_TRANSPOSE_3TO4X2_FLOAT(v0, v1, v2, v3);
    t0          = _mm256_add_ps(t0, v0);
    t1          = _mm256_add_ps(t1, v1);
    t2          = _mm256_add_ps(t2, v2);
    t3          = _mm256_add_ps(t3, v3);

    gmx_mm_maskstore_ps(base + align * offset[0], mask, _mm256_castps256_ps128(t0));
    gmx_mm_maskstore_ps(base + align * offset[1], mask, _mm256_castps256_ps128(t1));
    gmx_mm_maskstore_ps(base + align * offset[2], mask, _mm256_castps256_ps128(t2));
    gmx_mm_maskstore_ps(base + align * offset[3], mask, _mm256_castps256_ps128(t3));
    gmx_mm_maskstore_ps(base + align * offset[4], mask, _mm256_extractf128_ps(t0, 0x1));
    gmx_mm_maskstore_ps(base + align * offset[5], mask, _mm256_extractf128_ps(t1, 0x1));
    gmx_mm_maskstore_ps(base + align * offset[6], mask, _mm256_extractf128_ps(t2, 0x1));
    gmx_mm_maskstore_ps(base + align * offset[7], mask, _mm256_extractf128_ps(t3, 0x1));
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_decru_f_avx_256(float *              base,
                                           const gmx_int32_t    offset[],
                                           gmx_simd_float_t     v0,
                                           gmx_simd_float_t     v1,
                                           gmx_simd_float_t     v2)
{
    __m256  v3, t0, t1, t2, t3;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    t0  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[0], mask)),
                               gmx_mm_maskload_ps( base + align * offset[4], mask), 0x1);
    t1  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[1], mask)),
                               gmx_mm_maskload_ps( base + align * offset[5], mask), 0x1);
    t2  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[2], mask)),
                               gmx_mm_maskload_ps( base + align * offset[6], mask), 0x1);
    t3  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[3], mask)),
                               gmx_mm_maskload_ps( base + align * offset[7], mask), 0x1);

    GMX_MM_TRANSPOSE_3TO4X2_FLOAT(v0, v1, v2, v3);
    t0          = _mm256_sub_ps(t0, v0);
    t1          = _mm256_sub_ps(t1, v1);
    t2          = _mm256_sub_ps(t2, v2);
    t3          = _mm256_sub_ps(t3, v3);

    gmx_mm_maskstore_ps(base + align * offset[0], mask, _mm256_castps256_ps128(t0));
    gmx_mm_maskstore_ps(base + align * offset[1], mask, _mm256_castps256_ps128(t1));
    gmx_mm_maskstore_ps(base + align * offset[2], mask, _mm256_castps256_ps128(t2));
    gmx_mm_maskstore_ps(base + align * offset[3], mask, _mm256_castps256_ps128(t3));
    gmx_mm_maskstore_ps(base + align * offset[4], mask, _mm256_extractf128_ps(t0, 0x1));
    gmx_mm_maskstore_ps(base + align * offset[5], mask, _mm256_extractf128_ps(t1, 0x1));
    gmx_mm_maskstore_ps(base + align * offset[6], mask, _mm256_extractf128_ps(t2, 0x1));
    gmx_mm_maskstore_ps(base + align * offset[7], mask, _mm256_extractf128_ps(t3, 0x1));
}

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_f_avx_256(gmx_simd_float_t   scalar,
                                              gmx_simd_float_t  &triplets0,
                                              gmx_simd_float_t  &triplets1,
                                              gmx_simd_float_t  &triplets2)
{
    __m256 t0 = _mm256_permute2f128_ps(scalar, scalar, 0x21);
    __m256 t1 = _mm256_permute_ps(scalar, _MM_SHUFFLE(1, 0, 0, 0));
    __m256 t2 = _mm256_permute_ps(t0, _MM_SHUFFLE(2, 2, 1, 1));
    __m256 t3 = _mm256_permute_ps(scalar, _MM_SHUFFLE(3, 3, 3, 2));
    triplets0 = _mm256_blend_ps(t1, t2, 0xF0);
    triplets1 = _mm256_blend_ps(t3, t1, 0xF0);
    triplets2 = _mm256_blend_ps(t2, t3, 0xF0);
}


template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_f_avx_256(const float *       base,
                                                   gmx_simd_fint32_t   simdoffset,
                                                   gmx_simd_float_t   &v0,
                                                   gmx_simd_float_t   &v1,
                                                   gmx_simd_float_t   &v2,
                                                   gmx_simd_float_t   &v3)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;
    /* All x86 compilers we know that support SIMD also support GCC alignment attributes */
    int __attribute__ ((aligned (32))) offset[GMX_SIMD_FLOAT_WIDTH];

    /* Plain AVX does not support 256-bit integer shift operations,
     * so we always store to memory and multiply in the integer units instead.
     */
    _mm256_store_si256( (__m256i *)offset, simdoffset);

    t1  = _mm_load_ps( base + align * offset[0] );
    t2  = _mm_load_ps( base + align * offset[1] );
    t3  = _mm_load_ps( base + align * offset[2] );
    t4  = _mm_load_ps( base + align * offset[3] );
    t5  = _mm_load_ps( base + align * offset[4] );
    t6  = _mm_load_ps( base + align * offset[5] );
    t7  = _mm_load_ps( base + align * offset[6] );
    t8  = _mm_load_ps( base + align * offset[7] );

    v0  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    v1  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    v2  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    v3  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(v0, v1);
    tB  = _mm256_unpacklo_ps(v2, v3);
    tC  = _mm256_unpackhi_ps(v0, v1);
    tD  = _mm256_unpackhi_ps(v2, v3);

    v0  = _mm256_shuffle_ps(tA, tB, _MM_SHUFFLE(1, 0, 1, 0));
    v1  = _mm256_shuffle_ps(tA, tB, _MM_SHUFFLE(3, 2, 3, 2));
    v2  = _mm256_shuffle_ps(tC, tD, _MM_SHUFFLE(1, 0, 1, 0));
    v3  = _mm256_shuffle_ps(tC, tD, _MM_SHUFFLE(3, 2, 3, 2));
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_f_avx_256(const float *       base,
                                                   gmx_simd_fint32_t   simdoffset,
                                                   gmx_simd_float_t   &v0,
                                                   gmx_simd_float_t   &v1)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;
    /* All x86 compilers we know that support SIMD also support GCC alignment attributes */
    int __attribute__ ((aligned (32))) offset[GMX_SIMD_FLOAT_WIDTH];

    /* Plain AVX does not support 256-bit integer shift operations,
     * so we always store to memory and multiply in the integer units instead.
     */
    _mm256_store_si256( (__m256i *)offset, simdoffset);

    t1  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[0] ) );
    t2  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[1] ) );
    t3  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[2] ) );
    t4  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[3] ) );
    t5  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[4] ) );
    t6  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[5] ) );
    t7  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[6] ) );
    t8  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[7] ) );

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(tA, tC);
    tB  = _mm256_unpacklo_ps(tB, tD);
    v0  = _mm256_unpacklo_ps(tA, tB);
    v1  = _mm256_unpackhi_ps(tA, tB);
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_loadu_bysimdint_transpose_f_avx_256(const float *       base,
                                                    gmx_simd_fint32_t   simdoffset,
                                                    gmx_simd_float_t   &v0,
                                                    gmx_simd_float_t   &v1)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;
    /* All x86 compilers we know that support SIMD also support GCC alignment attributes */
    int __attribute__ ((aligned (32))) offset[GMX_SIMD_FLOAT_WIDTH];

    /* Plain AVX does not support 256-bit integer shift operations,
     * so we always store to memory and multiply in the integer units instead.
     */
    _mm256_store_si256( (__m256i *)offset, simdoffset);

    t1  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[0] ) );
    t2  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[1] ) );
    t3  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[2] ) );
    t4  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[3] ) );
    t5  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[4] ) );
    t6  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[5] ) );
    t7  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[6] ) );
    t8  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[7] ) );

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(tA, tC);
    tB  = _mm256_unpacklo_ps(tB, tD);
    v0  = _mm256_unpacklo_ps(tA, tB);
    v1  = _mm256_unpackhi_ps(tA, tB);
}

static gmx_inline float gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_f_avx_256(float *           m,
                                            gmx_simd_float_t  v0,
                                            gmx_simd_float_t  v1,
                                            gmx_simd_float_t  v2,
                                            gmx_simd_float_t  v3)
{
    float  f;
    __m128 t0, t2;
    v0 = _mm256_hadd_ps(v0, v1);
    v2 = _mm256_hadd_ps(v2, v3);
    v0 = _mm256_hadd_ps(v0, v2);
    t0 = _mm_add_ps(_mm256_castps256_ps128(v0), _mm256_extractf128_ps(v0, 0x1));

    t2 = _mm_add_ps(t0, _mm_load_ps(m));
    _mm_store_ps(m, t2);

    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    _MM_EXTRACT_FLOAT(f, t0, 0);
    return f;
}

/******************************************************
 * Single precision half-simd-width utility functions *
 ******************************************************/
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_load_dual_hsimd_f_avx_256(const float * m0,
                                   const float * m1)
{
    return _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_load_ps(m0)), _mm_load_ps(m1), 0x1);
}

static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_loaddup_hsimd_f_avx_256(const float * m)
{
    return _mm256_broadcast_ps((const __m128 *)m);
}

static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_load1_dual_hsimd_f_avx_256(const float * m)
{
    __m128 t0, t1;
    t0 = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)m);
    t1 = _mm_permute_ps(t0, _MM_SHUFFLE(1, 1, 1, 1));
    t0 = _mm_permute_ps(t0, _MM_SHUFFLE(0, 0, 0, 0));
    return _mm256_insertf128_ps(_mm256_castps128_ps256(t0), t1, 0x1);
}


static gmx_inline void gmx_simdcall
gmx_simd_store_dual_hsimd_f_avx_256(float *           m0,
                                    float *           m1,
                                    gmx_simd_float_t  a)
{
    _mm_store_ps(m0, _mm256_castps256_ps128(a));
    _mm_store_ps(m1, _mm256_extractf128_ps(a, 0x1));
}

static gmx_inline void gmx_simdcall
gmx_simd_decr_hsimd_f_avx_256(float *          m,
                              gmx_simd_float_t a)
{
    __m128 asum = _mm_add_ps(_mm256_castps256_ps128(a), _mm256_extractf128_ps(a, 0x1));
    _mm_store_ps(m, _mm_sub_ps(_mm_load_ps(m), asum));
}


template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_transpose_hsimd_f_avx_256(const float *       base0,
                                               const float *       base1,
                                               gmx_int32_t         offset[],
                                               gmx_simd_float_t   &v0,
                                               gmx_simd_float_t   &v1)
{
    __m128 t0, t1, t2, t3, t4, t5, t6, t7;
    __m256 tA, tB, tC, tD;
    t0  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base0 + align * offset[0]));
    t1  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base0 + align * offset[1]));
    t2  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base0 + align * offset[2]));
    t3  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base0 + align * offset[3]));
    t4  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base1 + align * offset[0]));
    t5  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base1 + align * offset[1]));
    t6  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base1 + align * offset[2]));
    t7  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base1 + align * offset[3]));

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t0), t4, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);

    tA  = _mm256_unpacklo_ps(tA, tC);
    tB  = _mm256_unpacklo_ps(tB, tD);
    v0  = _mm256_unpacklo_ps(tA, tB);
    v1  = _mm256_unpackhi_ps(tA, tB);
}


static gmx_inline float gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_hsimd_f_avx_256(float *            m,
                                                  gmx_simd_float_t   v0,
                                                  gmx_simd_float_t   v1)
{
    __m128 t0, t1;
    float  f;

    v0  = _mm256_hadd_ps(v0, v1);
    t0  = _mm256_extractf128_ps(v0, 0x1);
    t0  = _mm_hadd_ps(_mm256_castps256_ps128(v0), t0);
    t0  = _mm_permute_ps(t0, _MM_SHUFFLE(3, 1, 2, 0));

    t1  = _mm_add_ps(t0, _mm_load_ps(m));
    _mm_store_ps(m, t1);

    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    _MM_EXTRACT_FLOAT(f, t0, 0);
    return f;
}
#endif

#endif /* GMX_SIMD_IMPL_X86_AVX_256_UTIL_FLOAT_H */
