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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_256_UTIL_DOUBLE_H

#include "config.h"

#include <immintrin.h>

#include "impl_x86_avx_256_common.h"
#include "impl_x86_avx_256_simd_double.h"

/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_d             gmx_simd_gather_load_transpose_d_avx_256
#define gmx_simd_best_pair_alignment_d               gmx_simd_best_pair_alignment_d_avx_256
#define gmx_simd_gather_loadu_transpose_d            gmx_simd_gather_loadu_transpose_d_avx_256
#define gmx_simd_transpose_scatter_storeu_d          gmx_simd_transpose_scatter_storeu_d_avx_256
#define gmx_simd_transpose_scatter_incru_d           gmx_simd_transpose_scatter_incru_d_avx_256
#define gmx_simd_transpose_scatter_decru_d           gmx_simd_transpose_scatter_decru_d_avx_256
#define gmx_simd_expand_scalars_to_triplets_d        gmx_simd_expand_scalars_to_triplets_d_avx_256
#define gmx_simd_gather_load_bysimdint_transpose_d   gmx_simd_gather_load_bysimdint_transpose_d_avx_256
#define gmx_simd_gather_loadu_bysimdint_transpose_d  gmx_simd_gather_loadu_bysimdint_transpose_d_avx_256
#define gmx_simd_reduce_incr_4_return_sum_d          gmx_simd_reduce_incr_4_return_sum_d_avx_256

/****************************************************
 * Double precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_transpose_d_avx_256(const double *        base,
                                         const gmx_int32_t     offset[],
                                         gmx_simd_double_t    &v0,
                                         gmx_simd_double_t    &v1,
                                         gmx_simd_double_t    &v2,
                                         gmx_simd_double_t    &v3)
{
    v0  = _mm256_load_pd( base + align * offset[0] );
    v1  = _mm256_load_pd( base + align * offset[1] );
    v2  = _mm256_load_pd( base + align * offset[2] );
    v3  = _mm256_load_pd( base + align * offset[3] );
    GMX_MM_TRANSPOSE4_DOUBLE(v0, v1, v2, v3);
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_transpose_d_avx_256(const double *        base,
                                         const gmx_int32_t     offset[],
                                         gmx_simd_double_t    &v0,
                                         gmx_simd_double_t    &v1)
{
    __m128d t1, t2, t3, t4;
    __m256d tA, tB;

    t1   = _mm_load_pd( base + align * offset[0] );
    t2   = _mm_load_pd( base + align * offset[1] );
    t3   = _mm_load_pd( base + align * offset[2] );
    t4   = _mm_load_pd( base + align * offset[3] );
    tA   = _mm256_insertf128_pd(_mm256_castpd128_pd256(t1), t3, 0x1);
    tB   = _mm256_insertf128_pd(_mm256_castpd128_pd256(t2), t4, 0x1);

    v0  = _mm256_unpacklo_pd(tA, tB);
    v1  = _mm256_unpackhi_pd(tA, tB);
}

static const int gmx_simd_best_pair_alignment_d_avx_256 = 2;

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_loadu_transpose_d_avx_256(const double *        base,
                                          const gmx_int32_t     offset[],
                                          gmx_simd_double_t    &v0,
                                          gmx_simd_double_t    &v1,
                                          gmx_simd_double_t    &v2)
{
    __m256i mask = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_set_epi32(-1, -1, -1, -1)), _mm_set_epi32(0, 0, -1, -1), 0x1);
    __m256d t1, t2, t3, t4, t5, t6, t7, t8;
    t1  = _mm256_maskload_pd(base + align * offset[0], mask);
    t2  = _mm256_maskload_pd(base + align * offset[1], mask);
    t3  = _mm256_maskload_pd(base + align * offset[2], mask);
    t4  = _mm256_maskload_pd(base + align * offset[3], mask);
    t5  = _mm256_unpacklo_pd(t1, t2);
    t6  = _mm256_unpackhi_pd(t1, t2);
    t7  = _mm256_unpacklo_pd(t3, t4);
    t8  = _mm256_unpackhi_pd(t3, t4);
    v0  = _mm256_permute2f128_pd(t5, t7, 0x20);
    v1  = _mm256_permute2f128_pd(t6, t8, 0x20);
    v2  = _mm256_permute2f128_pd(t5, t7, 0x31);
}



template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_storeu_d_avx_256(double *              base,
                                            const gmx_int32_t     offset[],
                                            gmx_simd_double_t     v0,
                                            gmx_simd_double_t     v1,
                                            gmx_simd_double_t     v2)
{
    __m256i mask = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_set_epi32(-1, -1, -1, -1)), _mm_set_epi32(0, 0, -1, -1), 0x1);
    __m256d v3   = _mm256_setzero_pd();

    GMX_MM_TRANSPOSE4_DOUBLE(v0, v1, v2, v3);

    _mm256_maskstore_pd(base + align * offset[0], mask, v0);
    _mm256_maskstore_pd(base + align * offset[1], mask, v1);
    _mm256_maskstore_pd(base + align * offset[2], mask, v2);
    _mm256_maskstore_pd(base + align * offset[3], mask, v3);
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_incru_d_avx_256(double *              base,
                                           const gmx_int32_t     offset[],
                                           gmx_simd_double_t     v0,
                                           gmx_simd_double_t     v1,
                                           gmx_simd_double_t     v2)
{
    __m256i mask = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_set_epi32(-1, -1, -1, -1)), _mm_set_epi32(0, 0, -1, -1), 0x1);
    __m256d v3   = _mm256_setzero_pd();
    __m256d t0, t1, t2, t3;

    t0  = _mm256_maskload_pd(base + align * offset[0], mask);
    t1  = _mm256_maskload_pd(base + align * offset[1], mask);
    t2  = _mm256_maskload_pd(base + align * offset[2], mask);
    t3  = _mm256_maskload_pd(base + align * offset[3], mask);

    GMX_MM_TRANSPOSE4_DOUBLE(v0, v1, v2, v3);

    t0  = _mm256_add_pd(t0, v0);
    t1  = _mm256_add_pd(t1, v1);
    t2  = _mm256_add_pd(t2, v2);
    t3  = _mm256_add_pd(t3, v3);

    _mm256_maskstore_pd(base + align * offset[0], mask, t0);
    _mm256_maskstore_pd(base + align * offset[1], mask, t1);
    _mm256_maskstore_pd(base + align * offset[2], mask, t2);
    _mm256_maskstore_pd(base + align * offset[3], mask, t3);
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_decru_d_avx_256(double *              base,
                                           const gmx_int32_t     offset[],
                                           gmx_simd_double_t     v0,
                                           gmx_simd_double_t     v1,
                                           gmx_simd_double_t     v2)
{
    __m256i mask = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_set_epi32(-1, -1, -1, -1)), _mm_set_epi32(0, 0, -1, -1), 0x1);
    __m256d v3   = _mm256_setzero_pd();
    __m256d t0, t1, t2, t3;

    t0  = _mm256_maskload_pd(base + align * offset[0], mask);
    t1  = _mm256_maskload_pd(base + align * offset[1], mask);
    t2  = _mm256_maskload_pd(base + align * offset[2], mask);
    t3  = _mm256_maskload_pd(base + align * offset[3], mask);

    GMX_MM_TRANSPOSE4_DOUBLE(v0, v1, v2, v3);

    t0  = _mm256_sub_pd(t0, v0);
    t1  = _mm256_sub_pd(t1, v1);
    t2  = _mm256_sub_pd(t2, v2);
    t3  = _mm256_sub_pd(t3, v3);

    _mm256_maskstore_pd(base + align * offset[0], mask, t0);
    _mm256_maskstore_pd(base + align * offset[1], mask, t1);
    _mm256_maskstore_pd(base + align * offset[2], mask, t2);
    _mm256_maskstore_pd(base + align * offset[3], mask, t3);
}

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_d_avx_256(gmx_simd_double_t   scalar,
                                              gmx_simd_double_t  &triplets0,
                                              gmx_simd_double_t  &triplets1,
                                              gmx_simd_double_t  &triplets2)
{
    __m256d t0 = _mm256_permute2f128_pd(scalar, scalar, 0x21);
    __m256d t1 = _mm256_permute_pd(scalar, 0b0000);
    __m256d t2 = _mm256_permute_pd(scalar, 0b1111);
    triplets0 = _mm256_blend_pd(t1, t0, 0b1100);
    triplets1 = _mm256_blend_pd(t2, t1, 0b1100);
    triplets2 = _mm256_blend_pd(t0, t2, 0b1100);
}


template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_d_avx_256(const double *       base,
                                                   gmx_simd_dint32_t    offset,
                                                   gmx_simd_double_t   &v0,
                                                   gmx_simd_double_t   &v1,
                                                   gmx_simd_double_t   &v2,
                                                   gmx_simd_double_t   &v3)
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
        v0  = _mm256_load_pd(base + _mm_extract_epi32(offset, 0));
        v1  = _mm256_load_pd(base + _mm_extract_epi32(offset, 1));
        v2  = _mm256_load_pd(base + _mm_extract_epi32(offset, 2));
        v3  = _mm256_load_pd(base + _mm_extract_epi32(offset, 3));
    }
    else
    {
        v0  = _mm256_load_pd(base + align * _mm_extract_epi32(offset, 0));
        v1  = _mm256_load_pd(base + align * _mm_extract_epi32(offset, 1));
        v2  = _mm256_load_pd(base + align * _mm_extract_epi32(offset, 2));
        v3  = _mm256_load_pd(base + align * _mm_extract_epi32(offset, 3));
    }
    GMX_MM_TRANSPOSE4_DOUBLE(v0, v1, v2, v3);
}


template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_d_avx_256(const double *       base,
                                                   gmx_simd_dint32_t    offset,
                                                   gmx_simd_double_t   &v0,
                                                   gmx_simd_double_t   &v1)
{
    __m128d t1, t2, t3, t4;
    __m256d tA, tB;

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
        t3  = _mm_load_pd(base + _mm_extract_epi32(offset, 2));
        t4  = _mm_load_pd(base + _mm_extract_epi32(offset, 3));
    }
    else
    {
        t1  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 0));
        t2  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 1));
        t3  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 2));
        t4  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 3));
    }

    tA  = _mm256_insertf128_pd(_mm256_castpd128_pd256(t1), t3, 0x1);
    tB  = _mm256_insertf128_pd(_mm256_castpd128_pd256(t2), t4, 0x1);
    v0  = _mm256_unpacklo_pd(tA, tB);
    v1  = _mm256_unpackhi_pd(tA, tB);
}


template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_loadu_bysimdint_transpose_d_avx_256(const double *       base,
                                                    gmx_simd_dint32_t    offset,
                                                    gmx_simd_double_t   &v0,
                                                    gmx_simd_double_t   &v1)
{
    __m128d t1, t2, t3, t4;
    __m256d tA, tB;

    /* Use optimized bit-shift multiply for the most common alignments. */
    /* Do nothing for align == 1 */
    if (align == 2)
    {
        offset = _mm_slli_epi32(offset, 1);
    }
    else if (align == 4)
    {
        offset = _mm_slli_epi32(offset, 2);
    }

    if (align == 1 || align == 2 || align == 4)
    {
        t1   = _mm_loadu_pd(base + _mm_extract_epi32(offset, 0));
        t2   = _mm_loadu_pd(base + _mm_extract_epi32(offset, 1));
        t3   = _mm_loadu_pd(base + _mm_extract_epi32(offset, 2));
        t4   = _mm_loadu_pd(base + _mm_extract_epi32(offset, 3));
    }
    else
    {
        t1   = _mm_loadu_pd(base + align * _mm_extract_epi32(offset, 0));
        t2   = _mm_loadu_pd(base + align * _mm_extract_epi32(offset, 1));
        t3   = _mm_loadu_pd(base + align * _mm_extract_epi32(offset, 2));
        t4   = _mm_loadu_pd(base + align * _mm_extract_epi32(offset, 3));
    }
    tA  = _mm256_insertf128_pd(_mm256_castpd128_pd256(t1), t3, 0x1);
    tB  = _mm256_insertf128_pd(_mm256_castpd128_pd256(t2), t4, 0x1);

    v0  = _mm256_unpacklo_pd(tA, tB);
    v1  = _mm256_unpackhi_pd(tA, tB);
}


static gmx_inline double gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_d_avx_256(double * m,
                                            gmx_simd_double_t v0, gmx_simd_double_t v1,
                                            gmx_simd_double_t v2, gmx_simd_double_t v3)
{
    __m256d t0, t1, t2;

    t0 = _mm256_hadd_pd(v0, v1);
    t1 = _mm256_hadd_pd(v2, v3);
    t2 = _mm256_permute2f128_pd(t0, t1, 0x21);
    t0 = _mm256_add_pd(t0, t2);
    t1 = _mm256_add_pd(t1, t2);
    t0 = _mm256_blend_pd(t0, t1, 0b1100);

    t1 = _mm256_add_pd(t0, _mm256_load_pd(m));
    _mm256_store_pd(m, t1);

    return gmx_simd_reduce_d_avx_256(t0);
}
#endif


#endif /* GMX_SIMD_IMPL_X86_AVX_256_UTIL_DOUBLE_H */
