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

#ifndef GMX_SIMD_IMPL_X86_AVX2_256_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX2_256_UTIL_FLOAT_H

#include "config.h"

#include <immintrin.h>

#include "impl_x86_avx2_256_common.h"

/* Higher-level utility functions */
#undef  gmx_simd_gather_load_bysimdint_transpose_f
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_avx2_256
#undef  gmx_simd_gather_loadu_bysimdint_transpose_f
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_avx2_256

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus
template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_f_avx2_256(const float *       base,
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

    /* AVX2 has 8-way integer simd (bitwise shift) support.
     * Use optimized bit-shift multiply for the most common alignments.
     */
    if (align == 4)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        simdoffset = _mm256_add_epi32(simdoffset, _mm256_slli_epi32(simdoffset, 1));
        simdoffset = _mm256_slli_epi32(simdoffset, 2);
    }
    else if (align == 16)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 4);
    }

    _mm256_store_si256( (__m256i *)offset, simdoffset);

    if (align == 4 || align == 8 || align == 12 || align == 16)
    {
        t1  = _mm_load_ps( base + offset[0] );
        t2  = _mm_load_ps( base + offset[1] );
        t3  = _mm_load_ps( base + offset[2] );
        t4  = _mm_load_ps( base + offset[3] );
        t5  = _mm_load_ps( base + offset[4] );
        t6  = _mm_load_ps( base + offset[5] );
        t7  = _mm_load_ps( base + offset[6] );
        t8  = _mm_load_ps( base + offset[7] );
    }
    else
    {
        t1  = _mm_load_ps( base + align * offset[0] );
        t2  = _mm_load_ps( base + align * offset[1] );
        t3  = _mm_load_ps( base + align * offset[2] );
        t4  = _mm_load_ps( base + align * offset[3] );
        t5  = _mm_load_ps( base + align * offset[4] );
        t6  = _mm_load_ps( base + align * offset[5] );
        t7  = _mm_load_ps( base + align * offset[6] );
        t8  = _mm_load_ps( base + align * offset[7] );
    }

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
gmx_simd_gather_load_bysimdint_transpose_f_avx2_256(const float *       base,
                                                    gmx_simd_fint32_t   simdoffset,
                                                    gmx_simd_float_t   &v0,
                                                    gmx_simd_float_t   &v1)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;
    /* All x86 compilers we know that support SIMD also support GCC alignment attributes */
    int __attribute__ ((aligned (32))) offset[GMX_SIMD_FLOAT_WIDTH];

    /* AVX2 has 8-way integer simd (bitwise shift) support.
     * Use optimized bit-shift multiply for the most common alignments.
     */
    if (align == 2)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 1);
    }
    else if (align == 4)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 2);
    }
    else if (align == 6)
    {
        /* multiply by 3, then by 2 */
        simdoffset = _mm256_add_epi32(simdoffset, _mm256_slli_epi32(simdoffset, 1));
        simdoffset = _mm256_slli_epi32(simdoffset, 1);
    }
    else if (align == 8)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        simdoffset = _mm256_add_epi32(simdoffset, _mm256_slli_epi32(simdoffset, 1));
        simdoffset = _mm256_slli_epi32(simdoffset, 2);
    }
    else if (align == 16)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 4);
    }

    _mm256_store_si256( (__m256i *)offset, simdoffset);

    if (align == 2 || align == 4 || align == 6 ||
        align == 8 || align == 12 || align == 16)
    {
        t1  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[0] ) );
        t2  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[1] ) );
        t3  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[2] ) );
        t4  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[3] ) );
        t5  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[4] ) );
        t6  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[5] ) );
        t7  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[6] ) );
        t8  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[7] ) );
    }
    else
    {
        t1  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[0] ) );
        t2  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[1] ) );
        t3  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[2] ) );
        t4  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[3] ) );
        t5  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[4] ) );
        t6  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[5] ) );
        t7  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[6] ) );
        t8  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[7] ) );
    }

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
gmx_simd_gather_loadu_bysimdint_transpose_f_avx2_256(const float *       base,
                                                     gmx_simd_fint32_t   simdoffset,
                                                     gmx_simd_float_t   &v0,
                                                     gmx_simd_float_t   &v1)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;
    /* All x86 compilers we know that support SIMD also support GCC alignment attributes */
    int __attribute__ ((aligned (32))) offset[GMX_SIMD_FLOAT_WIDTH];

    /* AVX2 has 8-way integer simd (bitwise shift) support.
     * Use optimized bit-shift multiply for the most common alignments.
     */
    /* Do nothing for align == 1 */
    if (align == 2)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 1);
    }
    else if (align == 4)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 2);
    }

    _mm256_store_si256( (__m256i *)offset, simdoffset);

    if (align == 1 || align == 2 || align == 4)
    {
        t1  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[0] ) );
        t2  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[1] ) );
        t3  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[2] ) );
        t4  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[3] ) );
        t5  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[4] ) );
        t6  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[5] ) );
        t7  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[6] ) );
        t8  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + offset[7] ) );
    }
    else
    {
        t1  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[0] ) );
        t2  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[1] ) );
        t3  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[2] ) );
        t4  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[3] ) );
        t5  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[4] ) );
        t6  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[5] ) );
        t7  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[6] ) );
        t8  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[7] ) );
    }

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(tA, tC);
    tB  = _mm256_unpacklo_ps(tB, tD);
    v0  = _mm256_unpacklo_ps(tA, tB);
    v1  = _mm256_unpackhi_ps(tA, tB);
}
#endif

#endif /* GMX_SIMD_IMPL_X86_AVX2_256_UTIL_FLOAT_H */
