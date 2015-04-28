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

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <immintrin.h>

#include "impl_x86_avx_256_simd_double.h"

namespace gmx
{

// Internal utility function: Full 4x4 transpose of __m256d
static inline void gmx_simdcall
avx256Transpose4By4(__m256d * v0,
                    __m256d * v1,
                    __m256d * v2,
                    __m256d * v3)
{
    __m256d t1 = _mm256_unpacklo_pd(*v0, *v1);
    __m256d t2 = _mm256_unpackhi_pd(*v0, *v1);
    __m256d t3 = _mm256_unpacklo_pd(*v2, *v3);
    __m256d t4 = _mm256_unpackhi_pd(*v2, *v3);
    *v0        = _mm256_permute2f128_pd(t1, t3, 0x20);
    *v1        = _mm256_permute2f128_pd(t2, t4, 0x20);
    *v2        = _mm256_permute2f128_pd(t1, t3, 0x31);
    *v3        = _mm256_permute2f128_pd(t2, t4, 0x31);
}

template <int align>
static inline void gmx_simdcall
simdGatherLoadTransposeD(const double *        base,
                         const std::int32_t    offset[],
                         SimdDouble *          v0,
                         SimdDouble *          v1,
                         SimdDouble *          v2,
                         SimdDouble *          v3)
{
    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 32 == 0);
    assert(align % 4 == 0);

    v0->r = _mm256_load_pd( base + align * offset[0] );
    v1->r = _mm256_load_pd( base + align * offset[1] );
    v2->r = _mm256_load_pd( base + align * offset[2] );
    v3->r = _mm256_load_pd( base + align * offset[3] );
    avx256Transpose4By4(&v0->r, &v1->r, &v2->r, &v3->r);
}

template <int align>
static inline void gmx_simdcall
simdGatherLoadTransposeD(const double *        base,
                         const std::int32_t    offset[],
                         SimdDouble *          v0,
                         SimdDouble *          v1)
{
    __m128d t1, t2, t3, t4;
    __m256d tA, tB;

    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 16 == 0);
    assert(align % 2 == 0);

    t1   = _mm_load_pd( base + align * offset[0] );
    t2   = _mm_load_pd( base + align * offset[1] );
    t3   = _mm_load_pd( base + align * offset[2] );
    t4   = _mm_load_pd( base + align * offset[3] );
    tA   = _mm256_insertf128_pd(_mm256_castpd128_pd256(t1), t3, 0x1);
    tB   = _mm256_insertf128_pd(_mm256_castpd128_pd256(t2), t4, 0x1);

    v0->r = _mm256_unpacklo_pd(tA, tB);
    v1->r = _mm256_unpackhi_pd(tA, tB);
}

static const int c_simdBestPairAlignmentD = 2;

template <int align>
static inline void gmx_simdcall
simdGatherLoadUTransposeD(const double *        base,
                          const std::int32_t    offset[],
                          SimdDouble *          v0,
                          SimdDouble *          v1,
                          SimdDouble *          v2)
{
    assert((size_t)offset % 16 == 0);

    __m256i mask = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_set_epi32(-1, -1, -1, -1)), _mm_set_epi32(0, 0, -1, -1), 0x1);
    __m256d t1, t2, t3, t4, t5, t6, t7, t8;
    t1    = _mm256_maskload_pd(base + align * offset[0], mask);
    t2    = _mm256_maskload_pd(base + align * offset[1], mask);
    t3    = _mm256_maskload_pd(base + align * offset[2], mask);
    t4    = _mm256_maskload_pd(base + align * offset[3], mask);
    t5    = _mm256_unpacklo_pd(t1, t2);
    t6    = _mm256_unpackhi_pd(t1, t2);
    t7    = _mm256_unpacklo_pd(t3, t4);
    t8    = _mm256_unpackhi_pd(t3, t4);
    v0->r = _mm256_permute2f128_pd(t5, t7, 0x20);
    v1->r = _mm256_permute2f128_pd(t6, t8, 0x20);
    v2->r = _mm256_permute2f128_pd(t5, t7, 0x31);
}

template <int align>
static inline void gmx_simdcall
simdTransposeScatterStoreUD(double *            base,
                            const std::int32_t  offset[],
                            SimdDouble          v0,
                            SimdDouble          v1,
                            SimdDouble          v2)
{
    __m256i mask = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_set_epi32(-1, -1, -1, -1)), _mm_set_epi32(0, 0, -1, -1), 0x1);
    __m256d tv3  = _mm256_setzero_pd();

    assert((size_t)offset % 16 == 0);

    avx256Transpose4By4(&v0.r, &v1.r, &v2.r, &tv3);

    _mm256_maskstore_pd(base + align * offset[0], mask, v0.r);
    _mm256_maskstore_pd(base + align * offset[1], mask, v1.r);
    _mm256_maskstore_pd(base + align * offset[2], mask, v2.r);
    _mm256_maskstore_pd(base + align * offset[3], mask, tv3);
}

template <int align>
static inline void gmx_simdcall
simdTransposeScatterIncrUD(double *            base,
                           const std::int32_t  offset[],
                           SimdDouble          v0,
                           SimdDouble          v1,
                           SimdDouble          v2)
{
    __m256i mask = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_set_epi32(-1, -1, -1, -1)), _mm_set_epi32(0, 0, -1, -1), 0x1);
    __m256d tv3  = _mm256_setzero_pd();
    __m256d t0, t1, t2, t3;

    assert((size_t)offset % 16 == 0);

    t0  = _mm256_maskload_pd(base + align * offset[0], mask);
    t1  = _mm256_maskload_pd(base + align * offset[1], mask);
    t2  = _mm256_maskload_pd(base + align * offset[2], mask);
    t3  = _mm256_maskload_pd(base + align * offset[3], mask);

    avx256Transpose4By4(&v0.r, &v1.r, &v2.r, &tv3);

    t0  = _mm256_add_pd(t0, v0.r);
    t1  = _mm256_add_pd(t1, v1.r);
    t2  = _mm256_add_pd(t2, v2.r);
    t3  = _mm256_add_pd(t3, tv3);

    _mm256_maskstore_pd(base + align * offset[0], mask, t0);
    _mm256_maskstore_pd(base + align * offset[1], mask, t1);
    _mm256_maskstore_pd(base + align * offset[2], mask, t2);
    _mm256_maskstore_pd(base + align * offset[3], mask, t3);
}

template <int align>
static inline void gmx_simdcall
simdTransposeScatterDecrUD(double *            base,
                           const std::int32_t  offset[],
                           SimdDouble          v0,
                           SimdDouble          v1,
                           SimdDouble          v2)
{
    __m256i mask = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_set_epi32(-1, -1, -1, -1)), _mm_set_epi32(0, 0, -1, -1), 0x1);
    __m256d tv3  = _mm256_setzero_pd();
    __m256d t0, t1, t2, t3;

    assert((size_t)offset % 16 == 0);

    t0  = _mm256_maskload_pd(base + align * offset[0], mask);
    t1  = _mm256_maskload_pd(base + align * offset[1], mask);
    t2  = _mm256_maskload_pd(base + align * offset[2], mask);
    t3  = _mm256_maskload_pd(base + align * offset[3], mask);

    avx256Transpose4By4(&v0.r, &v1.r, &v2.r, &tv3);

    t0  = _mm256_sub_pd(t0, v0.r);
    t1  = _mm256_sub_pd(t1, v1.r);
    t2  = _mm256_sub_pd(t2, v2.r);
    t3  = _mm256_sub_pd(t3, tv3);

    _mm256_maskstore_pd(base + align * offset[0], mask, t0);
    _mm256_maskstore_pd(base + align * offset[1], mask, t1);
    _mm256_maskstore_pd(base + align * offset[2], mask, t2);
    _mm256_maskstore_pd(base + align * offset[3], mask, t3);
}

static inline void gmx_simdcall
simdExpandScalarsToTripletsD(SimdDouble    scalar,
                             SimdDouble *  triplets0,
                             SimdDouble *  triplets1,
                             SimdDouble *  triplets2)
{
    __m256d t0 = _mm256_permute2f128_pd(scalar.r, scalar.r, 0x21);
    __m256d t1 = _mm256_permute_pd(scalar.r, 0b0000);
    __m256d t2 = _mm256_permute_pd(scalar.r, 0b1111);
    triplets0->r = _mm256_blend_pd(t1, t0, 0b1100);
    triplets1->r = _mm256_blend_pd(t2, t1, 0b1100);
    triplets2->r = _mm256_blend_pd(t0, t2, 0b1100);
}

template <int align>
static inline void gmx_simdcall
simdGatherLoadBySimdIntTransposeD(const double *  base,
                                  SimdDInt32      offset,
                                  SimdDouble *    v0,
                                  SimdDouble *    v1,
                                  SimdDouble *    v2,
                                  SimdDouble *    v3)
{
    assert((size_t)base % 32 == 0);
    assert(align % 4 == 0);

    // Use optimized bit-shift multiply for the most common alignments
    if (align == 4)
    {
        offset.i = _mm_slli_epi32(offset.i, 2);
    }
    else if (align == 8)
    {
        offset.i = _mm_slli_epi32(offset.i, 3);
    }
    else if (align == 12)
    {
        // multiply by 3, then by 4
        offset.i = _mm_add_epi32(offset.i, _mm_slli_epi32(offset.i, 1));
        offset.i = _mm_slli_epi32(offset.i, 2);
    }
    else if (align == 16)
    {
        offset.i = _mm_slli_epi32(offset.i, 4);
    }

    if (align == 4 || align == 8 || align == 12 || align == 16)
    {
        v0->r = _mm256_load_pd(base + _mm_extract_epi32(offset.i, 0));
        v1->r = _mm256_load_pd(base + _mm_extract_epi32(offset.i, 1));
        v2->r = _mm256_load_pd(base + _mm_extract_epi32(offset.i, 2));
        v3->r = _mm256_load_pd(base + _mm_extract_epi32(offset.i, 3));
    }
    else
    {
        v0->r = _mm256_load_pd(base + align * _mm_extract_epi32(offset.i, 0));
        v1->r = _mm256_load_pd(base + align * _mm_extract_epi32(offset.i, 1));
        v2->r = _mm256_load_pd(base + align * _mm_extract_epi32(offset.i, 2));
        v3->r = _mm256_load_pd(base + align * _mm_extract_epi32(offset.i, 3));
    }
    avx256Transpose4By4(&v0->r, &v1->r, &v2->r, &v3->r);
}

template <int align>
static inline void gmx_simdcall
simdGatherLoadBySimdIntTransposeD(const double *    base,
                                  SimdDInt32        offset,
                                  SimdDouble *      v0,
                                  SimdDouble *      v1)
{
    __m128d t1, t2, t3, t4;
    __m256d tA, tB;

    assert((size_t)base % 16 == 0);
    assert(align % 2 == 0);

    // Use optimized bit-shift multiply for the most common alignments
    if (align == 2)
    {
        offset.i = _mm_slli_epi32(offset.i, 1);
    }
    else if (align == 4)
    {
        offset.i = _mm_slli_epi32(offset.i, 2);
    }
    else if (align == 6)
    {
        // multiply by 3, then by 2
        offset.i = _mm_add_epi32(offset.i, _mm_slli_epi32(offset.i, 1));
        offset.i = _mm_slli_epi32(offset.i, 1);
    }
    else if (align == 8)
    {
        offset.i = _mm_slli_epi32(offset.i, 3);
    }
    else if (align == 12)
    {
        // multiply by 3, then by 4
        offset.i = _mm_add_epi32(offset.i, _mm_slli_epi32(offset.i, 1));
        offset.i = _mm_slli_epi32(offset.i, 2);
    }
    else if (align == 16)
    {
        offset.i = _mm_slli_epi32(offset.i, 4);
    }

    if (align == 2 || align == 4 || align == 6 ||
        align == 8 || align == 12 || align == 16)
    {
        t1  = _mm_load_pd(base + _mm_extract_epi32(offset.i, 0));
        t2  = _mm_load_pd(base + _mm_extract_epi32(offset.i, 1));
        t3  = _mm_load_pd(base + _mm_extract_epi32(offset.i, 2));
        t4  = _mm_load_pd(base + _mm_extract_epi32(offset.i, 3));
    }
    else
    {
        t1  = _mm_load_pd(base + align * _mm_extract_epi32(offset.i, 0));
        t2  = _mm_load_pd(base + align * _mm_extract_epi32(offset.i, 1));
        t3  = _mm_load_pd(base + align * _mm_extract_epi32(offset.i, 2));
        t4  = _mm_load_pd(base + align * _mm_extract_epi32(offset.i, 3));
    }

    tA    = _mm256_insertf128_pd(_mm256_castpd128_pd256(t1), t3, 0x1);
    tB    = _mm256_insertf128_pd(_mm256_castpd128_pd256(t2), t4, 0x1);
    v0->r = _mm256_unpacklo_pd(tA, tB);
    v1->r = _mm256_unpackhi_pd(tA, tB);
}

template <int align>
static inline void gmx_simdcall
simdGatherLoadUBySimdIntTransposeD(const double *  base,
                                   SimdDInt32      offset,
                                   SimdDouble *    v0,
                                   SimdDouble *    v1)
{
    __m128d t1, t2, t3, t4;
    __m256d tA, tB;

    // Use optimized bit-shift multiply for the most common alignments.

    // Do nothing for align == 1
    if (align == 2)
    {
        offset.i = _mm_slli_epi32(offset.i, 1);
    }
    else if (align == 4)
    {
        offset.i = _mm_slli_epi32(offset.i, 2);
    }

    if (align == 1 || align == 2 || align == 4)
    {
        t1   = _mm_loadu_pd(base + _mm_extract_epi32(offset.i, 0));
        t2   = _mm_loadu_pd(base + _mm_extract_epi32(offset.i, 1));
        t3   = _mm_loadu_pd(base + _mm_extract_epi32(offset.i, 2));
        t4   = _mm_loadu_pd(base + _mm_extract_epi32(offset.i, 3));
    }
    else
    {
        t1   = _mm_loadu_pd(base + align * _mm_extract_epi32(offset.i, 0));
        t2   = _mm_loadu_pd(base + align * _mm_extract_epi32(offset.i, 1));
        t3   = _mm_loadu_pd(base + align * _mm_extract_epi32(offset.i, 2));
        t4   = _mm_loadu_pd(base + align * _mm_extract_epi32(offset.i, 3));
    }
    tA  = _mm256_insertf128_pd(_mm256_castpd128_pd256(t1), t3, 0x1);
    tB  = _mm256_insertf128_pd(_mm256_castpd128_pd256(t2), t4, 0x1);

    v0->r = _mm256_unpacklo_pd(tA, tB);
    v1->r = _mm256_unpackhi_pd(tA, tB);
}

static inline double gmx_simdcall
simdReduceIncr4ReturnSumD(double *    m,
                          SimdDouble  v0,
                          SimdDouble  v1,
                          SimdDouble  v2,
                          SimdDouble  v3)
{
    __m256d t0, t1, t2;
    __m128d a0, a1;

    assert((size_t)m % 32 == 0);

    t0 = _mm256_hadd_pd(v0.r, v1.r);
    t1 = _mm256_hadd_pd(v2.r, v3.r);
    t2 = _mm256_permute2f128_pd(t0, t1, 0x21);
    t0 = _mm256_add_pd(t0, t2);
    t1 = _mm256_add_pd(t1, t2);
    t0 = _mm256_blend_pd(t0, t1, 0b1100);

    t1 = _mm256_add_pd(t0, _mm256_load_pd(m));
    _mm256_store_pd(m, t1);

    t0  = _mm256_add_pd(t0, _mm256_permute_pd(t0, 0b0101 ));
    a0  = _mm256_castpd256_pd128(t0);
    a1  = _mm256_extractf128_pd(t0, 0x1);
    a0  = _mm_add_sd(a0, a1);

    return *reinterpret_cast<double *>(&a0);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_256_UTIL_DOUBLE_H
