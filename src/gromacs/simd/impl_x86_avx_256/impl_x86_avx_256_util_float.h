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

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <immintrin.h>

#include "gromacs/utility/basedefinitions.h"

#include "impl_x86_avx_256_simd_float.h"

namespace gmx
{

/* This is an internal helper function used by the three functions storing,
 * incrementing, or decrementing data. Do NOT use it outside this file.
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
static inline void gmx_simdcall
avx256Transpose3By4InLanes(__m256 * v0,
                           __m256 * v1,
                           __m256 * v2,
                           __m256 * v3)
{
    __m256 t1 = _mm256_unpacklo_ps(*v0, *v1);
    __m256 t2 = _mm256_unpackhi_ps(*v0, *v1);
    *v0       = _mm256_shuffle_ps(t1, *v2, _MM_SHUFFLE(0, 0, 1, 0));
    *v1       = _mm256_shuffle_ps(t1, *v2, _MM_SHUFFLE(0, 1, 3, 2));
    *v3       = _mm256_shuffle_ps(t2, *v2, _MM_SHUFFLE(0, 3, 3, 2));
    *v2       = _mm256_shuffle_ps(t2, *v2, _MM_SHUFFLE(0, 2, 1, 0));
}

template <int align>
static inline void gmx_simdcall
simdGatherLoadTransposeF(const float *        base,
                         const std::int32_t   offset[],
                         SimdFloat *          v0,
                         SimdFloat *          v1,
                         SimdFloat *          v2,
                         SimdFloat *          v3)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;

    assert(std::size_t(offset) % 32 == 0);
    assert(std::size_t(base) % 16 == 0);
    assert(align % 4 == 0);

    t1  = _mm_load_ps( base + align * offset[0] );
    t2  = _mm_load_ps( base + align * offset[1] );
    t3  = _mm_load_ps( base + align * offset[2] );
    t4  = _mm_load_ps( base + align * offset[3] );
    t5  = _mm_load_ps( base + align * offset[4] );
    t6  = _mm_load_ps( base + align * offset[5] );
    t7  = _mm_load_ps( base + align * offset[6] );
    t8  = _mm_load_ps( base + align * offset[7] );

    v0->r = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    v1->r = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    v2->r = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    v3->r = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(v0->r, v1->r);
    tB  = _mm256_unpacklo_ps(v2->r, v3->r);
    tC  = _mm256_unpackhi_ps(v0->r, v1->r);
    tD  = _mm256_unpackhi_ps(v2->r, v3->r);

    v0->r = _mm256_shuffle_ps(tA, tB, _MM_SHUFFLE(1, 0, 1, 0));
    v1->r = _mm256_shuffle_ps(tA, tB, _MM_SHUFFLE(3, 2, 3, 2));
    v2->r = _mm256_shuffle_ps(tC, tD, _MM_SHUFFLE(1, 0, 1, 0));
    v3->r = _mm256_shuffle_ps(tC, tD, _MM_SHUFFLE(3, 2, 3, 2));
}

template <int align>
static inline void gmx_simdcall
simdGatherLoadTransposeF(const float *        base,
                         const std::int32_t   offset[],
                         SimdFloat *          v0,
                         SimdFloat *          v1)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;

    assert(std::size_t(offset) % 32 == 0);
    assert(std::size_t(base) % 8 == 0);
    assert(align % 2 == 0);

    t1  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[0] ) );
    t2  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[1] ) );
    t3  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[2] ) );
    t4  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[3] ) );
    t5  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[4] ) );
    t6  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[5] ) );
    t7  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[6] ) );
    t8  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[7] ) );

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA    = _mm256_unpacklo_ps(tA, tC);
    tB    = _mm256_unpacklo_ps(tB, tD);
    v0->r = _mm256_unpacklo_ps(tA, tB);
    v1->r = _mm256_unpackhi_ps(tA, tB);
}

static const int c_simdBestPairAlignmentF = 2;

template <int align>
static inline void gmx_simdcall
simdGatherLoadUTransposeF(const float *        base,
                          const std::int32_t   offset[],
                          SimdFloat *          v0,
                          SimdFloat *          v1,
                          SimdFloat *          v2)
{
    __m256  t1, t2, t3, t4, t5, t6, t7, t8;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    assert(std::size_t(offset) % 32 == 0);

    if ( (align & 0x3) == 0)
    {
        // With alignment 4 or better we can read a byte beyond triplets and use _mm_loadu_ps()
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
        // Arbitrary alignment
        t1  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_maskload_ps( base + align * offset[0], mask)),
                                   _mm_maskload_ps( base + align * offset[4], mask), 0x1);
        t2  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_maskload_ps( base + align * offset[1], mask)),
                                   _mm_maskload_ps( base + align * offset[5], mask), 0x1);
        t3  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_maskload_ps( base + align * offset[2], mask)),
                                   _mm_maskload_ps( base + align * offset[6], mask), 0x1);
        t4  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_maskload_ps( base + align * offset[3], mask)),
                                   _mm_maskload_ps( base + align * offset[7], mask), 0x1);
    }
    t5    = _mm256_unpacklo_ps(t1, t2);
    t6    = _mm256_unpacklo_ps(t3, t4);
    t7    = _mm256_unpackhi_ps(t1, t2);
    t8    = _mm256_unpackhi_ps(t3, t4);
    v0->r = _mm256_shuffle_ps(t5, t6, _MM_SHUFFLE(1, 0, 1, 0));
    v1->r = _mm256_shuffle_ps(t5, t6, _MM_SHUFFLE(3, 2, 3, 2));
    v2->r = _mm256_shuffle_ps(t7, t8, _MM_SHUFFLE(1, 0, 1, 0));
}

template <int align>
static inline void gmx_simdcall
simdTransposeScatterStoreUF(float *              base,
                            const std::int32_t   offset[],
                            SimdFloat            v0,
                            SimdFloat            v1,
                            SimdFloat            v2)
{
    __m256  tv3;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    assert(std::size_t(offset) % 32 == 0);

    avx256Transpose3By4InLanes(&v0.r, &v1.r, &v2.r, &tv3);
    _mm_maskstore_ps( base + align * offset[0], mask, _mm256_castps256_ps128(v0.r));
    _mm_maskstore_ps( base + align * offset[1], mask, _mm256_castps256_ps128(v1.r));
    _mm_maskstore_ps( base + align * offset[2], mask, _mm256_castps256_ps128(v2.r));
    _mm_maskstore_ps( base + align * offset[3], mask, _mm256_castps256_ps128(tv3));
    _mm_maskstore_ps( base + align * offset[4], mask, _mm256_extractf128_ps(v0.r, 0x1));
    _mm_maskstore_ps( base + align * offset[5], mask, _mm256_extractf128_ps(v1.r, 0x1));
    _mm_maskstore_ps( base + align * offset[6], mask, _mm256_extractf128_ps(v2.r, 0x1));
    _mm_maskstore_ps( base + align * offset[7], mask, _mm256_extractf128_ps(tv3, 0x1));
}

template <int align>
static inline void gmx_simdcall
simdTransposeScatterIncrUF(float *              base,
                           const std::int32_t   offset[],
                           SimdFloat            v0,
                           SimdFloat            v1,
                           SimdFloat            v2)
{
    __m256  tv3, t0, t1, t2, t3;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    assert(std::size_t(offset) % 32 == 0);

    t0  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_maskload_ps( base + align * offset[0], mask)),
                               _mm_maskload_ps( base + align * offset[4], mask), 0x1);
    t1  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_maskload_ps( base + align * offset[1], mask)),
                               _mm_maskload_ps( base + align * offset[5], mask), 0x1);
    t2  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_maskload_ps( base + align * offset[2], mask)),
                               _mm_maskload_ps( base + align * offset[6], mask), 0x1);
    t3  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_maskload_ps( base + align * offset[3], mask)),
                               _mm_maskload_ps( base + align * offset[7], mask), 0x1);

    avx256Transpose3By4InLanes(&v0.r, &v1.r, &v2.r, &tv3);
    t0          = _mm256_add_ps(t0, v0.r);
    t1          = _mm256_add_ps(t1, v1.r);
    t2          = _mm256_add_ps(t2, v2.r);
    t3          = _mm256_add_ps(t3, tv3);

    _mm_maskstore_ps(base + align * offset[0], mask, _mm256_castps256_ps128(t0));
    _mm_maskstore_ps(base + align * offset[1], mask, _mm256_castps256_ps128(t1));
    _mm_maskstore_ps(base + align * offset[2], mask, _mm256_castps256_ps128(t2));
    _mm_maskstore_ps(base + align * offset[3], mask, _mm256_castps256_ps128(t3));
    _mm_maskstore_ps(base + align * offset[4], mask, _mm256_extractf128_ps(t0, 0x1));
    _mm_maskstore_ps(base + align * offset[5], mask, _mm256_extractf128_ps(t1, 0x1));
    _mm_maskstore_ps(base + align * offset[6], mask, _mm256_extractf128_ps(t2, 0x1));
    _mm_maskstore_ps(base + align * offset[7], mask, _mm256_extractf128_ps(t3, 0x1));
}

template <int align>
static inline void gmx_simdcall
simdTransposeScatterDecrUF(float *              base,
                           const std::int32_t   offset[],
                           SimdFloat            v0,
                           SimdFloat            v1,
                           SimdFloat            v2)
{
    __m256  tv3, t0, t1, t2, t3;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    assert(std::size_t(offset) % 32 == 0);

    t0  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_maskload_ps( base + align * offset[0], mask)),
                               _mm_maskload_ps( base + align * offset[4], mask), 0x1);
    t1  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_maskload_ps( base + align * offset[1], mask)),
                               _mm_maskload_ps( base + align * offset[5], mask), 0x1);
    t2  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_maskload_ps( base + align * offset[2], mask)),
                               _mm_maskload_ps( base + align * offset[6], mask), 0x1);
    t3  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_maskload_ps( base + align * offset[3], mask)),
                               _mm_maskload_ps( base + align * offset[7], mask), 0x1);

    avx256Transpose3By4InLanes(&v0.r, &v1.r, &v2.r, &tv3);
    t0          = _mm256_sub_ps(t0, v0.r);
    t1          = _mm256_sub_ps(t1, v1.r);
    t2          = _mm256_sub_ps(t2, v2.r);
    t3          = _mm256_sub_ps(t3, tv3);

    _mm_maskstore_ps(base + align * offset[0], mask, _mm256_castps256_ps128(t0));
    _mm_maskstore_ps(base + align * offset[1], mask, _mm256_castps256_ps128(t1));
    _mm_maskstore_ps(base + align * offset[2], mask, _mm256_castps256_ps128(t2));
    _mm_maskstore_ps(base + align * offset[3], mask, _mm256_castps256_ps128(t3));
    _mm_maskstore_ps(base + align * offset[4], mask, _mm256_extractf128_ps(t0, 0x1));
    _mm_maskstore_ps(base + align * offset[5], mask, _mm256_extractf128_ps(t1, 0x1));
    _mm_maskstore_ps(base + align * offset[6], mask, _mm256_extractf128_ps(t2, 0x1));
    _mm_maskstore_ps(base + align * offset[7], mask, _mm256_extractf128_ps(t3, 0x1));
}

static inline void gmx_simdcall
simdExpandScalarsToTripletsF(SimdFloat    scalar,
                             SimdFloat *  triplets0,
                             SimdFloat *  triplets1,
                             SimdFloat *  triplets2)
{
    __m256 t0 = _mm256_permute2f128_ps(scalar.r, scalar.r, 0x21);
    __m256 t1 = _mm256_permute_ps(scalar.r, _MM_SHUFFLE(1, 0, 0, 0));
    __m256 t2 = _mm256_permute_ps(t0, _MM_SHUFFLE(2, 2, 1, 1));
    __m256 t3 = _mm256_permute_ps(scalar.r, _MM_SHUFFLE(3, 3, 3, 2));
    triplets0->r = _mm256_blend_ps(t1, t2, 0xF0);
    triplets1->r = _mm256_blend_ps(t3, t1, 0xF0);
    triplets2->r = _mm256_blend_ps(t2, t3, 0xF0);
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
template <int align>
static inline void gmx_simdcall
simdGatherLoadBySimdIntTransposeF(const float *  base,
                                  SimdFInt32     simdoffset,
                                  SimdFloat *    v0,
                                  SimdFloat *    v1,
                                  SimdFloat *    v2,
                                  SimdFloat *    v3)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;
    GMX_ALIGNED(int, GMX_SIMD_FLOAT_WIDTH) offset[GMX_SIMD_FLOAT_WIDTH];

    assert(std::size_t(base) % 16 == 0);
    assert(align % 4 == 0);

    // Plain AVX does not support 256-bit integer shift operations,
    // so we always store to memory and multiply in the integer units instead.
    _mm256_store_si256( reinterpret_cast<__m256i *>(offset), simdoffset.i);

    t1  = _mm_load_ps( base + align * offset[0] );
    t2  = _mm_load_ps( base + align * offset[1] );
    t3  = _mm_load_ps( base + align * offset[2] );
    t4  = _mm_load_ps( base + align * offset[3] );
    t5  = _mm_load_ps( base + align * offset[4] );
    t6  = _mm_load_ps( base + align * offset[5] );
    t7  = _mm_load_ps( base + align * offset[6] );
    t8  = _mm_load_ps( base + align * offset[7] );

    v0->r = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    v1->r = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    v2->r = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    v3->r = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(v0->r, v1->r);
    tB  = _mm256_unpacklo_ps(v2->r, v3->r);
    tC  = _mm256_unpackhi_ps(v0->r, v1->r);
    tD  = _mm256_unpackhi_ps(v2->r, v3->r);

    v0->r = _mm256_shuffle_ps(tA, tB, _MM_SHUFFLE(1, 0, 1, 0));
    v1->r = _mm256_shuffle_ps(tA, tB, _MM_SHUFFLE(3, 2, 3, 2));
    v2->r = _mm256_shuffle_ps(tC, tD, _MM_SHUFFLE(1, 0, 1, 0));
    v3->r = _mm256_shuffle_ps(tC, tD, _MM_SHUFFLE(3, 2, 3, 2));
}

template <int align>
static inline void gmx_simdcall
simdGatherLoadBySimdIntTransposeF(const float *   base,
                                  SimdFInt32      simdoffset,
                                  SimdFloat *     v0,
                                  SimdFloat *     v1)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;
    GMX_ALIGNED(int, GMX_SIMD_FLOAT_WIDTH) offset[GMX_SIMD_FLOAT_WIDTH];

    assert(std::size_t(base) % 8 == 0);
    assert(align % 2 == 0);

    // Plain AVX does not support 256-bit integer shift operations,
    // so we always store to memory and multiply in the integer units instead.
    _mm256_store_si256( reinterpret_cast<__m256i *>(offset), simdoffset.i);

    t1  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[0] ) );
    t2  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[1] ) );
    t3  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[2] ) );
    t4  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[3] ) );
    t5  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[4] ) );
    t6  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[5] ) );
    t7  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[6] ) );
    t8  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[7] ) );

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA    = _mm256_unpacklo_ps(tA, tC);
    tB    = _mm256_unpacklo_ps(tB, tD);
    v0->r = _mm256_unpacklo_ps(tA, tB);
    v1->r = _mm256_unpackhi_ps(tA, tB);
}


template <int align>
static inline void gmx_simdcall
simdGatherLoadUBySimdIntTransposeF(const float *  base,
                                   SimdFInt32     simdoffset,
                                   SimdFloat *    v0,
                                   SimdFloat *    v1)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;
    GMX_ALIGNED(int, GMX_SIMD_FLOAT_WIDTH) offset[GMX_SIMD_FLOAT_WIDTH];

    // Plain AVX does not support 256-bit integer shift operations,
    // so we always store to memory and multiply in the integer units instead.
    _mm256_store_si256( reinterpret_cast<__m256i *>(offset), simdoffset.i);

    t1  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[0] ) );
    t2  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[1] ) );
    t3  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[2] ) );
    t4  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[3] ) );
    t5  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[4] ) );
    t6  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[5] ) );
    t7  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[6] ) );
    t8  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>( base + align * offset[7] ) );

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA    = _mm256_unpacklo_ps(tA, tC);
    tB    = _mm256_unpacklo_ps(tB, tD);
    v0->r = _mm256_unpacklo_ps(tA, tB);
    v1->r = _mm256_unpackhi_ps(tA, tB);
}
#endif // GMX_SIMD_X86_AVX_256

static inline float gmx_simdcall
simdReduceIncr4ReturnSumF(float *    m,
                          SimdFloat  v0,
                          SimdFloat  v1,
                          SimdFloat  v2,
                          SimdFloat  v3)
{
    __m128 t0, t2;

    assert(std::size_t(m) % 16 == 0);

    v0.r = _mm256_hadd_ps(v0.r, v1.r);
    v2.r = _mm256_hadd_ps(v2.r, v3.r);
    v0.r = _mm256_hadd_ps(v0.r, v2.r);
    t0   = _mm_add_ps(_mm256_castps256_ps128(v0.r), _mm256_extractf128_ps(v0.r, 0x1));

    t2 = _mm_add_ps(t0, _mm_load_ps(m));
    _mm_store_ps(m, t2);

    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float *>(&t0);
}


/*************************************
 * Half-simd-width utility functions *
 *************************************/
static inline SimdFloat gmx_simdcall
simdLoadDualHsimdF(const float * m0,
                   const float * m1)
{
    assert((size_t)m0 % 16 == 0);
    assert((size_t)m1 % 16 == 0);

    return {
               _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_load_ps(m0)), _mm_load_ps(m1), 0x1)
    };
}

static inline SimdFloat gmx_simdcall
simdLoadDupHsimdF(const float * m)
{
    assert((size_t)m % 16 == 0);

    return {
               _mm256_broadcast_ps(reinterpret_cast<const __m128 *>(m))
    };
}

static inline SimdFloat gmx_simdcall
simdLoad1DualHsimdF(const float * m)
{
    __m128 t0, t1;
    t0 = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(m));
    t1 = _mm_permute_ps(t0, _MM_SHUFFLE(1, 1, 1, 1));
    t0 = _mm_permute_ps(t0, _MM_SHUFFLE(0, 0, 0, 0));

    return {
               _mm256_insertf128_ps(_mm256_castps128_ps256(t0), t1, 0x1)
    };
}


static inline void gmx_simdcall
simdStoreDualHsimdF(float *     m0,
                    float *     m1,
                    SimdFloat   a)
{
    assert((size_t)m0 % 16 == 0);
    assert((size_t)m1 % 16 == 0);
    _mm_store_ps(m0, _mm256_castps256_ps128(a.r));
    _mm_store_ps(m1, _mm256_extractf128_ps(a.r, 0x1));
}

static inline void gmx_simdcall
simdDecrHsimdF(float *    m,
               SimdFloat  a)
{
    assert((size_t)m % 16 == 0);
    __m128 asum = _mm_add_ps(_mm256_castps256_ps128(a.r), _mm256_extractf128_ps(a.r, 0x1));
    _mm_store_ps(m, _mm_sub_ps(_mm_load_ps(m), asum));
}


template <int align>
static inline void gmx_simdcall
simdGatherLoadTransposeHsimdF(const float *    base0,
                              const float *    base1,
                              std::int32_t     offset[],
                              SimdFloat *      v0,
                              SimdFloat *      v1)
{
    __m128 t0, t1, t2, t3, t4, t5, t6, t7;
    __m256 tA, tB, tC, tD;

    assert((size_t)offset % 16 == 0);
    assert((size_t)base0 % 8 == 0);
    assert((size_t)base1 % 8 == 0);
    assert(align % 2 == 0);

    t0  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base0 + align * offset[0]));
    t1  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base0 + align * offset[1]));
    t2  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base0 + align * offset[2]));
    t3  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base0 + align * offset[3]));
    t4  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base1 + align * offset[0]));
    t5  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base1 + align * offset[1]));
    t6  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base1 + align * offset[2]));
    t7  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base1 + align * offset[3]));

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t0), t4, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);

    tA    = _mm256_unpacklo_ps(tA, tC);
    tB    = _mm256_unpacklo_ps(tB, tD);
    v0->r = _mm256_unpacklo_ps(tA, tB);
    v1->r = _mm256_unpackhi_ps(tA, tB);
}


static inline float gmx_simdcall
simdReduceIncr4ReturnSumHsimdF(float *     m,
                               SimdFloat   v0,
                               SimdFloat   v1)
{
    __m128 t0, t1;

    v0.r = _mm256_hadd_ps(v0.r, v1.r);
    t0   = _mm256_extractf128_ps(v0.r, 0x1);
    t0   = _mm_hadd_ps(_mm256_castps256_ps128(v0.r), t0);
    t0   = _mm_permute_ps(t0, _MM_SHUFFLE(3, 1, 2, 0));

    assert((size_t)m % 16 == 0);

    t1   = _mm_add_ps(t0, _mm_load_ps(m));
    _mm_store_ps(m, t1);

    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float *>(&t0);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_256_UTIL_FLOAT_H
