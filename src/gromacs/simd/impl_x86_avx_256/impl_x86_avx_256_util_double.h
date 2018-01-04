/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
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

#include "gromacs/utility/basedefinitions.h"

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
gatherLoadTranspose(const double *        base,
                    const std::int32_t    offset[],
                    SimdDouble *          v0,
                    SimdDouble *          v1,
                    SimdDouble *          v2,
                    SimdDouble *          v3)
{
    assert(std::size_t(offset) % 16 == 0);
    assert(std::size_t(base) % 32 == 0);
    assert(align % 4 == 0);

    v0->simdInternal_ = _mm256_load_pd( base + align * offset[0] );
    v1->simdInternal_ = _mm256_load_pd( base + align * offset[1] );
    v2->simdInternal_ = _mm256_load_pd( base + align * offset[2] );
    v3->simdInternal_ = _mm256_load_pd( base + align * offset[3] );
    avx256Transpose4By4(&v0->simdInternal_, &v1->simdInternal_, &v2->simdInternal_, &v3->simdInternal_);
}

template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const double *        base,
                    const std::int32_t    offset[],
                    SimdDouble *          v0,
                    SimdDouble *          v1)
{
    __m128d t1, t2, t3, t4;
    __m256d tA, tB;

    assert(std::size_t(offset) % 16 == 0);
    assert(std::size_t(base) % 16 == 0);
    assert(align % 2 == 0);

    t1   = _mm_load_pd( base + align * offset[0] );
    t2   = _mm_load_pd( base + align * offset[1] );
    t3   = _mm_load_pd( base + align * offset[2] );
    t4   = _mm_load_pd( base + align * offset[3] );
    tA   = _mm256_insertf128_pd(_mm256_castpd128_pd256(t1), t3, 0x1);
    tB   = _mm256_insertf128_pd(_mm256_castpd128_pd256(t2), t4, 0x1);

    v0->simdInternal_ = _mm256_unpacklo_pd(tA, tB);
    v1->simdInternal_ = _mm256_unpackhi_pd(tA, tB);
}

static const int c_simdBestPairAlignmentDouble = 2;

template <int align>
static inline void gmx_simdcall
gatherLoadUTranspose(const double *        base,
                     const std::int32_t    offset[],
                     SimdDouble *          v0,
                     SimdDouble *          v1,
                     SimdDouble *          v2)
{
    assert(std::size_t(offset) % 16 == 0);

    __m256d t1, t2, t3, t4, t5, t6, t7, t8;
    if (align % 4 == 0)
    {
        t1                = _mm256_load_pd(base + align * offset[0]);
        t2                = _mm256_load_pd(base + align * offset[1]);
        t3                = _mm256_load_pd(base + align * offset[2]);
        t4                = _mm256_load_pd(base + align * offset[3]);
    }
    else
    {
        t1                = _mm256_loadu_pd(base + align * offset[0]);
        t2                = _mm256_loadu_pd(base + align * offset[1]);
        t3                = _mm256_loadu_pd(base + align * offset[2]);
        t4                = _mm256_loadu_pd(base + align * offset[3]);
    }
    t5                = _mm256_unpacklo_pd(t1, t2);
    t6                = _mm256_unpackhi_pd(t1, t2);
    t7                = _mm256_unpacklo_pd(t3, t4);
    t8                = _mm256_unpackhi_pd(t3, t4);
    v0->simdInternal_ = _mm256_permute2f128_pd(t5, t7, 0x20);
    v1->simdInternal_ = _mm256_permute2f128_pd(t6, t8, 0x20);
    v2->simdInternal_ = _mm256_permute2f128_pd(t5, t7, 0x31);
}

template <int align>
static inline void gmx_simdcall
transposeScatterStoreU(double *            base,
                       const std::int32_t  offset[],
                       SimdDouble          v0,
                       SimdDouble          v1,
                       SimdDouble          v2)
{
    __m256d t0, t1, t2;


    assert(std::size_t(offset) % 16 == 0);

    // v0: x0 x1 | x2 x3
    // v1: y0 y1 | y2 y3
    // v2: z0 z1 | z2 z3

    t0 = _mm256_unpacklo_pd(v0.simdInternal_, v1.simdInternal_); // x0 y0 | x2 y2
    t1 = _mm256_unpackhi_pd(v0.simdInternal_, v1.simdInternal_); // x1 y1 | x3 y3
    t2 = _mm256_unpackhi_pd(v2.simdInternal_, v2.simdInternal_); // z1 z1 | z3 z3

    _mm_storeu_pd(base + align * offset[0], _mm256_castpd256_pd128(t0));
    _mm_storeu_pd(base + align * offset[1], _mm256_castpd256_pd128(t1));
    _mm_storeu_pd(base + align * offset[2], _mm256_extractf128_pd(t0, 0x1));
    _mm_storeu_pd(base + align * offset[3], _mm256_extractf128_pd(t1, 0x1));
    _mm_store_sd(base + align * offset[0] + 2, _mm256_castpd256_pd128(v2.simdInternal_));
    _mm_store_sd(base + align * offset[1] + 2, _mm256_castpd256_pd128(t2));
    _mm_store_sd(base + align * offset[2] + 2, _mm256_extractf128_pd(v2.simdInternal_, 0x1));
    _mm_store_sd(base + align * offset[3] + 2, _mm256_extractf128_pd(t2, 0x1));
}

template <int align>
static inline void gmx_simdcall
transposeScatterIncrU(double *            base,
                      const std::int32_t  offset[],
                      SimdDouble          v0,
                      SimdDouble          v1,
                      SimdDouble          v2)
{
    __m256d t0, t1;
    __m128d t2, tA, tB;

    assert(std::size_t(offset) % 16 == 0);

    if (align % 4 == 0)
    {
        // we can use aligned load/store
        t0 = _mm256_setzero_pd();
        avx256Transpose4By4(&v0.simdInternal_, &v1.simdInternal_, &v2.simdInternal_, &t0);
        _mm256_store_pd(base + align * offset[0], _mm256_add_pd(_mm256_load_pd(base + align * offset[0]), v0.simdInternal_));
        _mm256_store_pd(base + align * offset[1], _mm256_add_pd(_mm256_load_pd(base + align * offset[1]), v1.simdInternal_));
        _mm256_store_pd(base + align * offset[2], _mm256_add_pd(_mm256_load_pd(base + align * offset[2]), v2.simdInternal_));
        _mm256_store_pd(base + align * offset[3], _mm256_add_pd(_mm256_load_pd(base + align * offset[3]), t0));
    }
    else
    {
        // v0: x0 x1 | x2 x3
        // v1: y0 y1 | y2 y3
        // v2: z0 z1 | z2 z3

        t0 = _mm256_unpacklo_pd(v0.simdInternal_, v1.simdInternal_); // x0 y0 | x2 y2
        t1 = _mm256_unpackhi_pd(v0.simdInternal_, v1.simdInternal_); // x1 y1 | x3 y3
        t2 = _mm256_extractf128_pd(v2.simdInternal_, 0x1);           // z2 z3

        tA = _mm_loadu_pd(base + align * offset[0]);
        tB = _mm_load_sd(base + align * offset[0] + 2);
        tA = _mm_add_pd(tA, _mm256_castpd256_pd128(t0));
        tB = _mm_add_pd(tB, _mm256_castpd256_pd128(v2.simdInternal_));
        _mm_storeu_pd(base + align * offset[0], tA);
        _mm_store_sd(base + align * offset[0] + 2, tB);

        tA = _mm_loadu_pd(base + align * offset[1]);
        tB = _mm_loadh_pd(_mm_setzero_pd(), base + align * offset[1] + 2);
        tA = _mm_add_pd(tA, _mm256_castpd256_pd128(t1));
        tB = _mm_add_pd(tB, _mm256_castpd256_pd128(v2.simdInternal_));
        _mm_storeu_pd(base + align * offset[1], tA);
        _mm_storeh_pd(base + align * offset[1] + 2, tB);

        tA = _mm_loadu_pd(base + align * offset[2]);
        tB = _mm_load_sd(base + align * offset[2] + 2);
        tA = _mm_add_pd(tA, _mm256_extractf128_pd(t0, 0x1));
        tB = _mm_add_pd(tB, t2);
        _mm_storeu_pd(base + align * offset[2], tA);
        _mm_store_sd(base + align * offset[2] + 2, tB);

        tA = _mm_loadu_pd(base + align * offset[3]);
        tB = _mm_loadh_pd(_mm_setzero_pd(), base + align * offset[3] + 2);
        tA = _mm_add_pd(tA, _mm256_extractf128_pd(t1, 0x1));
        tB = _mm_add_pd(tB, t2);
        _mm_storeu_pd(base + align * offset[3], tA);
        _mm_storeh_pd(base + align * offset[3] + 2, tB);
    }
}
template <int align>
static inline void gmx_simdcall
transposeScatterDecrU(double *            base,
                      const std::int32_t  offset[],
                      SimdDouble          v0,
                      SimdDouble          v1,
                      SimdDouble          v2)
{
    __m256d t0, t1;
    __m128d t2, tA, tB;

    assert(std::size_t(offset) % 16 == 0);

    if (align % 4 == 0)
    {
        // we can use aligned load/store
        t0 = _mm256_setzero_pd();
        avx256Transpose4By4(&v0.simdInternal_, &v1.simdInternal_, &v2.simdInternal_, &t0);
        _mm256_store_pd(base + align * offset[0], _mm256_sub_pd(_mm256_load_pd(base + align * offset[0]), v0.simdInternal_));
        _mm256_store_pd(base + align * offset[1], _mm256_sub_pd(_mm256_load_pd(base + align * offset[1]), v1.simdInternal_));
        _mm256_store_pd(base + align * offset[2], _mm256_sub_pd(_mm256_load_pd(base + align * offset[2]), v2.simdInternal_));
        _mm256_store_pd(base + align * offset[3], _mm256_sub_pd(_mm256_load_pd(base + align * offset[3]), t0));
    }
    else
    {
        // v0: x0 x1 | x2 x3
        // v1: y0 y1 | y2 y3
        // v2: z0 z1 | z2 z3

        t0 = _mm256_unpacklo_pd(v0.simdInternal_, v1.simdInternal_); // x0 y0 | x2 y2
        t1 = _mm256_unpackhi_pd(v0.simdInternal_, v1.simdInternal_); // x1 y1 | x3 y3
        t2 = _mm256_extractf128_pd(v2.simdInternal_, 0x1);           // z2 z3

        tA = _mm_loadu_pd(base + align * offset[0]);
        tB = _mm_load_sd(base + align * offset[0] + 2);
        tA = _mm_sub_pd(tA, _mm256_castpd256_pd128(t0));
        tB = _mm_sub_pd(tB, _mm256_castpd256_pd128(v2.simdInternal_));
        _mm_storeu_pd(base + align * offset[0], tA);
        _mm_store_sd(base + align * offset[0] + 2, tB);

        tA = _mm_loadu_pd(base + align * offset[1]);
        tB = _mm_loadh_pd(_mm_setzero_pd(), base + align * offset[1] + 2);
        tA = _mm_sub_pd(tA, _mm256_castpd256_pd128(t1));
        tB = _mm_sub_pd(tB, _mm256_castpd256_pd128(v2.simdInternal_));
        _mm_storeu_pd(base + align * offset[1], tA);
        _mm_storeh_pd(base + align * offset[1] + 2, tB);

        tA = _mm_loadu_pd(base + align * offset[2]);
        tB = _mm_load_sd(base + align * offset[2] + 2);
        tA = _mm_sub_pd(tA, _mm256_extractf128_pd(t0, 0x1));
        tB = _mm_sub_pd(tB, t2);
        _mm_storeu_pd(base + align * offset[2], tA);
        _mm_store_sd(base + align * offset[2] + 2, tB);

        tA = _mm_loadu_pd(base + align * offset[3]);
        tB = _mm_loadh_pd(_mm_setzero_pd(), base + align * offset[3] + 2);
        tA = _mm_sub_pd(tA, _mm256_extractf128_pd(t1, 0x1));
        tB = _mm_sub_pd(tB, t2);
        _mm_storeu_pd(base + align * offset[3], tA);
        _mm_storeh_pd(base + align * offset[3] + 2, tB);
    }
}

static inline void gmx_simdcall
expandScalarsToTriplets(SimdDouble    scalar,
                        SimdDouble *  triplets0,
                        SimdDouble *  triplets1,
                        SimdDouble *  triplets2)
{
    __m256d t0 = _mm256_permute2f128_pd(scalar.simdInternal_, scalar.simdInternal_, 0x21);
    __m256d t1 = _mm256_permute_pd(scalar.simdInternal_, 0b0000);
    __m256d t2 = _mm256_permute_pd(scalar.simdInternal_, 0b1111);
    triplets0->simdInternal_ = _mm256_blend_pd(t1, t0, 0b1100);
    triplets1->simdInternal_ = _mm256_blend_pd(t2, t1, 0b1100);
    triplets2->simdInternal_ = _mm256_blend_pd(t0, t2, 0b1100);
}

template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const double *  base,
                             SimdDInt32      offset,
                             SimdDouble *    v0,
                             SimdDouble *    v1,
                             SimdDouble *    v2,
                             SimdDouble *    v3)
{
    assert(std::size_t(base) % 32 == 0);
    assert(align % 4 == 0);

    alignas(GMX_SIMD_ALIGNMENT) std::int32_t ioffset[GMX_SIMD_DINT32_WIDTH];
    _mm_store_si128( reinterpret_cast<__m128i *>(ioffset), offset.simdInternal_);

    v0->simdInternal_ = _mm256_load_pd(base + align * ioffset[0]);
    v1->simdInternal_ = _mm256_load_pd(base + align * ioffset[1]);
    v2->simdInternal_ = _mm256_load_pd(base + align * ioffset[2]);
    v3->simdInternal_ = _mm256_load_pd(base + align * ioffset[3]);

    avx256Transpose4By4(&v0->simdInternal_, &v1->simdInternal_, &v2->simdInternal_, &v3->simdInternal_);
}

template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const double *    base,
                             SimdDInt32        offset,
                             SimdDouble *      v0,
                             SimdDouble *      v1)
{
    __m128d t1, t2, t3, t4;
    __m256d tA, tB;

    assert(std::size_t(base) % 16 == 0);
    assert(align % 2 == 0);

    alignas(GMX_SIMD_ALIGNMENT) std::int32_t  ioffset[GMX_SIMD_DINT32_WIDTH];
    _mm_store_si128( reinterpret_cast<__m128i *>(ioffset), offset.simdInternal_);

    t1  = _mm_load_pd(base + align * ioffset[0]);
    t2  = _mm_load_pd(base + align * ioffset[1]);
    t3  = _mm_load_pd(base + align * ioffset[2]);
    t4  = _mm_load_pd(base + align * ioffset[3]);

    tA                = _mm256_insertf128_pd(_mm256_castpd128_pd256(t1), t3, 0x1);
    tB                = _mm256_insertf128_pd(_mm256_castpd128_pd256(t2), t4, 0x1);
    v0->simdInternal_ = _mm256_unpacklo_pd(tA, tB);
    v1->simdInternal_ = _mm256_unpackhi_pd(tA, tB);
}

template <int align>
static inline void gmx_simdcall
gatherLoadUBySimdIntTranspose(const double *  base,
                              SimdDInt32      offset,
                              SimdDouble *    v0,
                              SimdDouble *    v1)
{
    __m128d t1, t2, t3, t4;
    __m256d tA, tB;

    alignas(GMX_SIMD_ALIGNMENT) std::int32_t ioffset[GMX_SIMD_DINT32_WIDTH];
    _mm_store_si128( reinterpret_cast<__m128i *>(ioffset), offset.simdInternal_);

    t1   = _mm_loadu_pd(base + align * ioffset[0]);
    t2   = _mm_loadu_pd(base + align * ioffset[1]);
    t3   = _mm_loadu_pd(base + align * ioffset[2]);
    t4   = _mm_loadu_pd(base + align * ioffset[3]);

    tA  = _mm256_insertf128_pd(_mm256_castpd128_pd256(t1), t3, 0x1);
    tB  = _mm256_insertf128_pd(_mm256_castpd128_pd256(t2), t4, 0x1);

    v0->simdInternal_ = _mm256_unpacklo_pd(tA, tB);
    v1->simdInternal_ = _mm256_unpackhi_pd(tA, tB);
}

static inline double gmx_simdcall
reduceIncr4ReturnSum(double *    m,
                     SimdDouble  v0,
                     SimdDouble  v1,
                     SimdDouble  v2,
                     SimdDouble  v3)
{
    __m256d t0, t1, t2;
    __m128d a0, a1;

    assert(std::size_t(m) % 32 == 0);

    t0 = _mm256_hadd_pd(v0.simdInternal_, v1.simdInternal_);
    t1 = _mm256_hadd_pd(v2.simdInternal_, v3.simdInternal_);
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
