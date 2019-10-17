/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2019, by the GROMACS development team, led by
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
#ifndef GMX_SIMD_IMPL_X86_SSE2_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_X86_SSE2_UTIL_DOUBLE_H

#include "config.h"

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <emmintrin.h>

#include "impl_x86_sse2_simd_double.h"

namespace gmx
{

template<int align>
static inline void gmx_simdcall gatherLoadTranspose(const double*      base,
                                                    const std::int32_t offset[],
                                                    SimdDouble*        v0,
                                                    SimdDouble*        v1,
                                                    SimdDouble*        v2,
                                                    SimdDouble*        v3)
{
    __m128d t1, t2, t3, t4;

    assert(std::size_t(base + align * offset[0]) % 16 == 0);
    assert(std::size_t(base + align * offset[1]) % 16 == 0);

    t1                = _mm_load_pd(base + align * offset[0]);
    t2                = _mm_load_pd(base + align * offset[1]);
    t3                = _mm_load_pd(base + align * offset[0] + 2);
    t4                = _mm_load_pd(base + align * offset[1] + 2);
    v0->simdInternal_ = _mm_unpacklo_pd(t1, t2);
    v1->simdInternal_ = _mm_unpackhi_pd(t1, t2);
    v2->simdInternal_ = _mm_unpacklo_pd(t3, t4);
    v3->simdInternal_ = _mm_unpackhi_pd(t3, t4);
}

template<int align>
static inline void gmx_simdcall
                   gatherLoadTranspose(const double* base, const std::int32_t offset[], SimdDouble* v0, SimdDouble* v1)
{
    __m128d t1, t2;

    assert(std::size_t(base + align * offset[0]) % 16 == 0);
    assert(std::size_t(base + align * offset[1]) % 16 == 0);

    t1                = _mm_load_pd(base + align * offset[0]);
    t2                = _mm_load_pd(base + align * offset[1]);
    v0->simdInternal_ = _mm_unpacklo_pd(t1, t2);
    v1->simdInternal_ = _mm_unpackhi_pd(t1, t2);
}

static const int c_simdBestPairAlignmentDouble = 2;

template<int align>
static inline void gmx_simdcall gatherLoadUTranspose(const double*      base,
                                                     const std::int32_t offset[],
                                                     SimdDouble*        v0,
                                                     SimdDouble*        v1,
                                                     SimdDouble*        v2)
{
    __m128d t1, t2, t3, t4;
    t1                = _mm_loadu_pd(base + align * offset[0]);
    t2                = _mm_loadu_pd(base + align * offset[1]);
    t3                = _mm_load_sd(base + align * offset[0] + 2);
    t4                = _mm_load_sd(base + align * offset[1] + 2);
    v0->simdInternal_ = _mm_unpacklo_pd(t1, t2);
    v1->simdInternal_ = _mm_unpackhi_pd(t1, t2);
    v2->simdInternal_ = _mm_unpacklo_pd(t3, t4);
}

template<int align>
static inline void gmx_simdcall transposeScatterStoreU(double*            base,
                                                       const std::int32_t offset[],
                                                       SimdDouble         v0,
                                                       SimdDouble         v1,
                                                       SimdDouble         v2)
{
    __m128d t1, t2;
    t1 = _mm_unpacklo_pd(v0.simdInternal_, v1.simdInternal_);
    t2 = _mm_unpackhi_pd(v0.simdInternal_, v1.simdInternal_);
    _mm_storeu_pd(base + align * offset[0], t1);
    _mm_store_sd(base + align * offset[0] + 2, v2.simdInternal_);
    _mm_storeu_pd(base + align * offset[1], t2);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + align * offset[1] + 2), _mm_castpd_ps(v2.simdInternal_));
}

template<int align>
static inline void gmx_simdcall
                   transposeScatterIncrU(double* base, const std::int32_t offset[], SimdDouble v0, SimdDouble v1, SimdDouble v2)
{
    __m128d t1, t2, t3, t4, t5, t6, t7;

    t5 = _mm_unpacklo_pd(v0.simdInternal_, v1.simdInternal_);
    t6 = _mm_unpackhi_pd(v0.simdInternal_, v1.simdInternal_);
    t7 = _mm_unpackhi_pd(v2.simdInternal_, v2.simdInternal_);

    t1 = _mm_loadu_pd(base + align * offset[0]);
    t2 = _mm_load_sd(base + align * offset[0] + 2);
    t1 = _mm_add_pd(t1, t5);
    t2 = _mm_add_sd(t2, v2.simdInternal_);
    _mm_storeu_pd(base + align * offset[0], t1);
    _mm_store_sd(base + align * offset[0] + 2, t2);

    t3 = _mm_loadu_pd(base + align * offset[1]);
    t4 = _mm_load_sd(base + align * offset[1] + 2);
    t3 = _mm_add_pd(t3, t6);
    t4 = _mm_add_sd(t4, t7);
    _mm_storeu_pd(base + align * offset[1], t3);
    _mm_store_sd(base + align * offset[1] + 2, t4);
}

template<int align>
static inline void gmx_simdcall
                   transposeScatterDecrU(double* base, const std::int32_t offset[], SimdDouble v0, SimdDouble v1, SimdDouble v2)
{
    // This implementation is identical to the increment version, apart from using subtraction instead
    __m128d t1, t2, t3, t4, t5, t6, t7;

    t5 = _mm_unpacklo_pd(v0.simdInternal_, v1.simdInternal_);
    t6 = _mm_unpackhi_pd(v0.simdInternal_, v1.simdInternal_);
    t7 = _mm_unpackhi_pd(v2.simdInternal_, v2.simdInternal_);

    t1 = _mm_loadu_pd(base + align * offset[0]);
    t2 = _mm_load_sd(base + align * offset[0] + 2);
    t1 = _mm_sub_pd(t1, t5);
    t2 = _mm_sub_sd(t2, v2.simdInternal_);
    _mm_storeu_pd(base + align * offset[0], t1);
    _mm_store_sd(base + align * offset[0] + 2, t2);

    t3 = _mm_loadu_pd(base + align * offset[1]);
    t4 = _mm_load_sd(base + align * offset[1] + 2);
    t3 = _mm_sub_pd(t3, t6);
    t4 = _mm_sub_sd(t4, t7);
    _mm_storeu_pd(base + align * offset[1], t3);
    _mm_store_sd(base + align * offset[1] + 2, t4);
}

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline void gmx_simdcall expandScalarsToTriplets(SimdDouble  scalar,
                                                        SimdDouble* triplets0,
                                                        SimdDouble* triplets1,
                                                        SimdDouble* triplets2)
{
    triplets0->simdInternal_ =
            _mm_shuffle_pd(scalar.simdInternal_, scalar.simdInternal_, _MM_SHUFFLE2(0, 0));
    triplets1->simdInternal_ =
            _mm_shuffle_pd(scalar.simdInternal_, scalar.simdInternal_, _MM_SHUFFLE2(1, 0));
    triplets2->simdInternal_ =
            _mm_shuffle_pd(scalar.simdInternal_, scalar.simdInternal_, _MM_SHUFFLE2(1, 1));
}
#endif


template<int align>
static inline void gmx_simdcall gatherLoadBySimdIntTranspose(const double* base,
                                                             SimdDInt32    offset,
                                                             SimdDouble*   v0,
                                                             SimdDouble*   v1,
                                                             SimdDouble*   v2,
                                                             SimdDouble*   v3)
{
    __m128d t1, t2, t3, t4;
    // Use optimized bit-shift multiply for the most common alignments
    if (align == 4)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 2);
    }
    else if (align == 8)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset.simdInternal_ =
                _mm_add_epi32(offset.simdInternal_, _mm_slli_epi32(offset.simdInternal_, 1));
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 2);
    }
    else if (align == 16)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 4);
    }

    if (align == 4 || align == 8 || align == 12 || align == 16)
    {
        assert(std::size_t(base + extract<0>(offset)) % 16 == 0);
        assert(std::size_t(base + extract<1>(offset)) % 16 == 0);

        t1 = _mm_load_pd(base + extract<0>(offset));
        t2 = _mm_load_pd(base + extract<1>(offset));
        t3 = _mm_load_pd(base + extract<0>(offset) + 2);
        t4 = _mm_load_pd(base + extract<1>(offset) + 2);
    }
    else
    {
        assert(std::size_t(base + align * extract<0>(offset)) % 16 == 0);
        assert(std::size_t(base + align * extract<1>(offset)) % 16 == 0);

        t1 = _mm_load_pd(base + align * extract<0>(offset));
        t2 = _mm_load_pd(base + align * extract<1>(offset));
        t3 = _mm_load_pd(base + align * extract<0>(offset) + 2);
        t4 = _mm_load_pd(base + align * extract<1>(offset) + 2);
    }
    v0->simdInternal_ = _mm_unpacklo_pd(t1, t2);
    v1->simdInternal_ = _mm_unpackhi_pd(t1, t2);
    v2->simdInternal_ = _mm_unpacklo_pd(t3, t4);
    v3->simdInternal_ = _mm_unpackhi_pd(t3, t4);
}

template<int align>
static inline void gmx_simdcall
                   gatherLoadBySimdIntTranspose(const double* base, SimdDInt32 offset, SimdDouble* v0, SimdDouble* v1)
{
    __m128d t1, t2;

    // Use optimized bit-shift multiply for the most common alignments
    if (align == 2)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 1);
    }
    else if (align == 4)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 2);
    }
    else if (align == 6)
    {
        // multiply by 3, then by 2
        offset.simdInternal_ =
                _mm_add_epi32(offset.simdInternal_, _mm_slli_epi32(offset.simdInternal_, 1));
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 1);
    }
    else if (align == 8)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 3);
    }
    else if (align == 12)
    {
        // multiply by 3, then by 4
        offset.simdInternal_ =
                _mm_add_epi32(offset.simdInternal_, _mm_slli_epi32(offset.simdInternal_, 1));
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 2);
    }
    else if (align == 16)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 4);
    }

    if (align == 2 || align == 4 || align == 6 || align == 8 || align == 12 || align == 16)
    {
        assert(std::size_t(base + extract<0>(offset)) % 16 == 0);
        assert(std::size_t(base + extract<1>(offset)) % 16 == 0);

        t1 = _mm_load_pd(base + extract<0>(offset));
        t2 = _mm_load_pd(base + extract<1>(offset));
    }
    else
    {
        assert(std::size_t(base + align * extract<0>(offset)) % 16 == 0);
        assert(std::size_t(base + align * extract<1>(offset)) % 16 == 0);

        t1 = _mm_load_pd(base + align * extract<0>(offset));
        t2 = _mm_load_pd(base + align * extract<1>(offset));
    }
    v0->simdInternal_ = _mm_unpacklo_pd(t1, t2);
    v1->simdInternal_ = _mm_unpackhi_pd(t1, t2);
}


template<int align>
static inline void gmx_simdcall
                   gatherLoadUBySimdIntTranspose(const double* base, SimdDInt32 offset, SimdDouble* v0, SimdDouble* v1)
{
    __m128d t1, t2;
    // Use optimized bit-shift multiply for the most common alignments.

    // Do nothing for align == 1
    if (align == 2)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 1);
    }
    else if (align == 4)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 2);
    }

    if (align == 1 || align == 2 || align == 4)
    {
        t1 = _mm_loadu_pd(base + extract<0>(offset));
        t2 = _mm_loadu_pd(base + extract<1>(offset));
    }
    else
    {
        t1 = _mm_loadu_pd(base + align * extract<0>(offset));
        t2 = _mm_loadu_pd(base + align * extract<1>(offset));
    }
    v0->simdInternal_ = _mm_unpacklo_pd(t1, t2);
    v1->simdInternal_ = _mm_unpackhi_pd(t1, t2);
}

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline double gmx_simdcall
                     reduceIncr4ReturnSum(double* m, SimdDouble v0, SimdDouble v1, SimdDouble v2, SimdDouble v3)
{
    __m128d t1, t2, t3, t4;

    t1 = _mm_unpacklo_pd(v0.simdInternal_, v1.simdInternal_);
    t2 = _mm_unpackhi_pd(v0.simdInternal_, v1.simdInternal_);
    t3 = _mm_unpacklo_pd(v2.simdInternal_, v3.simdInternal_);
    t4 = _mm_unpackhi_pd(v2.simdInternal_, v3.simdInternal_);

    t1 = _mm_add_pd(t1, t2);
    t3 = _mm_add_pd(t3, t4);

    t2 = _mm_add_pd(t1, _mm_load_pd(m));
    t4 = _mm_add_pd(t3, _mm_load_pd(m + 2));

    assert(std::size_t(m) % 16 == 0);

    _mm_store_pd(m, t2);
    _mm_store_pd(m + 2, t4);

    t1 = _mm_add_pd(t1, t3);

    t2 = _mm_add_sd(t1, _mm_shuffle_pd(t1, t1, _MM_SHUFFLE2(1, 1)));
    return *reinterpret_cast<double*>(&t2);
}
#endif

} // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE2_UTIL_DOUBLE_H
