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
#ifndef GMX_SIMD_IMPL_X86_SSE2_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_X86_SSE2_UTIL_FLOAT_H

#include "config.h"

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <emmintrin.h>

#include "impl_x86_sse2_simd_float.h"

namespace gmx
{

template <int align>
static inline void gmx_simdcall
simdGatherLoadTransposeF(const float *        base,
                         const std::int32_t   offset[],
                         SimdFloat *          v0,
                         SimdFloat *          v1,
                         SimdFloat *          v2,
                         SimdFloat *          v3)
{
    assert(std::size_t(base + align * offset[0]) % 16 == 0);
    assert(std::size_t(base + align * offset[1]) % 16 == 0);
    assert(std::size_t(base + align * offset[2]) % 16 == 0);
    assert(std::size_t(base + align * offset[3]) % 16 == 0);

    v0->r = _mm_load_ps( base + align * offset[0] );
    v1->r = _mm_load_ps( base + align * offset[1] );
    v2->r = _mm_load_ps( base + align * offset[2] );
    v3->r = _mm_load_ps( base + align * offset[3] );

    _MM_TRANSPOSE4_PS(v0->r, v1->r, v2->r, v3->r);
}

template <int align>
static inline void gmx_simdcall
simdGatherLoadTransposeF(const float *        base,
                         const std::int32_t   offset[],
                         SimdFloat *          v0,
                         SimdFloat *          v1)
{
    __m128 t1, t2;

    v0->r = _mm_castpd_ps(_mm_load_sd( reinterpret_cast<const double *>( base + align * offset[0] ) ));
    v1->r = _mm_castpd_ps(_mm_load_sd( reinterpret_cast<const double *>( base + align * offset[1] ) ));
    t1    = _mm_castpd_ps(_mm_load_sd( reinterpret_cast<const double *>( base + align * offset[2] ) ));
    t2    = _mm_castpd_ps(_mm_load_sd( reinterpret_cast<const double *>( base + align * offset[3] ) ));
    t1    = _mm_unpacklo_ps(v0->r, t1);
    t2    = _mm_unpacklo_ps(v1->r, t2);
    v0->r = _mm_unpacklo_ps(t1, t2);
    v1->r = _mm_unpackhi_ps(t1, t2);
}

static const int c_simdBestPairAlignmentF = 2;

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
template <int align>
static inline void gmx_simdcall
simdGatherLoadUTransposeF(const float *        base,
                          const std::int32_t   offset[],
                          SimdFloat *          v0,
                          SimdFloat *          v1,
                          SimdFloat *          v2)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;

    // The template conditional should be optimized away at compile time
    if ( (align & 0x3) == 0)
    {
        // With alignment 4 or better we can read a byte beyond triplets and use _mm_loadu_ps()
        t1    = _mm_loadu_ps( base + align * offset[0] );
        t2    = _mm_loadu_ps( base + align * offset[1] );
        t3    = _mm_loadu_ps( base + align * offset[2] );
        t4    = _mm_loadu_ps( base + align * offset[3] );
        t5    = _mm_unpacklo_ps(t1, t2);
        t6    = _mm_unpacklo_ps(t3, t4);
        t7    = _mm_unpackhi_ps(t1, t2);
        t8    = _mm_unpackhi_ps(t3, t4);
        v0->r = _mm_movelh_ps(t5, t6);
        v1->r = _mm_movehl_ps(t6, t5);
        v2->r = _mm_movelh_ps(t7, t8);
    }
    else
    {
        t1    = _mm_castpd_ps(_mm_load_sd( reinterpret_cast<const double *>( base + align * offset[0] ) ));
        t2    = _mm_castpd_ps(_mm_load_sd( reinterpret_cast<const double *>( base + align * offset[1] ) ));
        t3    = _mm_castpd_ps(_mm_load_sd( reinterpret_cast<const double *>( base + align * offset[2] ) ));
        t4    = _mm_castpd_ps(_mm_load_sd( reinterpret_cast<const double *>( base + align * offset[3] ) ));
        t5    = _mm_load_ss( base + align * offset[0] + 2 );
        t6    = _mm_load_ss( base + align * offset[1] + 2 );
        t7    = _mm_load_ss( base + align * offset[2] + 2 );
        t8    = _mm_load_ss( base + align * offset[3] + 2 );
        t1    = _mm_unpacklo_ps(t1, t2);
        t3    = _mm_unpacklo_ps(t3, t4);
        v0->r = _mm_movelh_ps(t1, t3);
        v1->r = _mm_movehl_ps(t3, t1);
        t5    = _mm_unpacklo_ps(t5, t6);
        t7    = _mm_unpacklo_ps(t7, t8);
        v2->r = _mm_movelh_ps(t5, t7);
    }
}
#endif

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
template <int align>
static inline void gmx_simdcall
simdTransposeScatterStoreUF(float *              base,
                            const std::int32_t   offset[],
                            SimdFloat            v0,
                            SimdFloat            v1,
                            SimdFloat            v2)
{
    __m128 t1, t2;

    t1   = _mm_unpacklo_ps(v0.r, v1.r);
    t2   = _mm_unpackhi_ps(v0.r, v1.r);
    _mm_storel_pi( reinterpret_cast< __m64 *>( base + align * offset[0] ), t1);
    _mm_store_ss(base + align * offset[0] + 2, v2.r);
    _mm_storeh_pi( reinterpret_cast< __m64 *>( base + align * offset[1] ), t1);
    _mm_store_ss(base + align * offset[1] + 2, _mm_shuffle_ps(v2.r, v2.r, _MM_SHUFFLE(1, 1, 1, 1)));
    _mm_storel_pi( reinterpret_cast< __m64 *>( base + align * offset[2] ), t2);
    _mm_store_ss(base + align * offset[2] + 2, _mm_shuffle_ps(v2.r, v2.r, _MM_SHUFFLE(2, 2, 2, 2)));
    _mm_storeh_pi( reinterpret_cast< __m64 *>( base + align * offset[3] ), t2);
    _mm_store_ss(base + align * offset[3] + 2, _mm_shuffle_ps(v2.r, v2.r, _MM_SHUFFLE(3, 3, 3, 3)));
}
#endif

template <int align>
static inline void gmx_simdcall
simdTransposeScatterIncrUF(float *              base,
                           const std::int32_t   offset[],
                           SimdFloat            v0,
                           SimdFloat            v1,
                           SimdFloat            v2)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    t5          = _mm_unpacklo_ps(v1.r, v2.r);
    t6          = _mm_unpackhi_ps(v1.r, v2.r);
    t7          = _mm_shuffle_ps(v0.r, t5, _MM_SHUFFLE(1, 0, 0, 0));
    t8          = _mm_shuffle_ps(v0.r, t5, _MM_SHUFFLE(3, 2, 0, 1));
    t9          = _mm_shuffle_ps(v0.r, t6, _MM_SHUFFLE(1, 0, 0, 2));
    t10         = _mm_shuffle_ps(v0.r, t6, _MM_SHUFFLE(3, 2, 0, 3));
    t1          = _mm_load_ss(base + align * offset[0]);
    t1          = _mm_loadh_pi(t1, reinterpret_cast< __m64 *>(base + align * offset[0] + 1));
    t1          = _mm_add_ps(t1, t7);
    _mm_store_ss(base + align * offset[0], t1);
    _mm_storeh_pi(reinterpret_cast< __m64 *>(base + align * offset[0] + 1), t1);
    t2          = _mm_load_ss(base + align * offset[1]);
    t2          = _mm_loadh_pi(t2, reinterpret_cast< __m64 *>(base + align * offset[1] + 1));
    t2          = _mm_add_ps(t2, t8);
    _mm_store_ss(base + align * offset[1], t2);
    _mm_storeh_pi(reinterpret_cast< __m64 *>(base + align * offset[1] + 1), t2);
    t3          = _mm_load_ss(base + align * offset[2]);
    t3          = _mm_loadh_pi(t3, reinterpret_cast< __m64 *>(base + align * offset[2] + 1));
    t3          = _mm_add_ps(t3, t9);
    _mm_store_ss(base + align * offset[2], t3);
    _mm_storeh_pi(reinterpret_cast< __m64 *>(base + align * offset[2] + 1), t3);
    t4          = _mm_load_ss(base + align * offset[3]);
    t4          = _mm_loadh_pi(t4, reinterpret_cast< __m64 *>(base + align * offset[3] + 1));
    t4          = _mm_add_ps(t4, t10);
    _mm_store_ss(base + align * offset[3], t4);
    _mm_storeh_pi(reinterpret_cast< __m64 *>(base + align * offset[3] + 1), t4);
}

template <int align>
static inline void gmx_simdcall
simdTransposeScatterDecrUF(float *              base,
                           const std::int32_t   offset[],
                           SimdFloat            v0,
                           SimdFloat            v1,
                           SimdFloat            v2)
{
    // This implementation is identical to the increment version, apart from using subtraction instead
    __m128 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    t5          = _mm_unpacklo_ps(v1.r, v2.r);
    t6          = _mm_unpackhi_ps(v1.r, v2.r);
    t7          = _mm_shuffle_ps(v0.r, t5, _MM_SHUFFLE(1, 0, 0, 0));
    t8          = _mm_shuffle_ps(v0.r, t5, _MM_SHUFFLE(3, 2, 0, 1));
    t9          = _mm_shuffle_ps(v0.r, t6, _MM_SHUFFLE(1, 0, 0, 2));
    t10         = _mm_shuffle_ps(v0.r, t6, _MM_SHUFFLE(3, 2, 0, 3));
    t1          = _mm_load_ss(base + align * offset[0]);
    t1          = _mm_loadh_pi(t1, reinterpret_cast< __m64 *>(base + align * offset[0] + 1));
    t1          = _mm_sub_ps(t1, t7);
    _mm_store_ss(base + align * offset[0], t1);
    _mm_storeh_pi(reinterpret_cast< __m64 *>(base + align * offset[0] + 1), t1);
    t2          = _mm_load_ss(base + align * offset[1]);
    t2          = _mm_loadh_pi(t2, reinterpret_cast< __m64 *>(base + align * offset[1] + 1));
    t2          = _mm_sub_ps(t2, t8);
    _mm_store_ss(base + align * offset[1], t2);
    _mm_storeh_pi(reinterpret_cast< __m64 *>(base + align * offset[1] + 1), t2);
    t3          = _mm_load_ss(base + align * offset[2]);
    t3          = _mm_loadh_pi(t3, reinterpret_cast< __m64 *>(base + align * offset[2] + 1));
    t3          = _mm_sub_ps(t3, t9);
    _mm_store_ss(base + align * offset[2], t3);
    _mm_storeh_pi(reinterpret_cast< __m64 *>(base + align * offset[2] + 1), t3);
    t4          = _mm_load_ss(base + align * offset[3]);
    t4          = _mm_loadh_pi(t4, reinterpret_cast< __m64 *>(base + align * offset[3] + 1));
    t4          = _mm_sub_ps(t4, t10);
    _mm_store_ss(base + align * offset[3], t4);
    _mm_storeh_pi(reinterpret_cast< __m64 *>(base + align * offset[3] + 1), t4);
}

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline void gmx_simdcall
simdExpandScalarsToTripletsF(SimdFloat    scalar,
                             SimdFloat *  triplets0,
                             SimdFloat *  triplets1,
                             SimdFloat *  triplets2)
{
    triplets0->r = _mm_shuffle_ps(scalar.r, scalar.r, _MM_SHUFFLE(1, 0, 0, 0));
    triplets1->r = _mm_shuffle_ps(scalar.r, scalar.r, _MM_SHUFFLE(2, 2, 1, 1));
    triplets2->r = _mm_shuffle_ps(scalar.r, scalar.r, _MM_SHUFFLE(3, 3, 3, 2));
}
#endif

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
template <int align>
static inline void gmx_simdcall
simdGatherLoadBySimdIntTransposeF(const float *  base,
                                  SimdFInt32     offset,
                                  SimdFloat *    v0,
                                  SimdFloat *    v1,
                                  SimdFloat *    v2,
                                  SimdFloat *    v3)
{
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
        assert(std::size_t(base + simdExtractFI<0>(offset)) % 16 == 0);
        assert(std::size_t(base + simdExtractFI<1>(offset)) % 16 == 0);
        assert(std::size_t(base + simdExtractFI<2>(offset)) % 16 == 0);
        assert(std::size_t(base + simdExtractFI<3>(offset)) % 16 == 0);

        v0->r  = _mm_load_ps(base + simdExtractFI<0>(offset) );
        v1->r  = _mm_load_ps(base + simdExtractFI<1>(offset) );
        v2->r  = _mm_load_ps(base + simdExtractFI<2>(offset) );
        v3->r  = _mm_load_ps(base + simdExtractFI<3>(offset) );
    }
    else
    {
        assert(std::size_t(base + align * simdExtractFI<0>(offset)) % 16 == 0);
        assert(std::size_t(base + align * simdExtractFI<1>(offset)) % 16 == 0);
        assert(std::size_t(base + align * simdExtractFI<2>(offset)) % 16 == 0);
        assert(std::size_t(base + align * simdExtractFI<3>(offset)) % 16 == 0);

        v0->r  = _mm_load_ps(base + align * simdExtractFI<0>(offset) );
        v1->r  = _mm_load_ps(base + align * simdExtractFI<1>(offset) );
        v2->r  = _mm_load_ps(base + align * simdExtractFI<2>(offset) );
        v3->r  = _mm_load_ps(base + align * simdExtractFI<3>(offset) );
    }
    _MM_TRANSPOSE4_PS(v0->r, v1->r, v2->r, v3->r);
}

template <int align>
static inline void gmx_simdcall
simdGatherLoadBySimdIntTransposeF(const float *   base,
                                  SimdFInt32      offset,
                                  SimdFloat *     v0,
                                  SimdFloat *     v1)
{
    __m128 t1, t2;

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
        v0->r  = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + simdExtractFI<0>(offset))));
        v1->r  = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + simdExtractFI<1>(offset))));
        t1     = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + simdExtractFI<2>(offset))));
        t2     = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + simdExtractFI<3>(offset))));
    }
    else
    {
        v0->r  = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + align * simdExtractFI<0>(offset))));
        v1->r  = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + align * simdExtractFI<1>(offset))));
        t1     = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + align * simdExtractFI<2>(offset))));
        t2     = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + align * simdExtractFI<3>(offset))));
    }
    t1     = _mm_unpacklo_ps(v0->r, t1);
    t2     = _mm_unpacklo_ps(v1->r, t2);
    v0->r  = _mm_unpacklo_ps(t1, t2);
    v1->r  = _mm_unpackhi_ps(t1, t2);
}
#endif

template <int align>
static inline void gmx_simdcall
simdGatherLoadUBySimdIntTransposeF(const float *  base,
                                   SimdFInt32     offset,
                                   SimdFloat *    v0,
                                   SimdFloat *    v1)
{
    __m128 t1, t2;

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
        v0->r  = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + simdExtractFI<0>(offset))));
        v1->r  = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + simdExtractFI<1>(offset))));
        t1     = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + simdExtractFI<2>(offset))));
        t2     = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + simdExtractFI<3>(offset))));
    }
    else
    {
        v0->r  = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + align * simdExtractFI<0>(offset))));
        v1->r  = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + align * simdExtractFI<1>(offset))));
        t1     = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + align * simdExtractFI<2>(offset))));
        t2     = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double *>(base + align * simdExtractFI<3>(offset))));
    }
    t1    = _mm_unpacklo_ps(v0->r, t1);
    t2    = _mm_unpacklo_ps(v1->r, t2);
    v0->r = _mm_unpacklo_ps(t1, t2);
    v1->r = _mm_unpackhi_ps(t1, t2);
}

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline float gmx_simdcall
simdReduceIncr4ReturnSumF(float *    m,
                          SimdFloat  v0,
                          SimdFloat  v1,
                          SimdFloat  v2,
                          SimdFloat  v3)
{
    _MM_TRANSPOSE4_PS(v0.r, v1.r, v2.r, v3.r);
    v0.r = _mm_add_ps(v0.r, v1.r);
    v2.r = _mm_add_ps(v2.r, v3.r);
    v0.r = _mm_add_ps(v0.r, v2.r);
    v2.r = _mm_add_ps(v0.r, _mm_load_ps(m));

    assert(std::size_t(m) % 16 == 0);
    _mm_store_ps(m, v2.r);

    __m128 b = _mm_add_ps(v0.r, _mm_shuffle_ps(v0.r, v0.r, _MM_SHUFFLE(1, 0, 3, 2)));
    b = _mm_add_ss(b, _mm_shuffle_ps(b, b, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float *>(&b);
}
#endif

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE2_UTIL_FLOAT_H
