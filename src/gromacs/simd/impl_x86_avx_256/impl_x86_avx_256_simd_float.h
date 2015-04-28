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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_256_SIMD_FLOAT_H

#include "config.h"

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <immintrin.h>

namespace gmx
{

struct SimdFloat
{
    __m256 r;
};

struct SimdFInt32
{
    __m256i i;
};

struct SimdFBool
{
    __m256 b;
};

static inline SimdFloat gmx_simdcall
simdLoadF(const float *m)
{
    assert(std::size_t(m) % 32 == 0);
    return {
               _mm256_load_ps(m)
    };
}

static inline SimdFloat gmx_simdcall
simdLoad1F(const float *m)
{
    return {
               _mm256_broadcast_ss(m)
    };
}

static inline SimdFloat gmx_simdcall
simdSet1F(float r)
{
    return {
               _mm256_set1_ps(r)
    };
}

static inline SimdFloat gmx_simdcall
simdSetZeroF()
{
    return {
               _mm256_setzero_ps()
    };
}

static inline void gmx_simdcall
simdStoreF(float *m, SimdFloat a)
{
    assert(std::size_t(m) % 32 == 0);
    _mm256_store_ps(m, a.r);
}

static inline SimdFloat gmx_simdcall
simdLoadUF(const float *m)
{
    return {
               _mm256_loadu_ps(m)
    };
}

static inline void gmx_simdcall
simdStoreUF(float *m, SimdFloat a) { _mm256_storeu_ps(m, a.r); }

static inline SimdFInt32 gmx_simdcall
simdLoadFI(const std::int32_t * m)
{
    assert(std::size_t(m) % 32 == 0);
    return {
               _mm256_load_si256(reinterpret_cast<const __m256i *>(m))
    };
}

static inline SimdFInt32 gmx_simdcall
simdSet1FI(std::int32_t b)
{
    return {
               _mm256_set1_epi32(b)
    };
}

static inline SimdFInt32 gmx_simdcall
simdSetZeroFI()
{
    return {
               _mm256_setzero_si256()
    };
}

static inline void gmx_simdcall
simdStoreFI(std::int32_t * m, SimdFInt32 a)
{
    assert(std::size_t(m) % 32 == 0);
    _mm256_store_si256(reinterpret_cast<__m256i *>(m), a.i);
}

static inline SimdFInt32 gmx_simdcall
simdLoadUFI(const std::int32_t *m)
{
    return {
               _mm256_loadu_si256(reinterpret_cast<const __m256i *>(m))
    };
}

static inline void gmx_simdcall
simdStoreUFI(std::int32_t * m, SimdFInt32 a)
{
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(m), a.i);
}

template<int index> gmx_simdcall
static inline std::int32_t
simdExtractFI(SimdFInt32 a)
{
    return _mm_extract_epi32(_mm256_extractf128_si256(a.i, index>>2), index & 0x3);
}

static inline SimdFloat gmx_simdcall
simdAndF(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_and_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdAndNotF(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_andnot_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdOrF(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_or_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdXorF(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_xor_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdAddF(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_add_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdSubF(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_sub_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdMulF(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_mul_ps(a.r, b.r)
    };
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
static inline SimdFloat gmx_simdcall
simdFmaddF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_add_ps(_mm256_mul_ps(a.r, b.r), c.r)
    };
}

static inline SimdFloat gmx_simdcall
simdFmsubF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_sub_ps(_mm256_mul_ps(a.r, b.r), c.r)
    };
}

static inline SimdFloat gmx_simdcall
simdFnmaddF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_sub_ps(c.r, _mm256_mul_ps(a.r, b.r))
    };
}

static inline SimdFloat gmx_simdcall
simdFnmsubF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_sub_ps(_mm256_setzero_ps(), _mm256_add_ps(_mm256_mul_ps(a.r, b.r), c.r))
    };
}
#endif

static inline SimdFloat gmx_simdcall
simdRsqrtF(SimdFloat x)
{
    return {
               _mm256_rsqrt_ps(x.r)
    };
}

static inline SimdFloat gmx_simdcall
simdRcpF(SimdFloat x)
{
    return {
               _mm256_rcp_ps(x.r)
    };
}

static inline SimdFloat gmx_simdcall
simdMulMaskF(SimdFloat a, SimdFloat b, SimdFBool m)
{
    return {
               _mm256_and_ps(_mm256_mul_ps(a.r, b.r), m.b)
    };
}

static inline SimdFloat
simdFmaddMaskF(SimdFloat a, SimdFloat b, SimdFloat c, SimdFBool m)
{
    return {
               _mm256_and_ps(_mm256_add_ps(_mm256_mul_ps(a.r, b.r), c.r), m.b)
    };
}

static inline SimdFloat
simdRsqrtMaskF(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.r = _mm256_blendv_ps(_mm256_set1_ps(1.0f), x.r, m.b);
#endif
    return {
               _mm256_and_ps(_mm256_rsqrt_ps(x.r), m.b)
    };
}

static inline SimdFloat
simdRcpMaskF(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.r = _mm256_blendv_ps(_mm256_set1_ps(1.0f), x.r, m.b);
#endif
    return {
               _mm256_and_ps(_mm256_rcp_ps(x.r), m.b)
    };
}

static inline SimdFloat gmx_simdcall
simdAbsF(SimdFloat x)
{
    return {
               _mm256_andnot_ps( _mm256_set1_ps(GMX_FLOAT_NEGZERO), x.r )
    };
}

static inline SimdFloat gmx_simdcall
simdNegF(SimdFloat x)
{
    return {
               _mm256_xor_ps(x.r, _mm256_set1_ps(GMX_FLOAT_NEGZERO))
    };
}

static inline SimdFloat gmx_simdcall
simdMaxF(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_max_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdMinF(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_min_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdRoundF(SimdFloat x)
{
    return {
               _mm256_round_ps(x.r, _MM_FROUND_NINT)
    };
}

static inline SimdFloat gmx_simdcall
simdTruncF(SimdFloat x)
{
    return {
               _mm256_round_ps(x.r, _MM_FROUND_TRUNC)
    };
}

static inline SimdFloat gmx_simdcall
simdFractionF(SimdFloat x)
{
    return {
               _mm256_sub_ps(x.r, _mm256_round_ps(x.r, _MM_FROUND_TRUNC))
    };
}

static inline float gmx_simdcall
simdReduceF(SimdFloat a)
{
    __m128 t0;
    t0 = _mm_add_ps(_mm256_castps256_ps128(a.r), _mm256_extractf128_ps(a.r, 0x1));
    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float *>(&t0);
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
static inline SimdFloat gmx_simdcall
simdGetExponentF(SimdFloat x)
{
    const __m256  exponentMask      = _mm256_castsi256_ps(_mm256_set1_epi32(0x7F800000));
    const __m128i exponentBias      = _mm_set1_epi32(127);
    __m256i       iExponent;
    __m128i       iExponentLow, iExponentHigh;

    iExponent     = _mm256_castps_si256(_mm256_and_ps(x.r, exponentMask));
    iExponentHigh = _mm256_extractf128_si256(iExponent, 0x1);
    iExponentLow  = _mm256_castsi256_si128(iExponent);
    iExponentLow  = _mm_srli_epi32(iExponentLow, 23);
    iExponentHigh = _mm_srli_epi32(iExponentHigh, 23);
    iExponentLow  = _mm_sub_epi32(iExponentLow, exponentBias);
    iExponentHigh = _mm_sub_epi32(iExponentHigh, exponentBias);
    iExponent     = _mm256_castsi128_si256(iExponentLow);
    iExponent     = _mm256_insertf128_si256(iExponent, iExponentHigh, 0x1);
    return {
               _mm256_cvtepi32_ps(iExponent)
    };
}

static inline SimdFloat gmx_simdcall
simdSetExponentF(SimdFloat x)
{
    const __m128i exponentBias      = _mm_set1_epi32(127);
    __m256i       iExponent;
    __m128i       iExponentLow, iExponentHigh;

    iExponent     = _mm256_cvtps_epi32(x.r);
    iExponentHigh = _mm256_extractf128_si256(iExponent, 0x1);
    iExponentLow  = _mm256_castsi256_si128(iExponent);
    iExponentLow  = _mm_slli_epi32(_mm_add_epi32(iExponentLow, exponentBias), 23);
    iExponentHigh = _mm_slli_epi32(_mm_add_epi32(iExponentHigh, exponentBias), 23);
    iExponent     = _mm256_castsi128_si256(iExponentLow);
    iExponent     = _mm256_insertf128_si256(iExponent, iExponentHigh, 0x1);
    return {
               _mm256_castsi256_ps(iExponent)
    };
}
#endif

static inline SimdFloat gmx_simdcall
simdGetMantissaF(SimdFloat x)
{
    const __m256 mantissaMask = _mm256_castsi256_ps(_mm256_set1_epi32(0x007FFFFF));
    const __m256 one          = _mm256_set1_ps(1.0);

    return {
               _mm256_or_ps(_mm256_and_ps(x.r, mantissaMask), one)
    };
}

static inline SimdFBool gmx_simdcall
simdCmpEqF(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_cmp_ps(a.r, b.r, _CMP_EQ_OQ)
    };
}

static inline SimdFBool gmx_simdcall
simdCmpNzF(SimdFloat a)
{
    return {
               _mm256_cmp_ps(a.r, _mm256_setzero_ps(), _CMP_NEQ_OQ)
    };
}

static inline SimdFBool gmx_simdcall
simdCmpLtF(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_cmp_ps(a.r, b.r, _CMP_LT_OQ)
    };
}

static inline SimdFBool gmx_simdcall
simdCmpLeF(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_cmp_ps(a.r, b.r, _CMP_LE_OQ)
    };
}

static inline SimdFBool gmx_simdcall
simdAndFB(SimdFBool a, SimdFBool b)
{
    return {
               _mm256_and_ps(a.b, b.b)
    };
}

static inline SimdFBool gmx_simdcall
simdOrFB(SimdFBool a, SimdFBool b)
{
    return {
               _mm256_or_ps(a.b, b.b)
    };
}

static inline bool gmx_simdcall
simdAnyTrueFB(SimdFBool a) { return _mm256_movemask_ps(a.b) != 0; }

static inline SimdFloat gmx_simdcall
simdMaskF(SimdFloat a, SimdFBool mask)
{
    return {
               _mm256_and_ps(a.r, mask.b)
    };
}

static inline SimdFloat gmx_simdcall
simdMaskNotF(SimdFloat a, SimdFBool mask)
{
    return {
               _mm256_andnot_ps(mask.b, a.r)
    };
}

static inline SimdFloat gmx_simdcall
simdBlendF(SimdFloat a, SimdFloat b, SimdFBool sel)
{
    return {
               _mm256_blendv_ps(a.r, b.r, sel.b)
    };
}

static inline SimdFInt32 gmx_simdcall
simdCvtF2I(SimdFloat a)
{
    return {
               _mm256_cvtps_epi32(a.r)
    };
}

static inline SimdFInt32 gmx_simdcall
simdCvttF2I(SimdFloat a)
{
    return {
               _mm256_cvttps_epi32(a.r)
    };
}

static inline SimdFloat gmx_simdcall
simdCvtI2F(SimdFInt32 a)
{
    return {
               _mm256_cvtepi32_ps(a.i)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_256_SIMD_FLOAT_H
