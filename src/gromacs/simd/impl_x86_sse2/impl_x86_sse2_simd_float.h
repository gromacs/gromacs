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
#ifndef GMX_SIMD_IMPL_X86_SSE2_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_SSE2_SIMD_FLOAT_H

#include "config.h"

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <emmintrin.h>

namespace gmx
{

struct SimdFloat
{
    __m128 r;
};

struct SimdFInt32
{
    __m128i i;
};

struct SimdFBool
{
    __m128 b;
};

struct SimdFIBool
{
    __m128i b;
};

static inline SimdFloat gmx_simdcall
simdLoadF(const float *m)
{
    assert(std::size_t(m) % 16 == 0);
    return {
               _mm_load_ps(m)
    };
}

static inline SimdFloat gmx_simdcall
simdLoad1F(const float *m)
{
    return {
               _mm_load1_ps(m)
    };
}

static inline SimdFloat gmx_simdcall
simdSet1F(float r)
{
    return {
               _mm_set1_ps(r)
    };
}

static inline SimdFloat gmx_simdcall
simdSetZeroF()
{
    return {
               _mm_setzero_ps()
    };
}

static inline void gmx_simdcall
simdStoreF(float *m, SimdFloat a)
{
    assert(std::size_t(m) % 16 == 0);
    _mm_store_ps(m, a.r);
}

static inline SimdFloat gmx_simdcall
simdLoadUF(const float *m)
{
    return {
               _mm_loadu_ps(m)
    };
}

static inline void gmx_simdcall
simdStoreUF(float *m, SimdFloat a) { _mm_storeu_ps(m, a.r); }

static inline SimdFInt32 gmx_simdcall
simdLoadFI(const std::int32_t * m)
{
    assert(std::size_t(m) % 16 == 0);
    return {
               _mm_load_si128(reinterpret_cast<const __m128i *>(m))
    };
}

static inline SimdFInt32 gmx_simdcall
simdSet1FI(std::int32_t b)
{
    return {
               _mm_set1_epi32(b)
    };
}

static inline SimdFInt32 gmx_simdcall
simdSetZeroFI()
{
    return {
               _mm_setzero_si128()
    };
}

static inline void gmx_simdcall
simdStoreFI(std::int32_t * m, SimdFInt32 a)
{
    assert(std::size_t(m) % 16 == 0);
    _mm_store_si128(reinterpret_cast<__m128i *>(m), a.i);
}

static inline SimdFInt32 gmx_simdcall
simdLoadUFI(const std::int32_t *m)
{
    return {
               _mm_loadu_si128(reinterpret_cast<const __m128i *>(m))
    };
}


static inline void gmx_simdcall
simdStoreUFI(std::int32_t * m, SimdFInt32 a)
{
    _mm_storeu_si128(reinterpret_cast<__m128i *>(m), a.i);
}


// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
template<int index> gmx_simdcall
static inline std::int32_t
simdExtractFI(SimdFInt32 a)
{
    return _mm_cvtsi128_si32( _mm_srli_si128(a.i, 4 * index) );
}
#endif

static inline SimdFloat gmx_simdcall
simdAndF(SimdFloat a, SimdFloat b)
{
    return {
               _mm_and_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdAndNotF(SimdFloat a, SimdFloat b)
{
    return {
               _mm_andnot_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdOrF(SimdFloat a, SimdFloat b)
{
    return {
               _mm_or_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdXorF(SimdFloat a, SimdFloat b)
{
    return {
               _mm_xor_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdAddF(SimdFloat a, SimdFloat b)
{
    return {
               _mm_add_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdSubF(SimdFloat a, SimdFloat b)
{
    return {
               _mm_sub_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdMulF(SimdFloat a, SimdFloat b)
{
    return {
               _mm_mul_ps(a.r, b.r)
    };
}

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline SimdFloat gmx_simdcall
simdFmaddF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm_add_ps(_mm_mul_ps(a.r, b.r), c.r)
    };
}

static inline SimdFloat gmx_simdcall
simdFmsubF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm_sub_ps(_mm_mul_ps(a.r, b.r), c.r)
    };
}

static inline SimdFloat gmx_simdcall
simdFnmaddF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm_sub_ps(c.r, _mm_mul_ps(a.r, b.r))
    };
}

static inline SimdFloat gmx_simdcall
simdFnmsubF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm_sub_ps(_mm_setzero_ps(), _mm_add_ps(_mm_mul_ps(a.r, b.r), c.r))
    };
}
#endif

static inline SimdFloat gmx_simdcall
simdRsqrtF(SimdFloat x)
{
    return {
               _mm_rsqrt_ps(x.r)
    };
}

static inline SimdFloat gmx_simdcall
simdRcpF(SimdFloat x)
{
    return {
               _mm_rcp_ps(x.r)
    };
}

static inline SimdFloat gmx_simdcall
simdMulMaskF(SimdFloat a, SimdFloat b, SimdFBool m)
{
    return {
               _mm_and_ps(_mm_mul_ps(a.r, b.r), m.b)
    };
}

static inline SimdFloat gmx_simdcall
simdFmaddMaskF(SimdFloat a, SimdFloat b, SimdFloat c, SimdFBool m)
{
    return {
               _mm_and_ps(_mm_add_ps(_mm_mul_ps(a.r, b.r), c.r), m.b)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdFloat gmx_simdcall
simdRsqrtMaskF(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.r = _mm_or_ps(_mm_andnot_ps(m.b, _mm_set1_ps(1.0f)), _mm_and_ps(m.b, x.r));
#endif
    return {
               _mm_and_ps(_mm_rsqrt_ps(x.r), m.b)
    };
}

static inline SimdFloat gmx_simdcall
simdRcpMaskF(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.r = _mm_or_ps(_mm_andnot_ps(m.b, _mm_set1_ps(1.0f)), _mm_and_ps(m.b, x.r));
#endif
    return {
               _mm_and_ps(_mm_rcp_ps(x.r), m.b)
    };
}
#endif

static inline SimdFloat gmx_simdcall
simdAbsF(SimdFloat x)
{
    return {
               _mm_andnot_ps( _mm_set1_ps(GMX_FLOAT_NEGZERO), x.r )
    };
}

static inline SimdFloat gmx_simdcall
simdNegF(SimdFloat x)
{
    return {
               _mm_xor_ps(x.r, _mm_set1_ps(GMX_FLOAT_NEGZERO))
    };
}

static inline SimdFloat gmx_simdcall
simdMaxF(SimdFloat a, SimdFloat b)
{
    return {
               _mm_max_ps(a.r, b.r)
    };
}

static inline SimdFloat gmx_simdcall
simdMinF(SimdFloat a, SimdFloat b)
{
    return {
               _mm_min_ps(a.r, b.r)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdFloat gmx_simdcall
simdRoundF(SimdFloat x)
{
    return {
               _mm_cvtepi32_ps( _mm_cvtps_epi32(x.r) )
    };
}

static inline SimdFloat gmx_simdcall
simdTruncF(SimdFloat x)
{
    return {
               _mm_cvtepi32_ps( _mm_cvttps_epi32(x.r) )
    };
}

static inline SimdFloat gmx_simdcall
simdFractionF(SimdFloat x)
{
    return {
               _mm_sub_ps(x.r, _mm_cvtepi32_ps( _mm_cvttps_epi32(x.r) ))
    };
}
#endif

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline float gmx_simdcall
simdReduceF(SimdFloat a)
{
    // Shuffle has latency 1/throughput 1, followed by add with latency 3, t-put 1.
    // This is likely faster than using _mm_hadd_ps, which has latency 5, t-put 2.
    a.r = _mm_add_ps(a.r, _mm_shuffle_ps(a.r, a.r, _MM_SHUFFLE(1, 0, 3, 2)));
    a.r = _mm_add_ss(a.r, _mm_shuffle_ps(a.r, a.r, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float *>(&a);
}
#endif

static inline SimdFloat gmx_simdcall
simdGetExponentF(SimdFloat x)
{
    const __m128  exponentMask      = _mm_castsi128_ps(_mm_set1_epi32(0x7F800000));
    const __m128i exponentBias      = _mm_set1_epi32(127);
    __m128i       iExponent;

    iExponent = _mm_castps_si128(_mm_and_ps(x.r, exponentMask));
    iExponent = _mm_sub_epi32(_mm_srli_epi32(iExponent, 23), exponentBias);

    return {
               _mm_cvtepi32_ps(iExponent)
    };
}

static inline SimdFloat gmx_simdcall
simdGetMantissaF(SimdFloat x)
{
    const __m128 mantissaMask = _mm_castsi128_ps(_mm_set1_epi32(0x007FFFFF));
    const __m128 one          = _mm_set1_ps(1.0f);

    return {
               _mm_or_ps( _mm_and_ps(x.r, mantissaMask), one)
    };
}

static inline SimdFloat gmx_simdcall
simdSetExponentF(SimdFloat x)
{
    const __m128i exponentBias = _mm_set1_epi32(127);
    __m128i       iExponent    = _mm_cvtps_epi32(x.r);

    iExponent = _mm_slli_epi32( _mm_add_epi32(iExponent, exponentBias), 23);

    return {
               _mm_castsi128_ps(iExponent)
    };
}

static inline SimdFBool gmx_simdcall
simdCmpEqF(SimdFloat a, SimdFloat b)
{
    return {
               _mm_cmpeq_ps(a.r, b.r)
    };
}

static inline SimdFBool gmx_simdcall
simdCmpNzF(SimdFloat a)
{
    return {
               _mm_cmpneq_ps(a.r, _mm_setzero_ps())
    };
}

static inline SimdFBool gmx_simdcall
simdCmpLtF(SimdFloat a, SimdFloat b)
{
    return {
               _mm_cmplt_ps(a.r, b.r)
    };
}

static inline SimdFBool gmx_simdcall
simdCmpLeF(SimdFloat a, SimdFloat b)
{
    return {
               _mm_cmple_ps(a.r, b.r)
    };
}

static inline SimdFBool gmx_simdcall
simdAndFB(SimdFBool a, SimdFBool b)
{
    return {
               _mm_and_ps(a.b, b.b)
    };
}

static inline SimdFBool gmx_simdcall
simdOrFB(SimdFBool a, SimdFBool b)
{
    return {
               _mm_or_ps(a.b, b.b)
    };
}

static inline bool gmx_simdcall
simdAnyTrueFB(SimdFBool a) { return _mm_movemask_ps(a.b) != 0; }

static inline SimdFloat gmx_simdcall
simdMaskF(SimdFloat a, SimdFBool mask)
{
    return {
               _mm_and_ps(a.r, mask.b)
    };
}

static inline SimdFloat gmx_simdcall
simdMaskNotF(SimdFloat a, SimdFBool mask)
{
    return {
               _mm_andnot_ps(mask.b, a.r)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdFloat gmx_simdcall
simdBlendF(SimdFloat a, SimdFloat b, SimdFBool sel)
{
    return {
               _mm_or_ps(_mm_andnot_ps(sel.b, a.r), _mm_and_ps(sel.b, b.r))
    };
}
#endif

static inline SimdFInt32 gmx_simdcall
simdSlliFI(SimdFInt32 a, int n)
{
    return {
               _mm_slli_epi32(a.i, n)
    };
}

static inline SimdFInt32 gmx_simdcall
simdSrliFI(SimdFInt32 a, int n)
{
    return {
               _mm_srli_epi32(a.i, n)
    };
}

static inline SimdFInt32 gmx_simdcall
simdAndFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_and_si128(a.i, b.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdAndNotFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_andnot_si128(a.i, b.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdOrFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_or_si128(a.i, b.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdXorFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_xor_si128(a.i, b.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdAddFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_add_epi32(a.i, b.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdSubFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_sub_epi32(a.i, b.i)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdFInt32 gmx_simdcall
simdMulFI(SimdFInt32 a, SimdFInt32 b)
{
    __m128i a1 = _mm_srli_si128(a.i, 4); // - a[3] a[2] a[1]
    __m128i b1 = _mm_srli_si128(b.i, 4); // - b[3] b[2] b[1]
    __m128i c  = _mm_mul_epu32(a.i, b.i);
    __m128i c1 = _mm_mul_epu32(a1, b1);

    c  = _mm_shuffle_epi32(c, _MM_SHUFFLE(3, 1, 2, 0));  // - - a[2]*b[2] a[0]*b[0]
    c1 = _mm_shuffle_epi32(c1, _MM_SHUFFLE(3, 1, 2, 0)); // - - a[3]*b[3] a[1]*b[1]

    return {
               _mm_unpacklo_epi32(c, c1)
    };
}
#endif

static inline SimdFIBool gmx_simdcall
simdCmpEqFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_cmpeq_epi32(a.i, b.i)
    };
}

static inline SimdFIBool gmx_simdcall
simdCmpLtFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_cmplt_epi32(a.i, b.i)
    };
}

static inline SimdFIBool gmx_simdcall
simdAndFIB(SimdFIBool a, SimdFIBool b)
{
    return {
               _mm_and_si128(a.b, b.b)
    };
}

static inline SimdFIBool gmx_simdcall
simdOrFIB(SimdFIBool a, SimdFIBool b)
{
    return {
               _mm_or_si128(a.b, b.b)
    };
}

static inline bool gmx_simdcall
simdAnyTrueFIB(SimdFIBool a) { return _mm_movemask_epi8(a.b) != 0; }

static inline SimdFInt32 gmx_simdcall
simdMaskFI(SimdFInt32 a, SimdFIBool mask)
{
    return {
               _mm_and_si128(a.i, mask.b)
    };
}

static inline SimdFInt32 gmx_simdcall
simdMaskNotFI(SimdFInt32 a, SimdFIBool mask)
{
    return {
               _mm_andnot_si128(mask.b, a.i)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdFInt32 gmx_simdcall
simdBlendFI(SimdFInt32 a, SimdFInt32 b, SimdFIBool sel)
{
    return {
               _mm_or_si128(_mm_andnot_si128(sel.b, a.i), _mm_and_si128(sel.b, b.i))
    };
}
#endif

static inline SimdFInt32 gmx_simdcall
simdCvtF2I(SimdFloat a)
{
    return {
               _mm_cvtps_epi32(a.r)
    };
}

static inline SimdFInt32 gmx_simdcall
simdCvttF2I(SimdFloat a)
{
    return {
               _mm_cvttps_epi32(a.r)
    };
}

static inline SimdFloat gmx_simdcall
simdCvtI2F(SimdFInt32 a)
{
    return {
               _mm_cvtepi32_ps(a.i)
    };
}

static inline SimdFIBool gmx_simdcall
simdCvtFB2FIB(SimdFBool a)
{
    return {
               _mm_castps_si128(a.b)
    };
}

static inline SimdFBool gmx_simdcall
simdCvtFIB2FB(SimdFIBool a)
{
    return {
               _mm_castsi128_ps(a.b)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE2_SIMD_FLOAT_H
