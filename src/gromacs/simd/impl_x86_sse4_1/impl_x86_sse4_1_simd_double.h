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

#ifndef GMX_SIMD_IMPL_X86_SSE4_1_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_X86_SSE4_1_SIMD_DOUBLE_H

#include "config.h"

#include <smmintrin.h>

#include "gromacs/simd/impl_x86_sse2/impl_x86_sse2_simd_double.h"

namespace gmx
{

template<int index>
static inline std::int32_t gmx_simdcall
extract(SimdDInt32 a)
{
    return _mm_extract_epi32(a.simdInternal_, index);
}

static inline SimdDouble
maskzRsqrt(SimdDouble x, SimdDBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = _mm_blendv_pd(_mm_set1_pd(1.0), x.simdInternal_, m.simdInternal_);
#endif
    return {
               _mm_and_pd(_mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(x.simdInternal_))), m.simdInternal_)
    };
}

static inline SimdDouble
maskzRcp(SimdDouble x, SimdDBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = _mm_blendv_pd(_mm_set1_pd(1.0), x.simdInternal_, m.simdInternal_);
#endif
    return {
               _mm_and_pd(_mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(x.simdInternal_))), m.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
round(SimdDouble x)
{
    return {
               _mm_round_pd(x.simdInternal_, _MM_FROUND_NINT)
    };
}

static inline SimdDouble gmx_simdcall
trunc(SimdDouble x)
{
    return {
               _mm_round_pd(x.simdInternal_, _MM_FROUND_TRUNC)
    };
}

static inline SimdDBool gmx_simdcall
testBits(SimdDouble a)
{
    __m128i ia  = _mm_castpd_si128(a.simdInternal_);
    __m128i res = _mm_andnot_si128( _mm_cmpeq_epi64(ia, _mm_setzero_si128()), _mm_cmpeq_epi64(ia, ia));

    return {
               _mm_castsi128_pd(res)
    };
}

static inline SimdDouble gmx_simdcall
blend(SimdDouble a, SimdDouble b, SimdDBool sel)
{
    return {
               _mm_blendv_pd(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator*(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_mullo_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
blend(SimdDInt32 a, SimdDInt32 b, SimdDIBool sel)
{
    return {
               _mm_blendv_epi8(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

template <MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble
ldexp(SimdDouble value, SimdDInt32 exponent)
{
    const __m128i  exponentBias = _mm_set1_epi32(1023);
    __m128i        iExponent    = _mm_add_epi32(exponent.simdInternal_, exponentBias);

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponent = _mm_max_epi32(iExponent, _mm_setzero_si128());
    }

    // After conversion integers will be in slot 0,1. Move them to 0,2 so
    // we can do a 64-bit shift and get them to the dp exponents.
    iExponent = _mm_shuffle_epi32(iExponent, _MM_SHUFFLE(3, 1, 2, 0));
    iExponent = _mm_slli_epi64(iExponent, 52);

    return {
               _mm_mul_pd(value.simdInternal_, _mm_castsi128_pd(iExponent))
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE4_1_SIMD_DOUBLE_H
