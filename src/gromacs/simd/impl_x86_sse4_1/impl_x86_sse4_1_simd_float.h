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

#ifndef GMX_SIMD_IMPL_X86_SSE4_1_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_SSE4_1_SIMD_FLOAT_H

#include "config.h"

#include <smmintrin.h>

#include "gromacs/simd/impl_x86_sse2/impl_x86_sse2_simd_float.h"

namespace gmx
{

template<int index>
static inline std::int32_t gmx_simdcall
extract(SimdFInt32 a)
{
    return _mm_extract_epi32(a.simdInternal_, index);
}

static inline SimdFloat
maskzRsqrt(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = _mm_blendv_ps(_mm_set1_ps(1.0f), x.simdInternal_, m.simdInternal_);
#endif
    return {
               _mm_and_ps(_mm_rsqrt_ps(x.simdInternal_), m.simdInternal_)
    };
}

static inline SimdFloat
maskzRcp(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = _mm_blendv_ps(_mm_set1_ps(1.0f), x.simdInternal_, m.simdInternal_);
#endif
    return {
               _mm_and_ps(_mm_rcp_ps(x.simdInternal_), m.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
round(SimdFloat x)
{
    return {
               _mm_round_ps(x.simdInternal_, _MM_FROUND_NINT)
    };
}

static inline SimdFloat gmx_simdcall
trunc(SimdFloat x)
{
    return {
               _mm_round_ps(x.simdInternal_, _MM_FROUND_TRUNC)
    };
}

static inline SimdFloat gmx_simdcall
blend(SimdFloat a, SimdFloat b, SimdFBool sel)
{
    return {
               _mm_blendv_ps(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator*(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_mullo_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
blend(SimdFInt32 a, SimdFInt32 b, SimdFIBool sel)
{
    return {
               _mm_blendv_epi8(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

template <MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall
ldexp(SimdFloat value, SimdFInt32 exponent)
{
    const __m128i exponentBias = _mm_set1_epi32(127);
    __m128i       iExponent;

    iExponent = _mm_add_epi32(exponent.simdInternal_, exponentBias);

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponent = _mm_max_epi32(iExponent, _mm_setzero_si128());
    }

    iExponent = _mm_slli_epi32( iExponent, 23);

    return {
               _mm_mul_ps(value.simdInternal_, _mm_castsi128_ps(iExponent))
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE4_1_SIMD_FLOAT_H
