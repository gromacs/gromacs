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

#ifndef GMX_SIMD_IMPL_X86_AVX2_256_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX2_256_SIMD_DOUBLE_H

#include "config.h"

#include <immintrin.h>

#include "gromacs/math/utilities.h"
#include "gromacs/simd/impl_x86_avx_256/impl_x86_avx_256_simd_double.h"

namespace gmx
{

static inline SimdDouble gmx_simdcall
fma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm256_fmadd_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm256_fmsub_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fnma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm256_fnmadd_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fnms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm256_fnmsub_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdDBool gmx_simdcall
testBits(SimdDouble a)
{
    __m256i ia  = _mm256_castpd_si256(a.simdInternal_);
    __m256i res = _mm256_andnot_si256( _mm256_cmpeq_epi64(ia, _mm256_setzero_si256()), _mm256_cmpeq_epi64(ia, ia));

    return {
               _mm256_castsi256_pd(res)
    };
}

static inline SimdDouble
frexp(SimdDouble value, SimdDInt32 * exponent)
{
    const __m256d  exponentMask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FF0000000000000LL));
    const __m256d  mantissaMask = _mm256_castsi256_pd( _mm256_set1_epi64x(0x800FFFFFFFFFFFFFLL));
    const __m256i  exponentBias = _mm256_set1_epi64x(1022LL); // add 1 to make our definition identical to frexp()
    const __m256d  half         = _mm256_set1_pd(0.5);
    __m256i        iExponent;
    __m128i        iExponent128;

    iExponent = _mm256_castpd_si256(_mm256_and_pd(value.simdInternal_, exponentMask));
    iExponent = _mm256_sub_epi64(_mm256_srli_epi64(iExponent, 52), exponentBias);
    iExponent = _mm256_shuffle_epi32(iExponent, _MM_SHUFFLE(3, 1, 2, 0));

    iExponent128             = _mm256_extractf128_si256(iExponent, 1);
    exponent->simdInternal_  = _mm_unpacklo_epi64(_mm256_castsi256_si128(iExponent), iExponent128);

    return {
               _mm256_or_pd(_mm256_and_pd(value.simdInternal_, mantissaMask), half)
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

    __m256i iExponent256 = _mm256_slli_epi64(_mm256_cvtepi32_epi64(iExponent), 52);
    return {
               _mm256_mul_pd(value.simdInternal_, _mm256_castsi256_pd(iExponent256))
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX2_256_SIMD_DOUBLE_H
