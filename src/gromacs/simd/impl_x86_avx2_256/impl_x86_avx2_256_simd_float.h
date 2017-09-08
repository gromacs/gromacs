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

#ifndef GMX_SIMD_IMPL_X86_AVX2_256_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX2_256_SIMD_FLOAT_H

#include "config.h"

#include <immintrin.h>

#include "gromacs/math/utilities.h"
#include "gromacs/simd/impl_x86_avx_256/impl_x86_avx_256_simd_float.h"

namespace gmx
{

class SimdFIBool
{
    public:
        SimdFIBool() {}

        SimdFIBool(bool b) : simdInternal_(_mm256_set1_epi32( b ? 0xFFFFFFFF : 0)) {}

        // Internal utility constructor to simplify return statements
        SimdFIBool(__m256i simd) : simdInternal_(simd) {}

        __m256i  simdInternal_;
};

static inline SimdFloat gmx_simdcall
fma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_fmadd_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_fmsub_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fnma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_fnmadd_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fnms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_fnmsub_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
testBits(SimdFloat a)
{
    __m256i ia  = _mm256_castps_si256(a.simdInternal_);
    __m256i res = _mm256_andnot_si256( _mm256_cmpeq_epi32(ia, _mm256_setzero_si256()), _mm256_cmpeq_epi32(ia, ia));

    return {
               _mm256_castsi256_ps(res)
    };
}

static inline SimdFloat gmx_simdcall
frexp(SimdFloat value, SimdFInt32 * exponent)
{
    const __m256  exponentMask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7F800000));
    const __m256  mantissaMask = _mm256_castsi256_ps(_mm256_set1_epi32(0x807FFFFF));
    const __m256i exponentBias = _mm256_set1_epi32(126); // add 1 to make our definition identical to frexp()
    const __m256  half         = _mm256_set1_ps(0.5);
    __m256i       iExponent;

    iExponent               = _mm256_castps_si256(_mm256_and_ps(value.simdInternal_, exponentMask));
    exponent->simdInternal_ = _mm256_sub_epi32(_mm256_srli_epi32(iExponent, 23), exponentBias);

    return {
               _mm256_or_ps(_mm256_and_ps(value.simdInternal_, mantissaMask), half)
    };
}

template <MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall
ldexp(SimdFloat value, SimdFInt32 exponent)
{
    const __m256i  exponentBias = _mm256_set1_epi32(127);
    __m256i        iExponent    = _mm256_add_epi32(exponent.simdInternal_, exponentBias);

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponent = _mm256_max_epi32(iExponent, _mm256_setzero_si256());
    }

    iExponent = _mm256_slli_epi32(iExponent, 23);
    return {
               _mm256_mul_ps(value.simdInternal_, _mm256_castsi256_ps(iExponent))
    };
}

static inline SimdFInt32 gmx_simdcall
operator<<(SimdFInt32 a, int n)
{
    return {
               _mm256_slli_epi32(a.simdInternal_, n)
    };
}

static inline SimdFInt32 gmx_simdcall
operator>>(SimdFInt32 a, int n)
{
    return {
               _mm256_srli_epi32(a.simdInternal_, n)
    };
}

static inline SimdFInt32 gmx_simdcall
operator&(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_and_si256(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
andNot(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_andnot_si256(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator|(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_or_si256(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator^(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_xor_si256(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator+(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_add_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator-(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_sub_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator*(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_mullo_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator==(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_cmpeq_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
testBits(SimdFInt32 a)
{
    return {
               _mm256_andnot_si256(_mm256_cmpeq_epi32(a.simdInternal_, _mm256_setzero_si256()),
                                   _mm256_cmpeq_epi32(a.simdInternal_, a.simdInternal_))
    };
}

static inline SimdFIBool gmx_simdcall
operator<(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_cmpgt_epi32(b.simdInternal_, a.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator&&(SimdFIBool a, SimdFIBool b)
{
    return {
               _mm256_and_si256(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator||(SimdFIBool a, SimdFIBool b)
{
    return {
               _mm256_or_si256(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdFIBool a) { return _mm256_movemask_epi8(a.simdInternal_) != 0; }

static inline SimdFInt32 gmx_simdcall
selectByMask(SimdFInt32 a, SimdFIBool mask)
{
    return {
               _mm256_and_si256(a.simdInternal_, mask.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
selectByNotMask(SimdFInt32 a, SimdFIBool mask)
{
    return {
               _mm256_andnot_si256(mask.simdInternal_, a.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
blend(SimdFInt32 a, SimdFInt32 b, SimdFIBool sel)
{
    return {
               _mm256_blendv_epi8(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
cvtB2IB(SimdFBool a)
{
    return {
               _mm256_castps_si256(a.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
cvtIB2B(SimdFIBool a)
{
    return {
               _mm256_castsi256_ps(a.simdInternal_)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX2_256_SIMD_FLOAT_H
