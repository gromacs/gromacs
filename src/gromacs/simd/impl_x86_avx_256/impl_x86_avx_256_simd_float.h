/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "gromacs/math/utilities.h"

namespace gmx
{

class SimdFloat
{
    public:
        SimdFloat() {}

        SimdFloat(float f) : simdInternal_(_mm256_set1_ps(f)) {}

        // Internal utility constructor to simplify return statements
        SimdFloat(__m256 simd) : simdInternal_(simd) {}

        __m256  simdInternal_;
};

class SimdFInt32
{
    public:
        SimdFInt32() {}

        SimdFInt32(std::int32_t i) : simdInternal_(_mm256_set1_epi32(i)) {}

        // Internal utility constructor to simplify return statements
        SimdFInt32(__m256i simd) : simdInternal_(simd) {}

        __m256i  simdInternal_;
};

class SimdFBool
{
    public:
        SimdFBool() {}

        SimdFBool(bool b) : simdInternal_(_mm256_castsi256_ps(_mm256_set1_epi32( b ? 0xFFFFFFFF : 0))) {}

        // Internal utility constructor to simplify return statements
        SimdFBool(__m256 simd) : simdInternal_(simd) {}

        __m256  simdInternal_;
};

static inline SimdFloat gmx_simdcall
simdLoad(const float *m)
{
    assert(std::size_t(m) % 32 == 0);
    return {
               _mm256_load_ps(m)
    };
}

static inline void gmx_simdcall
store(float *m, SimdFloat a)
{
    assert(std::size_t(m) % 32 == 0);
    _mm256_store_ps(m, a.simdInternal_);
}

static inline SimdFloat gmx_simdcall
simdLoadU(const float *m)
{
    return {
               _mm256_loadu_ps(m)
    };
}

static inline void gmx_simdcall
storeU(float *m, SimdFloat a)
{
    _mm256_storeu_ps(m, a.simdInternal_);
}

static inline SimdFloat gmx_simdcall
setZeroF()
{
    return {
               _mm256_setzero_ps()
    };
}

static inline SimdFInt32 gmx_simdcall
simdLoadFI(const std::int32_t * m)
{
    assert(std::size_t(m) % 32 == 0);
    return {
               _mm256_load_si256(reinterpret_cast<const __m256i *>(m))
    };
}

static inline void gmx_simdcall
store(std::int32_t * m, SimdFInt32 a)
{
    assert(std::size_t(m) % 32 == 0);
    _mm256_store_si256(reinterpret_cast<__m256i *>(m), a.simdInternal_);
}

static inline SimdFInt32 gmx_simdcall
simdLoadUFI(const std::int32_t *m)
{
    return {
               _mm256_loadu_si256(reinterpret_cast<const __m256i *>(m))
    };
}

static inline void gmx_simdcall
storeU(std::int32_t * m, SimdFInt32 a)
{
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(m), a.simdInternal_);
}

static inline SimdFInt32 gmx_simdcall
setZeroFI()
{
    return {
               _mm256_setzero_si256()
    };
}

template<int index>
static inline std::int32_t gmx_simdcall
extract(SimdFInt32 a)
{
    return _mm_extract_epi32(_mm256_extractf128_si256(a.simdInternal_, index>>2), index & 0x3);
}

static inline SimdFloat gmx_simdcall
operator&(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_and_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
andNot(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_andnot_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator|(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_or_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator^(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_xor_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator+(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_add_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator-(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_sub_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator-(SimdFloat x)
{
    return {
               _mm256_xor_ps(x.simdInternal_, _mm256_set1_ps(GMX_FLOAT_NEGZERO))
    };
}

static inline SimdFloat gmx_simdcall
operator*(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_mul_ps(a.simdInternal_, b.simdInternal_)
    };
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
static inline SimdFloat gmx_simdcall
fma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_add_ps(_mm256_mul_ps(a.simdInternal_, b.simdInternal_), c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_sub_ps(_mm256_mul_ps(a.simdInternal_, b.simdInternal_), c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fnma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_sub_ps(c.simdInternal_, _mm256_mul_ps(a.simdInternal_, b.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
fnms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_sub_ps(_mm256_setzero_ps(), _mm256_add_ps(_mm256_mul_ps(a.simdInternal_, b.simdInternal_), c.simdInternal_))
    };
}
#endif

static inline SimdFloat gmx_simdcall
rsqrt(SimdFloat x)
{
    return {
               _mm256_rsqrt_ps(x.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
rcp(SimdFloat x)
{
    return {
               _mm256_rcp_ps(x.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
maskAdd(SimdFloat a, SimdFloat b, SimdFBool m)
{
    return {
               _mm256_add_ps(a.simdInternal_, _mm256_and_ps(b.simdInternal_, m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
maskzMul(SimdFloat a, SimdFloat b, SimdFBool m)
{
    return {
               _mm256_and_ps(_mm256_mul_ps(a.simdInternal_, b.simdInternal_), m.simdInternal_)
    };
}

static inline SimdFloat
maskzFma(SimdFloat a, SimdFloat b, SimdFloat c, SimdFBool m)
{
    return {
               _mm256_and_ps(_mm256_add_ps(_mm256_mul_ps(a.simdInternal_, b.simdInternal_), c.simdInternal_), m.simdInternal_)
    };
}

static inline SimdFloat
maskzRsqrt(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = _mm256_blendv_ps(_mm256_set1_ps(1.0f), x.simdInternal_, m.simdInternal_);
#endif
    return {
               _mm256_and_ps(_mm256_rsqrt_ps(x.simdInternal_), m.simdInternal_)
    };
}

static inline SimdFloat
maskzRcp(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = _mm256_blendv_ps(_mm256_set1_ps(1.0f), x.simdInternal_, m.simdInternal_);
#endif
    return {
               _mm256_and_ps(_mm256_rcp_ps(x.simdInternal_), m.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
abs(SimdFloat x)
{
    return {
               _mm256_andnot_ps( _mm256_set1_ps(GMX_FLOAT_NEGZERO), x.simdInternal_ )
    };
}

static inline SimdFloat gmx_simdcall
max(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_max_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
min(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_min_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
round(SimdFloat x)
{
    return {
               _mm256_round_ps(x.simdInternal_, _MM_FROUND_NINT)
    };
}

static inline SimdFloat gmx_simdcall
trunc(SimdFloat x)
{
    return {
               _mm256_round_ps(x.simdInternal_, _MM_FROUND_TRUNC)
    };
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
static inline SimdFloat gmx_simdcall
frexp(SimdFloat value, SimdFInt32 * exponent)
{
    const __m256  exponentMask      = _mm256_castsi256_ps(_mm256_set1_epi32(0x7F800000));
    const __m256  mantissaMask      = _mm256_castsi256_ps(_mm256_set1_epi32(0x807FFFFF));
    const __m256  half              = _mm256_set1_ps(0.5);
    const __m128i exponentBias      = _mm_set1_epi32(126);  // add 1 to make our definition identical to frexp()
    __m256i       iExponent;
    __m128i       iExponentLow, iExponentHigh;

    iExponent               = _mm256_castps_si256(_mm256_and_ps(value.simdInternal_, exponentMask));
    iExponentHigh           = _mm256_extractf128_si256(iExponent, 0x1);
    iExponentLow            = _mm256_castsi256_si128(iExponent);
    iExponentLow            = _mm_srli_epi32(iExponentLow, 23);
    iExponentHigh           = _mm_srli_epi32(iExponentHigh, 23);
    iExponentLow            = _mm_sub_epi32(iExponentLow, exponentBias);
    iExponentHigh           = _mm_sub_epi32(iExponentHigh, exponentBias);
    iExponent               = _mm256_castsi128_si256(iExponentLow);
    exponent->simdInternal_ = _mm256_insertf128_si256(iExponent, iExponentHigh, 0x1);

    return {
               _mm256_or_ps(_mm256_and_ps(value.simdInternal_, mantissaMask), half)
    };

}

template <MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall
ldexp(SimdFloat value, SimdFInt32 exponent)
{
    const __m128i exponentBias      = _mm_set1_epi32(127);
    __m256i       iExponent;
    __m128i       iExponentLow, iExponentHigh;

    iExponentHigh = _mm256_extractf128_si256(exponent.simdInternal_, 0x1);
    iExponentLow  = _mm256_castsi256_si128(exponent.simdInternal_);

    iExponentLow  = _mm_add_epi32(iExponentLow, exponentBias);
    iExponentHigh = _mm_add_epi32(iExponentHigh, exponentBias);

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponentLow  = _mm_max_epi32(iExponentLow, _mm_setzero_si128());
        iExponentHigh = _mm_max_epi32(iExponentHigh, _mm_setzero_si128());
    }

    iExponentLow  = _mm_slli_epi32(iExponentLow, 23);
    iExponentHigh = _mm_slli_epi32(iExponentHigh, 23);
    iExponent     = _mm256_castsi128_si256(iExponentLow);
    iExponent     = _mm256_insertf128_si256(iExponent, iExponentHigh, 0x1);
    return {
               _mm256_mul_ps(value.simdInternal_, _mm256_castsi256_ps(iExponent))
    };
}
#endif

static inline float gmx_simdcall
reduce(SimdFloat a)
{
    __m128 t0;
    t0 = _mm_add_ps(_mm256_castps256_ps128(a.simdInternal_), _mm256_extractf128_ps(a.simdInternal_, 0x1));
    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float *>(&t0);
}

static inline SimdFBool gmx_simdcall
operator==(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_cmp_ps(a.simdInternal_, b.simdInternal_, _CMP_EQ_OQ)
    };
}

static inline SimdFBool gmx_simdcall
operator!=(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_cmp_ps(a.simdInternal_, b.simdInternal_, _CMP_NEQ_OQ)
    };
}

static inline SimdFBool gmx_simdcall
operator<(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_cmp_ps(a.simdInternal_, b.simdInternal_, _CMP_LT_OQ)
    };
}

static inline SimdFBool gmx_simdcall
operator<=(SimdFloat a, SimdFloat b)
{
    return {
               _mm256_cmp_ps(a.simdInternal_, b.simdInternal_, _CMP_LE_OQ)
    };
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
static inline SimdFBool gmx_simdcall
testBits(SimdFloat a)
{
    __m256 tst = _mm256_cvtepi32_ps(_mm256_castps_si256(a.simdInternal_));

    return {
               _mm256_cmp_ps(tst, _mm256_setzero_ps(), _CMP_NEQ_OQ)
    };
}
#endif

static inline SimdFBool gmx_simdcall
operator&&(SimdFBool a, SimdFBool b)
{
    return {
               _mm256_and_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
operator||(SimdFBool a, SimdFBool b)
{
    return {
               _mm256_or_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdFBool a) { return _mm256_movemask_ps(a.simdInternal_) != 0; }

static inline SimdFloat gmx_simdcall
selectByMask(SimdFloat a, SimdFBool mask)
{
    return {
               _mm256_and_ps(a.simdInternal_, mask.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
selectByNotMask(SimdFloat a, SimdFBool mask)
{
    return {
               _mm256_andnot_ps(mask.simdInternal_, a.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
blend(SimdFloat a, SimdFloat b, SimdFBool sel)
{
    return {
               _mm256_blendv_ps(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
cvtR2I(SimdFloat a)
{
    return {
               _mm256_cvtps_epi32(a.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
cvttR2I(SimdFloat a)
{
    return {
               _mm256_cvttps_epi32(a.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
cvtI2R(SimdFInt32 a)
{
    return {
               _mm256_cvtepi32_ps(a.simdInternal_)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_256_SIMD_FLOAT_H
