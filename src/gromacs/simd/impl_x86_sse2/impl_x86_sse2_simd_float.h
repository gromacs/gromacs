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
#ifndef GMX_SIMD_IMPL_X86_SSE2_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_SSE2_SIMD_FLOAT_H

#include "config.h"

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <emmintrin.h>

#include "gromacs/math/utilities.h"

namespace gmx
{

class SimdFloat
{
    public:
        SimdFloat() {}

        SimdFloat(float f) : simdInternal_(_mm_set1_ps(f)) {}

        // Internal utility constructor to simplify return statements
        SimdFloat(__m128 simd) : simdInternal_(simd) {}

        __m128  simdInternal_;
};

class SimdFInt32
{
    public:
        SimdFInt32() {}

        SimdFInt32(std::int32_t i) : simdInternal_(_mm_set1_epi32(i)) {}

        // Internal utility constructor to simplify return statements
        SimdFInt32(__m128i simd) : simdInternal_(simd) {}

        __m128i  simdInternal_;
};

class SimdFBool
{
    public:
        SimdFBool() {}

        SimdFBool(bool b) : simdInternal_(_mm_castsi128_ps(_mm_set1_epi32( b ? 0xFFFFFFFF : 0))) {}

        // Internal utility constructor to simplify return statements
        SimdFBool(__m128 simd) : simdInternal_(simd) {}

        __m128  simdInternal_;
};

class SimdFIBool
{
    public:
        SimdFIBool() {}

        SimdFIBool(bool b) : simdInternal_(_mm_set1_epi32( b ? 0xFFFFFFFF : 0)) {}

        // Internal utility constructor to simplify return statements
        SimdFIBool(__m128i simd) : simdInternal_(simd) {}

        __m128i  simdInternal_;
};

static inline SimdFloat gmx_simdcall
simdLoad(const float *m)
{
    assert(std::size_t(m) % 16 == 0);
    return {
               _mm_load_ps(m)
    };
}

static inline void gmx_simdcall
store(float *m, SimdFloat a)
{
    assert(std::size_t(m) % 16 == 0);
    _mm_store_ps(m, a.simdInternal_);
}

static inline SimdFloat gmx_simdcall
simdLoadU(const float *m)
{
    return {
               _mm_loadu_ps(m)
    };
}

static inline void gmx_simdcall
storeU(float *m, SimdFloat a) { _mm_storeu_ps(m, a.simdInternal_); }

static inline SimdFloat gmx_simdcall
setZeroF()
{
    return {
               _mm_setzero_ps()
    };
}

static inline SimdFInt32 gmx_simdcall
simdLoadFI(const std::int32_t * m)
{
    assert(std::size_t(m) % 16 == 0);
    return {
               _mm_load_si128(reinterpret_cast<const __m128i *>(m))
    };
}

static inline void gmx_simdcall
store(std::int32_t * m, SimdFInt32 a)
{
    assert(std::size_t(m) % 16 == 0);
    _mm_store_si128(reinterpret_cast<__m128i *>(m), a.simdInternal_);
}

static inline SimdFInt32 gmx_simdcall
simdLoadUFI(const std::int32_t *m)
{
    return {
               _mm_loadu_si128(reinterpret_cast<const __m128i *>(m))
    };
}

static inline void gmx_simdcall
storeU(std::int32_t * m, SimdFInt32 a)
{
    _mm_storeu_si128(reinterpret_cast<__m128i *>(m), a.simdInternal_);
}

static inline SimdFInt32 gmx_simdcall
setZeroFI()
{
    return {
               _mm_setzero_si128()
    };
}


// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
template<int index>
static inline std::int32_t gmx_simdcall
extract(SimdFInt32 a)
{
    return _mm_cvtsi128_si32( _mm_srli_si128(a.simdInternal_, 4 * index) );
}
#endif

static inline SimdFloat gmx_simdcall
operator&(SimdFloat a, SimdFloat b)
{
    return {
               _mm_and_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
andNot(SimdFloat a, SimdFloat b)
{
    return {
               _mm_andnot_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator|(SimdFloat a, SimdFloat b)
{
    return {
               _mm_or_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator^(SimdFloat a, SimdFloat b)
{
    return {
               _mm_xor_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator+(SimdFloat a, SimdFloat b)
{
    return {
               _mm_add_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator-(SimdFloat a, SimdFloat b)
{
    return {
               _mm_sub_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator-(SimdFloat x)
{
    return {
               _mm_xor_ps(x.simdInternal_, _mm_set1_ps(GMX_FLOAT_NEGZERO))
    };
}

static inline SimdFloat gmx_simdcall
operator*(SimdFloat a, SimdFloat b)
{
    return {
               _mm_mul_ps(a.simdInternal_, b.simdInternal_)
    };
}

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline SimdFloat gmx_simdcall
fma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm_add_ps(_mm_mul_ps(a.simdInternal_, b.simdInternal_), c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm_sub_ps(_mm_mul_ps(a.simdInternal_, b.simdInternal_), c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fnma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm_sub_ps(c.simdInternal_, _mm_mul_ps(a.simdInternal_, b.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
fnms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm_sub_ps(_mm_setzero_ps(), _mm_add_ps(_mm_mul_ps(a.simdInternal_, b.simdInternal_), c.simdInternal_))
    };
}
#endif

static inline SimdFloat gmx_simdcall
rsqrt(SimdFloat x)
{
    return {
               _mm_rsqrt_ps(x.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
rcp(SimdFloat x)
{
    return {
               _mm_rcp_ps(x.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
maskAdd(SimdFloat a, SimdFloat b, SimdFBool m)
{
    return {
               _mm_add_ps(a.simdInternal_, _mm_and_ps(b.simdInternal_, m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
maskzMul(SimdFloat a, SimdFloat b, SimdFBool m)
{
    return {
               _mm_and_ps(_mm_mul_ps(a.simdInternal_, b.simdInternal_), m.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
maskzFma(SimdFloat a, SimdFloat b, SimdFloat c, SimdFBool m)
{
    return {
               _mm_and_ps(_mm_add_ps(_mm_mul_ps(a.simdInternal_, b.simdInternal_), c.simdInternal_), m.simdInternal_)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdFloat gmx_simdcall
maskzRsqrt(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = _mm_or_ps(_mm_andnot_ps(m.simdInternal_, _mm_set1_ps(1.0f)), _mm_and_ps(m.simdInternal_, x.simdInternal_));
#endif
    return {
               _mm_and_ps(_mm_rsqrt_ps(x.simdInternal_), m.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
maskzRcp(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = _mm_or_ps(_mm_andnot_ps(m.simdInternal_, _mm_set1_ps(1.0f)), _mm_and_ps(m.simdInternal_, x.simdInternal_));
#endif
    return {
               _mm_and_ps(_mm_rcp_ps(x.simdInternal_), m.simdInternal_)
    };
}
#endif

static inline SimdFloat gmx_simdcall
abs(SimdFloat x)
{
    return {
               _mm_andnot_ps( _mm_set1_ps(GMX_FLOAT_NEGZERO), x.simdInternal_ )
    };
}

static inline SimdFloat gmx_simdcall
max(SimdFloat a, SimdFloat b)
{
    return {
               _mm_max_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
min(SimdFloat a, SimdFloat b)
{
    return {
               _mm_min_ps(a.simdInternal_, b.simdInternal_)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdFloat gmx_simdcall
round(SimdFloat x)
{
    return {
               _mm_cvtepi32_ps( _mm_cvtps_epi32(x.simdInternal_) )
    };
}

static inline SimdFloat gmx_simdcall
trunc(SimdFloat x)
{
    return {
               _mm_cvtepi32_ps( _mm_cvttps_epi32(x.simdInternal_) )
    };
}

#endif

static inline SimdFloat gmx_simdcall
frexp(SimdFloat value, SimdFInt32 * exponent)
{
    const __m128  exponentMask   = _mm_castsi128_ps(_mm_set1_epi32(0x7F800000));
    const __m128  mantissaMask   = _mm_castsi128_ps(_mm_set1_epi32(0x807FFFFF));
    const __m128i exponentBias   = _mm_set1_epi32(126); // add 1 to make our definition identical to frexp()
    const __m128  half           = _mm_set1_ps(0.5f);
    __m128i       iExponent;

    iExponent               = _mm_castps_si128(_mm_and_ps(value.simdInternal_, exponentMask));
    iExponent               = _mm_sub_epi32(_mm_srli_epi32(iExponent, 23), exponentBias);
    exponent->simdInternal_ = iExponent;

    return {
               _mm_or_ps( _mm_and_ps(value.simdInternal_, mantissaMask), half)
    };
}

// Override for SSE4.1
#if GMX_SIMD_X86_SSE2
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
        iExponent = _mm_and_si128(iExponent, _mm_cmpgt_epi32(iExponent, _mm_setzero_si128()));
    }

    iExponent = _mm_slli_epi32( iExponent, 23);

    return {
               _mm_mul_ps(value.simdInternal_, _mm_castsi128_ps(iExponent))
    };
}
#endif

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline float gmx_simdcall
reduce(SimdFloat a)
{
    // Shuffle has latency 1/throughput 1, followed by add with latency 3, t-put 1.
    // This is likely faster than using _mm_hadd_ps, which has latency 5, t-put 2.
    a.simdInternal_ = _mm_add_ps(a.simdInternal_, _mm_shuffle_ps(a.simdInternal_, a.simdInternal_, _MM_SHUFFLE(1, 0, 3, 2)));
    a.simdInternal_ = _mm_add_ss(a.simdInternal_, _mm_shuffle_ps(a.simdInternal_, a.simdInternal_, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float *>(&a);
}
#endif

static inline SimdFBool gmx_simdcall
operator==(SimdFloat a, SimdFloat b)
{
    return {
               _mm_cmpeq_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
operator!=(SimdFloat a, SimdFloat b)
{
    return {
               _mm_cmpneq_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
operator<(SimdFloat a, SimdFloat b)
{
    return {
               _mm_cmplt_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
operator<=(SimdFloat a, SimdFloat b)
{
    return {
               _mm_cmple_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
testBits(SimdFloat a)
{
    __m128i ia  = _mm_castps_si128(a.simdInternal_);
    __m128i res = _mm_andnot_si128( _mm_cmpeq_epi32(ia, _mm_setzero_si128()), _mm_cmpeq_epi32(ia, ia));

    return {
               _mm_castsi128_ps(res)
    };
}

static inline SimdFBool gmx_simdcall
operator&&(SimdFBool a, SimdFBool b)
{
    return {
               _mm_and_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
operator||(SimdFBool a, SimdFBool b)
{
    return {
               _mm_or_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdFBool a) { return _mm_movemask_ps(a.simdInternal_) != 0; }

static inline SimdFloat gmx_simdcall
selectByMask(SimdFloat a, SimdFBool mask)
{
    return {
               _mm_and_ps(a.simdInternal_, mask.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
selectByNotMask(SimdFloat a, SimdFBool mask)
{
    return {
               _mm_andnot_ps(mask.simdInternal_, a.simdInternal_)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdFloat gmx_simdcall
blend(SimdFloat a, SimdFloat b, SimdFBool sel)
{
    return {
               _mm_or_ps(_mm_andnot_ps(sel.simdInternal_, a.simdInternal_), _mm_and_ps(sel.simdInternal_, b.simdInternal_))
    };
}
#endif

static inline SimdFInt32 gmx_simdcall
operator<<(SimdFInt32 a, int n)
{
    return {
               _mm_slli_epi32(a.simdInternal_, n)
    };
}

static inline SimdFInt32 gmx_simdcall
operator>>(SimdFInt32 a, int n)
{
    return {
               _mm_srli_epi32(a.simdInternal_, n)
    };
}

static inline SimdFInt32 gmx_simdcall
operator&(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_and_si128(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
andNot(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_andnot_si128(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator|(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_or_si128(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator^(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_xor_si128(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator+(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_add_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator-(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_sub_epi32(a.simdInternal_, b.simdInternal_)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdFInt32 gmx_simdcall
operator*(SimdFInt32 a, SimdFInt32 b)
{
    __m128i a1 = _mm_srli_si128(a.simdInternal_, 4); // - a[3] a[2] a[1]
    __m128i b1 = _mm_srli_si128(b.simdInternal_, 4); // - b[3] b[2] b[1]
    __m128i c  = _mm_mul_epu32(a.simdInternal_, b.simdInternal_);
    __m128i c1 = _mm_mul_epu32(a1, b1);

    c  = _mm_shuffle_epi32(c, _MM_SHUFFLE(3, 1, 2, 0));  // - - a[2]*b[2] a[0]*b[0]
    c1 = _mm_shuffle_epi32(c1, _MM_SHUFFLE(3, 1, 2, 0)); // - - a[3]*b[3] a[1]*b[1]

    return {
               _mm_unpacklo_epi32(c, c1)
    };
}
#endif

static inline SimdFIBool gmx_simdcall
operator==(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_cmpeq_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
testBits(SimdFInt32 a)
{
    __m128i x   = a.simdInternal_;
    __m128i res = _mm_andnot_si128( _mm_cmpeq_epi32(x, _mm_setzero_si128()), _mm_cmpeq_epi32(x, x));

    return {
               res
    };
}

static inline SimdFIBool gmx_simdcall
operator<(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_cmplt_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator&&(SimdFIBool a, SimdFIBool b)
{
    return {
               _mm_and_si128(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator||(SimdFIBool a, SimdFIBool b)
{
    return {
               _mm_or_si128(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdFIBool a) { return _mm_movemask_epi8(a.simdInternal_) != 0; }

static inline SimdFInt32 gmx_simdcall
selectByMask(SimdFInt32 a, SimdFIBool mask)
{
    return {
               _mm_and_si128(a.simdInternal_, mask.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
selectByNotMask(SimdFInt32 a, SimdFIBool mask)
{
    return {
               _mm_andnot_si128(mask.simdInternal_, a.simdInternal_)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdFInt32 gmx_simdcall
blend(SimdFInt32 a, SimdFInt32 b, SimdFIBool sel)
{
    return {
               _mm_or_si128(_mm_andnot_si128(sel.simdInternal_, a.simdInternal_), _mm_and_si128(sel.simdInternal_, b.simdInternal_))
    };
}
#endif

static inline SimdFInt32 gmx_simdcall
cvtR2I(SimdFloat a)
{
    return {
               _mm_cvtps_epi32(a.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
cvttR2I(SimdFloat a)
{
    return {
               _mm_cvttps_epi32(a.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
cvtI2R(SimdFInt32 a)
{
    return {
               _mm_cvtepi32_ps(a.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
cvtB2IB(SimdFBool a)
{
    return {
               _mm_castps_si128(a.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
cvtIB2B(SimdFIBool a)
{
    return {
               _mm_castsi128_ps(a.simdInternal_)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE2_SIMD_FLOAT_H
