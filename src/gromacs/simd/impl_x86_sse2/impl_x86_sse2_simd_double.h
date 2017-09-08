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
#ifndef GMX_SIMD_IMPL_X86_SSE2_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_X86_SSE2_SIMD_DOUBLE_H

#include "config.h"

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <emmintrin.h>

#include "gromacs/math/utilities.h"

#include "impl_x86_sse2_simd_float.h"

namespace gmx
{

class SimdDouble
{
    public:
        SimdDouble() {}

        SimdDouble(double d) : simdInternal_(_mm_set1_pd(d)) {}

        // Internal utility constructor to simplify return statements
        SimdDouble(__m128d simd) : simdInternal_(simd) {}

        __m128d  simdInternal_;
};

class SimdDInt32
{
    public:
        SimdDInt32() {}

        SimdDInt32(std::int32_t i) : simdInternal_(_mm_set1_epi32(i)) {}

        // Internal utility constructor to simplify return statements
        SimdDInt32(__m128i simd) : simdInternal_(simd) {}

        __m128i  simdInternal_;
};

class SimdDBool
{
    public:
        SimdDBool() {}

        SimdDBool(bool b) : simdInternal_(_mm_castsi128_pd(_mm_set1_epi32( b ? 0xFFFFFFFF : 0))) {}

        // Internal utility constructor to simplify return statements
        SimdDBool(__m128d simd) : simdInternal_(simd) {}

        __m128d  simdInternal_;
};

class SimdDIBool
{
    public:
        SimdDIBool() {}

        SimdDIBool(bool b) : simdInternal_(_mm_set1_epi32( b ? 0xFFFFFFFF : 0)) {}

        // Internal utility constructor to simplify return statements
        SimdDIBool(__m128i simd) : simdInternal_(simd) {}

        __m128i  simdInternal_;
};

static inline SimdDouble gmx_simdcall
simdLoad(const double *m)
{
    assert(std::size_t(m) % 16 == 0);
    return {
               _mm_load_pd(m)
    };
}

static inline void gmx_simdcall
store(double *m, SimdDouble a)
{
    assert(std::size_t(m) % 16 == 0);
    _mm_store_pd(m, a.simdInternal_);
}

static inline SimdDouble gmx_simdcall
simdLoadU(const double *m)
{
    return {
               _mm_loadu_pd(m)
    };
}

static inline void gmx_simdcall
storeU(double *m, SimdDouble a) { _mm_storeu_pd(m, a.simdInternal_); }

static inline SimdDouble gmx_simdcall
setZeroD()
{
    return {
               _mm_setzero_pd()
    };
}

static inline SimdDInt32 gmx_simdcall
simdLoadDI(const std::int32_t * m)
{
    assert(std::size_t(m) % 8 == 0);
    return {
               _mm_loadl_epi64(reinterpret_cast<const __m128i *>(m))
    };
}

static inline void gmx_simdcall
store(std::int32_t * m, SimdDInt32 a)
{
    assert(std::size_t(m) % 8 == 0);
    _mm_storel_epi64(reinterpret_cast<__m128i *>(m), a.simdInternal_);
}

static inline SimdDInt32 gmx_simdcall
simdLoadUDI(const std::int32_t *m)
{
    return {
               _mm_loadl_epi64(reinterpret_cast<const __m128i *>(m))
    };
}

static inline void gmx_simdcall
storeU(std::int32_t * m, SimdDInt32 a)
{
    _mm_storel_epi64(reinterpret_cast<__m128i *>(m), a.simdInternal_);
}

static inline SimdDInt32 gmx_simdcall
setZeroDI()
{
    return {
               _mm_setzero_si128()
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
template<int index>
static inline std::int32_t gmx_simdcall
extract(SimdDInt32 a)
{
    return _mm_cvtsi128_si32( _mm_srli_si128(a.simdInternal_, 4 * index) );
}
#endif

static inline SimdDouble gmx_simdcall
operator&(SimdDouble a, SimdDouble b)
{
    return {
               _mm_and_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
andNot(SimdDouble a, SimdDouble b)
{
    return {
               _mm_andnot_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator|(SimdDouble a, SimdDouble b)
{
    return {
               _mm_or_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator^(SimdDouble a, SimdDouble b)
{
    return {
               _mm_xor_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator+(SimdDouble a, SimdDouble b)
{
    return {
               _mm_add_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator-(SimdDouble a, SimdDouble b)
{
    return {
               _mm_sub_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator-(SimdDouble x)
{
    return {
               _mm_xor_pd(x.simdInternal_, _mm_set1_pd(GMX_DOUBLE_NEGZERO))
    };
}

static inline SimdDouble gmx_simdcall
operator*(SimdDouble a, SimdDouble b)
{
    return {
               _mm_mul_pd(a.simdInternal_, b.simdInternal_)
    };
}

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline SimdDouble gmx_simdcall
fma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm_add_pd(_mm_mul_pd(a.simdInternal_, b.simdInternal_), c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm_sub_pd(_mm_mul_pd(a.simdInternal_, b.simdInternal_), c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fnma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm_sub_pd(c.simdInternal_, _mm_mul_pd(a.simdInternal_, b.simdInternal_))
    };
}

static inline SimdDouble gmx_simdcall
fnms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm_sub_pd(_mm_setzero_pd(), _mm_add_pd(_mm_mul_pd(a.simdInternal_, b.simdInternal_), c.simdInternal_))
    };
}
#endif

static inline SimdDouble gmx_simdcall
rsqrt(SimdDouble x)
{
    return {
               _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(x.simdInternal_)))
    };
}

static inline SimdDouble gmx_simdcall
rcp(SimdDouble x)
{
    return {
               _mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(x.simdInternal_)))
    };
}

static inline SimdDouble gmx_simdcall
maskAdd(SimdDouble a, SimdDouble b, SimdDBool m)
{
    return {
               _mm_add_pd(a.simdInternal_, _mm_and_pd(b.simdInternal_, m.simdInternal_))
    };
}

static inline SimdDouble gmx_simdcall
maskzMul(SimdDouble a, SimdDouble b, SimdDBool m)
{
    return {
               _mm_and_pd(_mm_mul_pd(a.simdInternal_, b.simdInternal_), m.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
maskzFma(SimdDouble a, SimdDouble b, SimdDouble c, SimdDBool m)
{
    return {
               _mm_and_pd(_mm_add_pd(_mm_mul_pd(a.simdInternal_, b.simdInternal_), c.simdInternal_), m.simdInternal_)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdDouble gmx_simdcall
maskzRsqrt(SimdDouble x, SimdDBool m)
{
    // The result will always be correct since we mask the result with m, but
    // for debug builds we also want to make sure not to generate FP exceptions
#ifndef NDEBUG
    x.simdInternal_ = _mm_or_pd(_mm_andnot_pd(m.simdInternal_, _mm_set1_pd(1.0)), _mm_and_pd(m.simdInternal_, x.simdInternal_));
#endif
    return {
               _mm_and_pd(_mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(x.simdInternal_))), m.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
maskzRcp(SimdDouble x, SimdDBool m)
{
    // The result will always be correct since we mask the result with m, but
    // for debug builds we also want to make sure not to generate FP exceptions
#ifndef NDEBUG
    x.simdInternal_ = _mm_or_pd(_mm_andnot_pd(m.simdInternal_, _mm_set1_pd(1.0)), _mm_and_pd(m.simdInternal_, x.simdInternal_));
#endif
    return {
               _mm_and_pd(_mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(x.simdInternal_))), m.simdInternal_)
    };
}
#endif

static inline SimdDouble gmx_simdcall
abs(SimdDouble x)
{
    return {
               _mm_andnot_pd( _mm_set1_pd(GMX_DOUBLE_NEGZERO), x.simdInternal_ )
    };
}

static inline SimdDouble gmx_simdcall
max(SimdDouble a, SimdDouble b)
{
    return {
               _mm_max_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
min(SimdDouble a, SimdDouble b)
{
    return {
               _mm_min_pd(a.simdInternal_, b.simdInternal_)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdDouble gmx_simdcall
round(SimdDouble x)
{
    return {
               _mm_cvtepi32_pd( _mm_cvtpd_epi32(x.simdInternal_) )
    };
}

static inline SimdDouble gmx_simdcall
trunc(SimdDouble x)
{
    return {
               _mm_cvtepi32_pd( _mm_cvttpd_epi32(x.simdInternal_) )
    };
}

#endif

static inline SimdDouble
frexp(SimdDouble value, SimdDInt32 * exponent)
{
    // Don't use _mm_set1_epi64x() - on MSVC it is only supported for 64-bit builds
    const __m128d exponentMask = _mm_castsi128_pd( _mm_set_epi32(0x7FF00000, 0x00000000, 0x7FF00000, 0x00000000) );
    const __m128d mantissaMask = _mm_castsi128_pd( _mm_set_epi32(0x800FFFFF, 0xFFFFFFFF, 0x800FFFFF, 0xFFFFFFFF) );
    const __m128i exponentBias = _mm_set1_epi32(1022); // add 1 to make our definition identical to frexp()
    const __m128d half         = _mm_set1_pd(0.5);
    __m128i       iExponent;

    iExponent               = _mm_castpd_si128(_mm_and_pd(value.simdInternal_, exponentMask));
    iExponent               = _mm_sub_epi32(_mm_srli_epi64(iExponent, 52), exponentBias);
    iExponent               = _mm_shuffle_epi32(iExponent, _MM_SHUFFLE(3, 1, 2, 0) );
    exponent->simdInternal_ = iExponent;

    return {
               _mm_or_pd(_mm_and_pd(value.simdInternal_, mantissaMask), half)
    };
}

// Override for SSE4.1
#if GMX_SIMD_X86_SSE2
template <MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble
ldexp(SimdDouble value, SimdDInt32 exponent)
{
    const __m128i  exponentBias = _mm_set1_epi32(1023);
    __m128i        iExponent    = _mm_add_epi32(exponent.simdInternal_, exponentBias);

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponent = _mm_and_si128(iExponent, _mm_cmpgt_epi32(iExponent, _mm_setzero_si128()));
    }

    // After conversion integers will be in slot 0,1. Move them to 0,2 so
    // we can do a 64-bit shift and get them to the dp exponents.
    iExponent = _mm_shuffle_epi32(iExponent, _MM_SHUFFLE(3, 1, 2, 0));
    iExponent = _mm_slli_epi64(iExponent, 52);

    return {
               _mm_mul_pd(value.simdInternal_, _mm_castsi128_pd(iExponent))
    };
}
#endif

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline double gmx_simdcall
reduce(SimdDouble a)
{
    __m128d b = _mm_add_sd(a.simdInternal_, _mm_shuffle_pd(a.simdInternal_, a.simdInternal_, _MM_SHUFFLE2(1, 1)));
    return *reinterpret_cast<double *>(&b);
}
#endif

static inline SimdDBool gmx_simdcall
operator==(SimdDouble a, SimdDouble b)
{
    return {
               _mm_cmpeq_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDBool gmx_simdcall
operator!=(SimdDouble a, SimdDouble b)
{
    return {
               _mm_cmpneq_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDBool gmx_simdcall
operator<(SimdDouble a, SimdDouble b)
{
    return {
               _mm_cmplt_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDBool gmx_simdcall
operator<=(SimdDouble a, SimdDouble b)
{
    return {
               _mm_cmple_pd(a.simdInternal_, b.simdInternal_)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdDBool gmx_simdcall
testBits(SimdDouble a)
{
    __m128i ia  = _mm_castpd_si128(a.simdInternal_);
    __m128i res = _mm_andnot_si128( _mm_cmpeq_epi32(ia, _mm_setzero_si128()), _mm_cmpeq_epi32(ia, ia));

    // set each 64-bit element if low or high 32-bit part is set
    res = _mm_or_si128(res, _mm_shuffle_epi32(res, _MM_SHUFFLE(2, 3, 0, 1)));

    return {
               _mm_castsi128_pd(res)
    };
}
#endif

static inline SimdDBool gmx_simdcall
operator&&(SimdDBool a, SimdDBool b)
{
    return {
               _mm_and_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDBool gmx_simdcall
operator||(SimdDBool a, SimdDBool b)
{
    return {
               _mm_or_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdDBool a) { return _mm_movemask_pd(a.simdInternal_) != 0; }

static inline SimdDouble gmx_simdcall
selectByMask(SimdDouble a, SimdDBool mask)
{
    return {
               _mm_and_pd(a.simdInternal_, mask.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
selectByNotMask(SimdDouble a, SimdDBool mask)
{
    return {
               _mm_andnot_pd(mask.simdInternal_, a.simdInternal_)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdDouble gmx_simdcall
blend(SimdDouble a, SimdDouble b, SimdDBool sel)
{
    return {
               _mm_or_pd(_mm_andnot_pd(sel.simdInternal_, a.simdInternal_), _mm_and_pd(sel.simdInternal_, b.simdInternal_))
    };
}
#endif

static inline SimdDInt32 gmx_simdcall
operator<<(SimdDInt32 a, int n)
{
    return {
               _mm_slli_epi32(a.simdInternal_, n)
    };
}

static inline SimdDInt32 gmx_simdcall
operator>>(SimdDInt32 a, int n)
{
    return {
               _mm_srli_epi32(a.simdInternal_, n)
    };
}

static inline SimdDInt32 gmx_simdcall
operator&(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_and_si128(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
andNot(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_andnot_si128(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator|(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_or_si128(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator^(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_xor_si128(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator+(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_add_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator-(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_sub_epi32(a.simdInternal_, b.simdInternal_)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdDInt32 gmx_simdcall
operator*(SimdDInt32 a, SimdDInt32 b)
{

    __m128i tmpA = _mm_unpacklo_epi32(a.simdInternal_, _mm_setzero_si128()); // 0 a[1] 0 a[0]
    __m128i tmpB = _mm_unpacklo_epi32(b.simdInternal_, _mm_setzero_si128()); // 0 b[1] 0 b[0]

    __m128i tmpC  = _mm_mul_epu32(tmpA, tmpB);                               // 0 a[1]*b[1] 0 a[0]*b[0]

    return {
               _mm_shuffle_epi32(tmpC, _MM_SHUFFLE(3, 1, 2, 0))
    };
}
#endif

static inline SimdDIBool gmx_simdcall
operator==(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_cmpeq_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDIBool gmx_simdcall
testBits(SimdDInt32 a)
{
    __m128i x   = a.simdInternal_;
    __m128i res = _mm_andnot_si128( _mm_cmpeq_epi32(x, _mm_setzero_si128()), _mm_cmpeq_epi32(x, x));

    return {
               res
    };
}

static inline SimdDIBool gmx_simdcall
operator<(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_cmplt_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDIBool gmx_simdcall
operator&&(SimdDIBool a, SimdDIBool b)
{
    return {
               _mm_and_si128(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDIBool gmx_simdcall
operator||(SimdDIBool a, SimdDIBool b)
{
    return {
               _mm_or_si128(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdDIBool a)
{
    return _mm_movemask_epi8(_mm_shuffle_epi32(a.simdInternal_, _MM_SHUFFLE(1, 0, 1, 0))) != 0;
}

static inline SimdDInt32 gmx_simdcall
selectByMask(SimdDInt32 a, SimdDIBool mask)
{
    return {
               _mm_and_si128(a.simdInternal_, mask.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
selectByNotMask(SimdDInt32 a, SimdDIBool mask)
{
    return {
               _mm_andnot_si128(mask.simdInternal_, a.simdInternal_)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdDInt32 gmx_simdcall
blend(SimdDInt32 a, SimdDInt32 b, SimdDIBool sel)
{
    return {
               _mm_or_si128(_mm_andnot_si128(sel.simdInternal_, a.simdInternal_), _mm_and_si128(sel.simdInternal_, b.simdInternal_))
    };
}
#endif

static inline SimdDInt32 gmx_simdcall
cvtR2I(SimdDouble a)
{
    return {
               _mm_cvtpd_epi32(a.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
cvttR2I(SimdDouble a)
{
    return {
               _mm_cvttpd_epi32(a.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
cvtI2R(SimdDInt32 a)
{
    return {
               _mm_cvtepi32_pd(a.simdInternal_)
    };
}

static inline SimdDIBool gmx_simdcall
cvtB2IB(SimdDBool a)
{
    return {
               _mm_shuffle_epi32(_mm_castpd_si128(a.simdInternal_), _MM_SHUFFLE(2, 0, 2, 0))
    };
}

static inline SimdDBool gmx_simdcall
cvtIB2B(SimdDIBool a)
{
    return {
               _mm_castsi128_pd(_mm_shuffle_epi32(a.simdInternal_, _MM_SHUFFLE(1, 1, 0, 0)))
    };
}

static inline void gmx_simdcall
cvtF2DD(SimdFloat f, SimdDouble *d0, SimdDouble *d1)
{
    d0->simdInternal_ = _mm_cvtps_pd(f.simdInternal_);
    d1->simdInternal_ = _mm_cvtps_pd(_mm_movehl_ps(f.simdInternal_, f.simdInternal_));
}

static inline SimdFloat gmx_simdcall
cvtDD2F(SimdDouble d0, SimdDouble d1)
{
    return {
               _mm_movelh_ps(_mm_cvtpd_ps(d0.simdInternal_), _mm_cvtpd_ps(d1.simdInternal_))
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE2_SIMD_DOUBLE_H
