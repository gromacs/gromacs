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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_256_SIMD_DOUBLE_H

#include "config.h"

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <immintrin.h>

#include "gromacs/math/utilities.h"

#include "impl_x86_avx_256_simd_float.h"




namespace gmx
{

class SimdDouble
{
    public:
        SimdDouble() {}

        SimdDouble(double d) : simdInternal_(_mm256_set1_pd(d)) {}

        // Internal utility constructor to simplify return statements
        SimdDouble(__m256d simd) : simdInternal_(simd) {}

        __m256d  simdInternal_;
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

        SimdDBool(bool b) : simdInternal_(_mm256_castsi256_pd(_mm256_set1_epi32( b ? 0xFFFFFFFF : 0))) {}

        // Internal utility constructor to simplify return statements
        SimdDBool(__m256d simd) : simdInternal_(simd) {}

        __m256d  simdInternal_;
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
    assert(std::size_t(m) % 32 == 0);
    return {
               _mm256_load_pd(m)
    };
}

static inline void gmx_simdcall
store(double *m, SimdDouble a)
{
    assert(std::size_t(m) % 32 == 0);
    _mm256_store_pd(m, a.simdInternal_);
}

static inline SimdDouble gmx_simdcall
simdLoadU(const double *m)
{
    return {
               _mm256_loadu_pd(m)
    };
}

static inline void gmx_simdcall
storeU(double *m, SimdDouble a)
{
    _mm256_storeu_pd(m, a.simdInternal_);
}

static inline SimdDouble gmx_simdcall
setZeroD()
{
    return {
               _mm256_setzero_pd()
    };
}

static inline SimdDInt32 gmx_simdcall
simdLoadDI(const std::int32_t * m)
{
    assert(std::size_t(m) % 16 == 0);
    return {
               _mm_load_si128(reinterpret_cast<const __m128i *>(m))
    };
}

static inline void gmx_simdcall
store(std::int32_t * m, SimdDInt32 a)
{
    assert(std::size_t(m) % 16 == 0);
    _mm_store_si128(reinterpret_cast<__m128i *>(m), a.simdInternal_);
}

static inline SimdDInt32 gmx_simdcall
simdLoadUDI(const std::int32_t *m)
{
    return {
               _mm_loadu_si128(reinterpret_cast<const __m128i *>(m))
    };
}

static inline void gmx_simdcall
storeU(std::int32_t * m, SimdDInt32 a)
{
    _mm_storeu_si128(reinterpret_cast<__m128i *>(m), a.simdInternal_);
}

static inline SimdDInt32 gmx_simdcall
setZeroDI()
{
    return {
               _mm_setzero_si128()
    };
}

template<int index>
static inline std::int32_t gmx_simdcall
extract(SimdDInt32 a)
{
    return _mm_extract_epi32(a.simdInternal_, index);
}

static inline SimdDouble gmx_simdcall
operator&(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_and_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
andNot(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_andnot_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator|(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_or_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator^(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_xor_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator+(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_add_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator-(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_sub_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator-(SimdDouble x)
{
    return {
               _mm256_xor_pd(x.simdInternal_, _mm256_set1_pd(GMX_DOUBLE_NEGZERO))
    };
}

static inline SimdDouble gmx_simdcall
operator*(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_mul_pd(a.simdInternal_, b.simdInternal_)
    };
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
static inline SimdDouble gmx_simdcall
fma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm256_add_pd(_mm256_mul_pd(a.simdInternal_, b.simdInternal_), c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm256_sub_pd(_mm256_mul_pd(a.simdInternal_, b.simdInternal_), c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fnma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm256_sub_pd(c.simdInternal_, _mm256_mul_pd(a.simdInternal_, b.simdInternal_))
    };
}

static inline SimdDouble gmx_simdcall
fnms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm256_sub_pd(_mm256_setzero_pd(), _mm256_add_pd(_mm256_mul_pd(a.simdInternal_, b.simdInternal_), c.simdInternal_))
    };
}
#endif

static inline SimdDouble gmx_simdcall
rsqrt(SimdDouble x)
{
    return {
               _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x.simdInternal_)))
    };
}

static inline SimdDouble gmx_simdcall
rcp(SimdDouble x)
{
    return {
               _mm256_cvtps_pd(_mm_rcp_ps(_mm256_cvtpd_ps(x.simdInternal_)))
    };
}

static inline SimdDouble gmx_simdcall
maskAdd(SimdDouble a, SimdDouble b, SimdDBool m)
{
    return {
               _mm256_add_pd(a.simdInternal_, _mm256_and_pd(b.simdInternal_, m.simdInternal_))
    };
}

static inline SimdDouble gmx_simdcall
maskzMul(SimdDouble a, SimdDouble b, SimdDBool m)
{
    return {
               _mm256_and_pd(_mm256_mul_pd(a.simdInternal_, b.simdInternal_), m.simdInternal_)
    };
}

static inline SimdDouble
maskzFma(SimdDouble a, SimdDouble b, SimdDouble c, SimdDBool m)
{
    return {
               _mm256_and_pd(_mm256_add_pd(_mm256_mul_pd(a.simdInternal_, b.simdInternal_), c.simdInternal_), m.simdInternal_)
    };
}

static inline SimdDouble
maskzRsqrt(SimdDouble x, SimdDBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = _mm256_blendv_pd(_mm256_set1_pd(1.0), x.simdInternal_, m.simdInternal_);
#endif
    return {
               _mm256_and_pd(_mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x.simdInternal_))), m.simdInternal_)
    };
}

static inline SimdDouble
maskzRcp(SimdDouble x, SimdDBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = _mm256_blendv_pd(_mm256_set1_pd(1.0), x.simdInternal_, m.simdInternal_);
#endif
    return {
               _mm256_and_pd(_mm256_cvtps_pd(_mm_rcp_ps(_mm256_cvtpd_ps(x.simdInternal_))), m.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
abs(SimdDouble x)
{
    return {
               _mm256_andnot_pd( _mm256_set1_pd(GMX_DOUBLE_NEGZERO), x.simdInternal_ )
    };
}

static inline SimdDouble gmx_simdcall
max(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_max_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
min(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_min_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
round(SimdDouble x)
{
    return {
               _mm256_round_pd(x.simdInternal_, _MM_FROUND_NINT)
    };
}

static inline SimdDouble gmx_simdcall
trunc(SimdDouble x)
{
    return {
               _mm256_round_pd(x.simdInternal_, _MM_FROUND_TRUNC)
    };
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
static inline SimdDouble
frexp(SimdDouble value, SimdDInt32 * exponent)
{
    const __m256d exponentMask      = _mm256_castsi256_pd( _mm256_set1_epi64x(0x7FF0000000000000LL));
    const __m256d mantissaMask      = _mm256_castsi256_pd( _mm256_set1_epi64x(0x800FFFFFFFFFFFFFLL));
    const __m256d half              = _mm256_set1_pd(0.5);
    const __m128i exponentBias      = _mm_set1_epi32(1022); // add 1 to make our definition identical to frexp()
    __m256i       iExponent;
    __m128i       iExponentLow, iExponentHigh;

    iExponent               = _mm256_castpd_si256(_mm256_and_pd(value.simdInternal_, exponentMask));
    iExponentHigh           = _mm256_extractf128_si256(iExponent, 0x1);
    iExponentLow            = _mm256_castsi256_si128(iExponent);
    iExponentLow            = _mm_srli_epi64(iExponentLow, 52);
    iExponentHigh           = _mm_srli_epi64(iExponentHigh, 52);
    iExponentLow            = _mm_shuffle_epi32(iExponentLow, _MM_SHUFFLE(1, 1, 2, 0));
    iExponentHigh           = _mm_shuffle_epi32(iExponentHigh, _MM_SHUFFLE(2, 0, 1, 1));
    iExponentLow            = _mm_or_si128(iExponentLow, iExponentHigh);
    exponent->simdInternal_ = _mm_sub_epi32(iExponentLow, exponentBias);

    return {
               _mm256_or_pd(_mm256_and_pd(value.simdInternal_, mantissaMask), half)
    };
}

template <MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble
ldexp(SimdDouble value, SimdDInt32 exponent)
{
    const __m128i exponentBias = _mm_set1_epi32(1023);
    __m128i       iExponentLow, iExponentHigh;
    __m256d       fExponent;

    iExponentLow  = _mm_add_epi32(exponent.simdInternal_, exponentBias);

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponentLow  = _mm_max_epi32(iExponentLow, _mm_setzero_si128());
    }

    iExponentHigh = _mm_shuffle_epi32(iExponentLow, _MM_SHUFFLE(3, 3, 2, 2));
    iExponentLow  = _mm_shuffle_epi32(iExponentLow, _MM_SHUFFLE(1, 1, 0, 0));
    iExponentHigh = _mm_slli_epi64(iExponentHigh, 52);
    iExponentLow  = _mm_slli_epi64(iExponentLow, 52);
    fExponent     = _mm256_castsi256_pd(_mm256_insertf128_si256(_mm256_castsi128_si256(iExponentLow), iExponentHigh, 0x1));
    return {
               _mm256_mul_pd(value.simdInternal_, fExponent)
    };
}
#endif

static inline double gmx_simdcall
reduce(SimdDouble a)
{
    __m128d a0, a1;
    a.simdInternal_ = _mm256_add_pd(a.simdInternal_, _mm256_permute_pd(a.simdInternal_, 0b0101 ));
    a0              = _mm256_castpd256_pd128(a.simdInternal_);
    a1              = _mm256_extractf128_pd(a.simdInternal_, 0x1);
    a0              = _mm_add_sd(a0, a1);

    return *reinterpret_cast<double *>(&a0);
}

static inline SimdDBool gmx_simdcall
operator==(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_cmp_pd(a.simdInternal_, b.simdInternal_, _CMP_EQ_OQ)
    };
}

static inline SimdDBool gmx_simdcall
operator!=(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_cmp_pd(a.simdInternal_, b.simdInternal_, _CMP_NEQ_OQ)
    };
}

static inline SimdDBool gmx_simdcall
operator<(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_cmp_pd(a.simdInternal_, b.simdInternal_, _CMP_LT_OQ)
    };
}

static inline SimdDBool gmx_simdcall
operator<=(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_cmp_pd(a.simdInternal_, b.simdInternal_, _CMP_LE_OQ)
    };
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
static inline SimdDBool gmx_simdcall
testBits(SimdDouble a)
{
    // Do an or of the low/high 32 bits of each double (so the data is replicated),
    // and then use the same algorithm as we use for single precision.
    __m256 tst = _mm256_castpd_ps(a.simdInternal_);

    tst = _mm256_or_ps(tst, _mm256_permute_ps(tst, _MM_SHUFFLE(2, 3, 0, 1)));
    tst = _mm256_cvtepi32_ps(_mm256_castps_si256(tst));

    return {
               _mm256_castps_pd(_mm256_cmp_ps(tst, _mm256_setzero_ps(), _CMP_NEQ_OQ))
    };
}
#endif

static inline SimdDBool gmx_simdcall
operator&&(SimdDBool a, SimdDBool b)
{
    return {
               _mm256_and_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDBool gmx_simdcall
operator||(SimdDBool a, SimdDBool b)
{
    return {
               _mm256_or_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdDBool a) { return _mm256_movemask_pd(a.simdInternal_) != 0; }

static inline SimdDouble gmx_simdcall
selectByMask(SimdDouble a, SimdDBool mask)
{
    return {
               _mm256_and_pd(a.simdInternal_, mask.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
selectByNotMask(SimdDouble a, SimdDBool mask)
{
    return {
               _mm256_andnot_pd(mask.simdInternal_, a.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
blend(SimdDouble a, SimdDouble b, SimdDBool sel)
{
    return {
               _mm256_blendv_pd(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

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

static inline SimdDInt32 gmx_simdcall
operator*(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_mullo_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDIBool gmx_simdcall
operator==(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_cmpeq_epi32(a.simdInternal_, b.simdInternal_)
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
testBits(SimdDInt32 a)
{
    __m128i x   = a.simdInternal_;
    __m128i res = _mm_andnot_si128( _mm_cmpeq_epi32(x, _mm_setzero_si128()), _mm_cmpeq_epi32(x, x));

    return {
               res
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
anyTrue(SimdDIBool a) { return _mm_movemask_epi8(_mm_shuffle_epi32(a.simdInternal_, _MM_SHUFFLE(1, 0, 1, 0))) != 0; }

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

static inline SimdDInt32 gmx_simdcall
blend(SimdDInt32 a, SimdDInt32 b, SimdDIBool sel)
{
    return {
               _mm_blendv_epi8(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
cvtR2I(SimdDouble a)
{
    return {
               _mm256_cvtpd_epi32(a.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
cvttR2I(SimdDouble a)
{
    return {
               _mm256_cvttpd_epi32(a.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
cvtI2R(SimdDInt32 a)
{
    return {
               _mm256_cvtepi32_pd(a.simdInternal_)
    };
}

static inline SimdDIBool gmx_simdcall
cvtB2IB(SimdDBool a)
{
    __m128i a1 = _mm256_extractf128_si256(_mm256_castpd_si256(a.simdInternal_), 0x1);
    __m128i a0 = _mm256_castsi256_si128(_mm256_castpd_si256(a.simdInternal_));
    a0 = _mm_shuffle_epi32(a0, _MM_SHUFFLE(2, 0, 2, 0));
    a1 = _mm_shuffle_epi32(a1, _MM_SHUFFLE(2, 0, 2, 0));

    return {
               _mm_blend_epi16(a0, a1, 0xF0)
    };
}

static inline SimdDBool gmx_simdcall
cvtIB2B(SimdDIBool a)
{
    __m128d lo = _mm_castsi128_pd(_mm_unpacklo_epi32(a.simdInternal_, a.simdInternal_));
    __m128d hi = _mm_castsi128_pd(_mm_unpackhi_epi32(a.simdInternal_, a.simdInternal_));

    return {
               _mm256_insertf128_pd(_mm256_castpd128_pd256(lo), hi, 0x1)
    };
}

static inline void gmx_simdcall
cvtF2DD(SimdFloat f, SimdDouble *d0, SimdDouble *d1)
{
    d0->simdInternal_ = _mm256_cvtps_pd(_mm256_castps256_ps128(f.simdInternal_));
    d1->simdInternal_ = _mm256_cvtps_pd(_mm256_extractf128_ps(f.simdInternal_, 0x1));
}

static inline SimdFloat gmx_simdcall
cvtDD2F(SimdDouble d0, SimdDouble d1)
{
    __m128 f0 = _mm256_cvtpd_ps(d0.simdInternal_);
    __m128 f1 = _mm256_cvtpd_ps(d1.simdInternal_);
    return {
               _mm256_insertf128_ps(_mm256_castps128_ps256(f0), f1, 0x1)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_256_SIMD_DOUBLE_H
