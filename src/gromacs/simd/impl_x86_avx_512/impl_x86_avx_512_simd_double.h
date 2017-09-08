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

#ifndef GMX_SIMD_IMPL_X86_AVX_512_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_512_SIMD_DOUBLE_H

#include "config.h"

#include <cassert>
#include <cstdint>

#include <immintrin.h>

#include "gromacs/math/utilities.h"
#include "gromacs/utility/basedefinitions.h"

#include "impl_x86_avx_512_general.h"
#include "impl_x86_avx_512_simd_float.h"

namespace gmx
{

class SimdDouble
{
    public:
        SimdDouble() {}

        SimdDouble(double d) : simdInternal_(_mm512_set1_pd(d)) {}

        // Internal utility constructor to simplify return statements
        SimdDouble(__m512d simd) : simdInternal_(simd) {}

        __m512d  simdInternal_;
};

class SimdDInt32
{
    public:
        SimdDInt32() {}

        SimdDInt32(std::int32_t i) : simdInternal_(_mm256_set1_epi32(i)) {}

        // Internal utility constructor to simplify return statements
        SimdDInt32(__m256i simd) : simdInternal_(simd) {}

        __m256i  simdInternal_;
};

class SimdDBool
{
    public:
        SimdDBool() {}

        // Internal utility constructor to simplify return statements
        SimdDBool(__mmask8 simd) : simdInternal_(simd) {}

        __mmask8  simdInternal_;
};

class SimdDIBool
{
    public:
        SimdDIBool() {}

        // Internal utility constructor to simplify return statements
        SimdDIBool(__mmask16 simd) : simdInternal_(simd) {}

        __mmask16  simdInternal_;
};

static inline SimdDouble gmx_simdcall
simdLoad(const double *m)
{
    assert(std::size_t(m) % 64 == 0);
    return {
               _mm512_load_pd(m)
    };
}

static inline void gmx_simdcall
store(double *m, SimdDouble a)
{
    assert(std::size_t(m) % 64 == 0);
    _mm512_store_pd(m, a.simdInternal_);
}

static inline SimdDouble gmx_simdcall
simdLoadU(const double *m)
{
    return {
               _mm512_loadu_pd(m)
    };
}

static inline void gmx_simdcall
storeU(double *m, SimdDouble a)
{
    _mm512_storeu_pd(m, a.simdInternal_);
}

static inline SimdDouble gmx_simdcall
setZeroD()
{
    return {
               _mm512_setzero_pd()
    };
}

static inline SimdDInt32 gmx_simdcall
simdLoadDI(const std::int32_t * m)
{
    assert(std::size_t(m) % 32 == 0);
    return {
               _mm256_load_si256(reinterpret_cast<const __m256i *>(m))
    };
}

static inline void gmx_simdcall
store(std::int32_t * m, SimdDInt32 a)
{
    assert(std::size_t(m) % 32 == 0);
    _mm256_store_si256(reinterpret_cast<__m256i *>(m), a.simdInternal_);
}

static inline SimdDInt32 gmx_simdcall
simdLoadUDI(const std::int32_t *m)
{
    return {
               _mm256_loadu_si256(reinterpret_cast<const __m256i *>(m))
    };
}

static inline void gmx_simdcall
storeU(std::int32_t * m, SimdDInt32 a)
{
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(m), a.simdInternal_);
}

static inline SimdDInt32 gmx_simdcall
setZeroDI()
{
    return {
               _mm256_setzero_si256()
    };
}

static inline SimdDouble gmx_simdcall
operator&(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_castsi512_pd(_mm512_and_epi32(_mm512_castpd_si512(a.simdInternal_), _mm512_castpd_si512(b.simdInternal_)))
    };
}

static inline SimdDouble gmx_simdcall
andNot(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_castsi512_pd(_mm512_andnot_epi32(_mm512_castpd_si512(a.simdInternal_), _mm512_castpd_si512(b.simdInternal_)))
    };
}

static inline SimdDouble gmx_simdcall
operator|(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_castsi512_pd(_mm512_or_epi32(_mm512_castpd_si512(a.simdInternal_), _mm512_castpd_si512(b.simdInternal_)))
    };
}

static inline SimdDouble gmx_simdcall
operator^(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_castsi512_pd(_mm512_xor_epi32(_mm512_castpd_si512(a.simdInternal_), _mm512_castpd_si512(b.simdInternal_)))
    };
}

static inline SimdDouble gmx_simdcall
operator+(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_add_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator-(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_sub_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator-(SimdDouble x)
{
    return {
               _mm512_castsi512_pd(_mm512_xor_epi32(_mm512_castpd_si512(x.simdInternal_), _mm512_castpd_si512(_mm512_set1_pd(GMX_DOUBLE_NEGZERO))))
    };
}

static inline SimdDouble gmx_simdcall
operator*(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_mul_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm512_fmadd_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm512_fmsub_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fnma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm512_fnmadd_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fnms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm512_fnmsub_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

// Override for AVX-512-KNL
#if GMX_SIMD_X86_AVX_512
static inline SimdDouble gmx_simdcall
rsqrt(SimdDouble x)
{
    return {
               _mm512_rsqrt14_pd(x.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
rcp(SimdDouble x)
{
    return {
               _mm512_rcp14_pd(x.simdInternal_)
    };
}
#endif

static inline SimdDouble gmx_simdcall
maskAdd(SimdDouble a, SimdDouble b, SimdDBool m)
{
    return {
               _mm512_mask_add_pd(a.simdInternal_, m.simdInternal_, a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
maskzMul(SimdDouble a, SimdDouble b, SimdDBool m)
{
    return {
               _mm512_maskz_mul_pd(m.simdInternal_, a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
maskzFma(SimdDouble a, SimdDouble b, SimdDouble c, SimdDBool m)
{
    return {
               _mm512_maskz_fmadd_pd(m.simdInternal_, a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

// Override for AVX-512-KNL
#if GMX_SIMD_X86_AVX_512
static inline SimdDouble gmx_simdcall
maskzRsqrt(SimdDouble x, SimdDBool m)
{
    return {
               _mm512_maskz_rsqrt14_pd(m.simdInternal_, x.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
maskzRcp(SimdDouble x, SimdDBool m)
{
    return {
               _mm512_maskz_rcp14_pd(m.simdInternal_, x.simdInternal_)
    };
}
#endif

static inline SimdDouble gmx_simdcall
abs(SimdDouble x)
{
    return {
               _mm512_castsi512_pd(_mm512_andnot_epi32(_mm512_castpd_si512(_mm512_set1_pd(GMX_DOUBLE_NEGZERO)), _mm512_castpd_si512(x.simdInternal_)))
    };
}

static inline SimdDouble gmx_simdcall
max(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_max_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
min(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_min_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
round(SimdDouble x)
{
    return {
               _mm512_roundscale_pd(x.simdInternal_, 0)
    };
}

static inline SimdDouble gmx_simdcall
trunc(SimdDouble x)
{
#if defined(__INTEL_COMPILER) || defined(__ECC)
    return {
               _mm512_trunc_pd(x.simdInternal_)
    };
#else
    return {
               _mm512_cvtepi32_pd(_mm512_cvttpd_epi32(x.simdInternal_))
    };
#endif
}

static inline SimdDouble
frexp(SimdDouble value, SimdDInt32 * exponent)
{
    __m512d rExponent = _mm512_getexp_pd(value.simdInternal_);
    __m256i iExponent = _mm512_cvtpd_epi32(rExponent);

    exponent->simdInternal_ = _mm256_add_epi32(iExponent, _mm256_set1_epi32(1));

    return {
               _mm512_getmant_pd(value.simdInternal_, _MM_MANT_NORM_p5_1, _MM_MANT_SIGN_src)
    };
}

template <MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble
ldexp(SimdDouble value, SimdDInt32 exponent)
{
    const __m256i exponentBias = _mm256_set1_epi32(1023);
    __m256i       iExponent    = _mm256_add_epi32(exponent.simdInternal_, exponentBias);
    __m512i       iExponent512;

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponent = _mm256_max_epi32(iExponent, _mm256_setzero_si256());
    }

    iExponent512 = _mm512_permutexvar_epi32(_mm512_set_epi32(7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 0, 0), _mm512_castsi256_si512(iExponent));
    iExponent512 = _mm512_mask_slli_epi32(_mm512_setzero_epi32(), avx512Int2Mask(0xAAAA), iExponent512, 20);
    return _mm512_mul_pd(_mm512_castsi512_pd(iExponent512), value.simdInternal_);
}

static inline double gmx_simdcall
reduce(SimdDouble a)
{
    __m512d x = a.simdInternal_;
    x = _mm512_add_pd(x, _mm512_shuffle_f64x2(x, x, 0xEE));
    x = _mm512_add_pd(x, _mm512_shuffle_f64x2(x, x, 0x11));
    x = _mm512_add_pd(x, _mm512_permute_pd(x, 0x01));
    return *reinterpret_cast<double *>(&x);
}

static inline SimdDBool gmx_simdcall
operator==(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_cmp_pd_mask(a.simdInternal_, b.simdInternal_, _CMP_EQ_OQ)
    };
}

static inline SimdDBool gmx_simdcall
operator!=(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_cmp_pd_mask(a.simdInternal_, b.simdInternal_, _CMP_NEQ_OQ)
    };
}

static inline SimdDBool gmx_simdcall
operator<(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_cmp_pd_mask(a.simdInternal_, b.simdInternal_, _CMP_LT_OQ)
    };
}

static inline SimdDBool gmx_simdcall
operator<=(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_cmp_pd_mask(a.simdInternal_, b.simdInternal_, _CMP_LE_OQ)
    };
}

static inline SimdDBool gmx_simdcall
testBits(SimdDouble a)
{
    return {
               _mm512_test_epi64_mask(_mm512_castpd_si512(a.simdInternal_), _mm512_castpd_si512(a.simdInternal_))
    };
}

static inline SimdDBool gmx_simdcall
operator&&(SimdDBool a, SimdDBool b)
{
    return {
               static_cast<__mmask8>(_mm512_kand(a.simdInternal_, b.simdInternal_))
    };
}

static inline SimdDBool gmx_simdcall
operator||(SimdDBool a, SimdDBool b)
{
    return {
               static_cast<__mmask8>(_mm512_kor(a.simdInternal_, b.simdInternal_))
    };
}

static inline bool gmx_simdcall
anyTrue(SimdDBool a)
{
    return ( avx512Mask2Int(a.simdInternal_) != 0);
}

static inline SimdDouble gmx_simdcall
selectByMask(SimdDouble a, SimdDBool m)
{
    return {
               _mm512_mask_mov_pd(_mm512_setzero_pd(), m.simdInternal_, a.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
selectByNotMask(SimdDouble a, SimdDBool m)
{
    return {
               _mm512_mask_mov_pd(a.simdInternal_, m.simdInternal_, _mm512_setzero_pd())
    };
}

static inline SimdDouble gmx_simdcall
blend(SimdDouble a, SimdDouble b, SimdDBool sel)
{
    return {
               _mm512_mask_blend_pd(sel.simdInternal_, a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator<<(SimdDInt32 a, int n)
{
    return {
               _mm256_slli_epi32(a.simdInternal_, n)
    };
}

static inline SimdDInt32 gmx_simdcall
operator>>(SimdDInt32 a, int n)
{
    return {
               _mm256_srli_epi32(a.simdInternal_, n)
    };
}

static inline SimdDInt32 gmx_simdcall
operator&(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm256_and_si256(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
andNot(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm256_andnot_si256(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator|(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm256_or_si256(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator^(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm256_xor_si256(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator+(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm256_add_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator-(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm256_sub_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator*(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm256_mullo_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDIBool gmx_simdcall
operator==(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm512_mask_cmp_epi32_mask(avx512Int2Mask(0xFF), _mm512_castsi256_si512(a.simdInternal_), _mm512_castsi256_si512(b.simdInternal_), _MM_CMPINT_EQ)
    };
}

static inline SimdDIBool gmx_simdcall
testBits(SimdDInt32 a)
{
    return {
               _mm512_mask_test_epi32_mask(avx512Int2Mask(0xFF), _mm512_castsi256_si512(a.simdInternal_), _mm512_castsi256_si512(a.simdInternal_))
    };
}

static inline SimdDIBool gmx_simdcall
operator<(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm512_mask_cmp_epi32_mask(avx512Int2Mask(0xFF), _mm512_castsi256_si512(a.simdInternal_), _mm512_castsi256_si512(b.simdInternal_), _MM_CMPINT_LT)
    };
}

static inline SimdDIBool gmx_simdcall
operator&&(SimdDIBool a, SimdDIBool b)
{
    return {
               _mm512_kand(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDIBool gmx_simdcall
operator||(SimdDIBool a, SimdDIBool b)
{
    return {
               _mm512_kor(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdDIBool a)
{
    return ( avx512Mask2Int(a.simdInternal_) & 0xFF) != 0;
}

static inline SimdDInt32 gmx_simdcall
selectByMask(SimdDInt32 a, SimdDIBool m)
{
    return {
               _mm512_castsi512_si256(_mm512_mask_mov_epi32(_mm512_setzero_si512(), m.simdInternal_, _mm512_castsi256_si512(a.simdInternal_)))
    };
}

static inline SimdDInt32 gmx_simdcall
selectByNotMask(SimdDInt32 a, SimdDIBool m)
{
    return {
               _mm512_castsi512_si256(_mm512_mask_mov_epi32(_mm512_castsi256_si512(a.simdInternal_), m.simdInternal_, _mm512_setzero_si512()))
    };
}

static inline SimdDInt32 gmx_simdcall
blend(SimdDInt32 a, SimdDInt32 b, SimdDIBool sel)
{
    return {
               _mm512_castsi512_si256(_mm512_mask_blend_epi32(sel.simdInternal_, _mm512_castsi256_si512(a.simdInternal_), _mm512_castsi256_si512(b.simdInternal_)))
    };
}

static inline SimdDInt32 gmx_simdcall
cvtR2I(SimdDouble a)
{
    return {
               _mm512_cvtpd_epi32(a.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
cvttR2I(SimdDouble a)
{
    return {
               _mm512_cvttpd_epi32(a.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
cvtI2R(SimdDInt32 a)
{
    return {
               _mm512_cvtepi32_pd(a.simdInternal_)
    };
}

static inline SimdDIBool gmx_simdcall
cvtB2IB(SimdDBool a)
{
    return {
               a.simdInternal_
    };
}

static inline SimdDBool gmx_simdcall
cvtIB2B(SimdDIBool a)
{
    return {
               static_cast<__mmask8>(a.simdInternal_)
    };
}

static inline void gmx_simdcall
cvtF2DD(SimdFloat f, SimdDouble *d0, SimdDouble *d1)
{
    d0->simdInternal_ = _mm512_cvtps_pd(_mm512_castps512_ps256(f.simdInternal_));
    d1->simdInternal_ = _mm512_cvtps_pd(_mm512_castps512_ps256(_mm512_shuffle_f32x4(f.simdInternal_, f.simdInternal_, 0xEE)));
}

static inline SimdFloat gmx_simdcall
cvtDD2F(SimdDouble d0, SimdDouble d1)
{
    __m512 f0 = _mm512_castps256_ps512(_mm512_cvtpd_ps(d0.simdInternal_));
    __m512 f1 = _mm512_castps256_ps512(_mm512_cvtpd_ps(d1.simdInternal_));
    return {
               _mm512_shuffle_f32x4(f0, f1, 0x44)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_512_SIMD_DOUBLE_H
