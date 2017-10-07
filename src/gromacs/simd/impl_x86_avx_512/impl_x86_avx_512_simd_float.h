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

#ifndef GMX_SIMD_IMPL_X86_AVX_512_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_512_SIMD_FLOAT_H

#include "config.h"

#include <cassert>
#include <cstdint>

#include <immintrin.h>

#include "gromacs/math/utilities.h"
#include "gromacs/utility/real.h"

#include "impl_x86_avx_512_general.h"

namespace gmx
{

class SimdFloat
{
    public:
        SimdFloat() {}

        SimdFloat(float f) : simdInternal_(_mm512_set1_ps(f)) {}

        // Internal utility constructor to simplify return statements
        SimdFloat(__m512 simd) : simdInternal_(simd) {}

        __m512  simdInternal_;
};

class SimdFInt32
{
    public:
        SimdFInt32() {}

        SimdFInt32(std::int32_t i) : simdInternal_(_mm512_set1_epi32(i)) {}

        // Internal utility constructor to simplify return statements
        SimdFInt32(__m512i simd) : simdInternal_(simd) {}

        __m512i  simdInternal_;
};

class SimdFBool
{
    public:
        SimdFBool() {}

        // Internal utility constructor to simplify return statements
        SimdFBool(__mmask16 simd) : simdInternal_(simd) {}

        __mmask16  simdInternal_;
};

class SimdFIBool
{
    public:
        SimdFIBool() {}

        // Internal utility constructor to simplify return statements
        SimdFIBool(__mmask16 simd) : simdInternal_(simd) {}

        __mmask16  simdInternal_;
};

static inline SimdFloat gmx_simdcall
simdLoad(const float *m, SimdFloatTag = {})
{
    assert(std::size_t(m) % 64 == 0);
    return {
               _mm512_load_ps(m)
    };
}

static inline void gmx_simdcall
store(float *m, SimdFloat a)
{
    assert(std::size_t(m) % 64 == 0);
    _mm512_store_ps(m, a.simdInternal_);
}

static inline SimdFloat gmx_simdcall
simdLoadU(const float *m, SimdFloatTag = {})
{
    return {
               _mm512_loadu_ps(m)
    };
}

static inline void gmx_simdcall
storeU(float *m, SimdFloat a)
{
    _mm512_storeu_ps(m, a.simdInternal_);
}

static inline SimdFloat gmx_simdcall
setZeroF()
{
    return {
               _mm512_setzero_ps()
    };
}

static inline SimdFInt32 gmx_simdcall
simdLoad(const std::int32_t * m, SimdFInt32Tag)
{
    assert(std::size_t(m) % 64 == 0);
    return {
               _mm512_load_si512(m)
    };
}

static inline void gmx_simdcall
store(std::int32_t * m, SimdFInt32 a)
{
    assert(std::size_t(m) % 64 == 0);
    _mm512_store_si512(m, a.simdInternal_);
}

static inline SimdFInt32 gmx_simdcall
simdLoadU(const std::int32_t *m, SimdFInt32Tag)
{
    return {
               _mm512_loadu_si512(m)
    };
}

static inline void gmx_simdcall
storeU(std::int32_t * m, SimdFInt32 a)
{
    _mm512_storeu_si512(m, a.simdInternal_);
}

static inline SimdFInt32 gmx_simdcall
setZeroFI()
{
    return {
               _mm512_setzero_si512()
    };
}


static inline SimdFloat gmx_simdcall
operator&(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(a.simdInternal_), _mm512_castps_si512(b.simdInternal_)))
    };
}

static inline SimdFloat gmx_simdcall
andNot(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_castsi512_ps(_mm512_andnot_epi32(_mm512_castps_si512(a.simdInternal_), _mm512_castps_si512(b.simdInternal_)))
    };
}

static inline SimdFloat gmx_simdcall
operator|(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_castsi512_ps(_mm512_or_epi32(_mm512_castps_si512(a.simdInternal_), _mm512_castps_si512(b.simdInternal_)))
    };
}

static inline SimdFloat gmx_simdcall
operator^(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(a.simdInternal_), _mm512_castps_si512(b.simdInternal_)))
    };
}

static inline SimdFloat gmx_simdcall
operator+(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_add_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator-(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_sub_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator-(SimdFloat x)
{
    return {
               _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(x.simdInternal_), _mm512_castps_si512(_mm512_set1_ps(GMX_FLOAT_NEGZERO))))
    };
}

static inline SimdFloat gmx_simdcall
operator*(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_mul_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm512_fmadd_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm512_fmsub_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fnma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm512_fnmadd_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fnms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm512_fnmsub_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

// Override for AVX-512-KNL
#if GMX_SIMD_X86_AVX_512
static inline SimdFloat gmx_simdcall
rsqrt(SimdFloat x)
{
    return {
               _mm512_rsqrt14_ps(x.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
rcp(SimdFloat x)
{
    return {
               _mm512_rcp14_ps(x.simdInternal_)
    };
}
#endif

static inline SimdFloat gmx_simdcall
maskAdd(SimdFloat a, SimdFloat b, SimdFBool m)
{
    return {
               _mm512_mask_add_ps(a.simdInternal_, m.simdInternal_, a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
maskzMul(SimdFloat a, SimdFloat b, SimdFBool m)
{
    return {
               _mm512_maskz_mul_ps(m.simdInternal_, a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
maskzFma(SimdFloat a, SimdFloat b, SimdFloat c, SimdFBool m)
{
    return {
               _mm512_maskz_fmadd_ps(m.simdInternal_, a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

// Override for AVX-512-KNL
#if GMX_SIMD_X86_AVX_512
static inline SimdFloat gmx_simdcall
maskzRsqrt(SimdFloat x, SimdFBool m)
{
    return {
               _mm512_maskz_rsqrt14_ps(m.simdInternal_, x.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
maskzRcp(SimdFloat x, SimdFBool m)
{
    return {
               _mm512_maskz_rcp14_ps(m.simdInternal_, x.simdInternal_)
    };
}
#endif

static inline SimdFloat gmx_simdcall
abs(SimdFloat x)
{
    return {
               _mm512_castsi512_ps(_mm512_andnot_epi32(_mm512_castps_si512(_mm512_set1_ps(GMX_FLOAT_NEGZERO)), _mm512_castps_si512(x.simdInternal_)))
    };
}

static inline SimdFloat gmx_simdcall
max(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_max_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
min(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_min_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
round(SimdFloat x)
{
    return {
               _mm512_roundscale_ps(x.simdInternal_, 0)
    };
}

static inline SimdFloat gmx_simdcall
trunc(SimdFloat x)
{
#if defined(__INTEL_COMPILER) || defined(__ECC)
    return {
               _mm512_trunc_ps(x.simdInternal_)
    };
#else
    return {
               _mm512_cvtepi32_ps(_mm512_cvttps_epi32(x.simdInternal_))
    };
#endif
}

static inline SimdFloat gmx_simdcall
frexp(SimdFloat value, SimdFInt32 * exponent)
{
    __m512  rExponent = _mm512_getexp_ps(value.simdInternal_);
    __m512i iExponent =  _mm512_cvtps_epi32(rExponent);

    exponent->simdInternal_ = _mm512_add_epi32(iExponent, _mm512_set1_epi32(1));

    return {
               _mm512_getmant_ps(value.simdInternal_, _MM_MANT_NORM_p5_1, _MM_MANT_SIGN_src)
    };
}

template <MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall
ldexp(SimdFloat value, SimdFInt32 exponent)
{
    const __m512i exponentBias = _mm512_set1_epi32(127);
    __m512i       iExponent    =  _mm512_add_epi32(exponent.simdInternal_, exponentBias);

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponent = _mm512_max_epi32(iExponent, _mm512_setzero_epi32());
    }

    iExponent = _mm512_slli_epi32(iExponent, 23);

    return {
               _mm512_mul_ps(value.simdInternal_, _mm512_castsi512_ps(iExponent))
    };
}

static inline float gmx_simdcall
reduce(SimdFloat a)
{
    __m512 x = a.simdInternal_;
    x = _mm512_add_ps(x, _mm512_shuffle_f32x4(x, x, 0xEE));
    x = _mm512_add_ps(x, _mm512_shuffle_f32x4(x, x, 0x11));
    x = _mm512_add_ps(x, _mm512_permute_ps(x, 0xEE));
    x = _mm512_add_ps(x, _mm512_permute_ps(x, 0x11));
    return *reinterpret_cast<float *>(&x);
}

static inline SimdFBool gmx_simdcall
operator==(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_cmp_ps_mask(a.simdInternal_, b.simdInternal_, _CMP_EQ_OQ)
    };
}

static inline SimdFBool gmx_simdcall
operator!=(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_cmp_ps_mask(a.simdInternal_, b.simdInternal_, _CMP_NEQ_OQ)
    };
}

static inline SimdFBool gmx_simdcall
operator<(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_cmp_ps_mask(a.simdInternal_, b.simdInternal_, _CMP_LT_OQ)
    };
}

static inline SimdFBool gmx_simdcall
operator<=(SimdFloat a, SimdFloat b)
{
    return {
               _mm512_cmp_ps_mask(a.simdInternal_, b.simdInternal_, _CMP_LE_OQ)
    };
}

static inline SimdFBool gmx_simdcall
testBits(SimdFloat a)
{
    return {
               _mm512_test_epi32_mask( _mm512_castps_si512(a.simdInternal_), _mm512_castps_si512(a.simdInternal_) )
    };
}

static inline SimdFBool gmx_simdcall
operator&&(SimdFBool a, SimdFBool b)
{
    return {
               _mm512_kand(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
operator||(SimdFBool a, SimdFBool b)
{
    return {
               _mm512_kor(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdFBool a)
{
    return ( avx512Mask2Int(a.simdInternal_) != 0);
}

static inline SimdFloat gmx_simdcall
selectByMask(SimdFloat a, SimdFBool m)
{
    return {
               _mm512_mask_mov_ps(_mm512_setzero_ps(), m.simdInternal_, a.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
selectByNotMask(SimdFloat a, SimdFBool m)
{
    return {
               _mm512_mask_mov_ps(a.simdInternal_, m.simdInternal_, _mm512_setzero_ps())
    };
}

static inline SimdFloat gmx_simdcall
blend(SimdFloat a, SimdFloat b, SimdFBool sel)
{
    return {
               _mm512_mask_blend_ps(sel.simdInternal_, a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator<<(SimdFInt32 a, int n)
{
    return {
               _mm512_slli_epi32(a.simdInternal_, n)
    };
}

static inline SimdFInt32 gmx_simdcall
operator>>(SimdFInt32 a, int n)
{
    return {
               _mm512_srli_epi32(a.simdInternal_, n)
    };
}

static inline SimdFInt32 gmx_simdcall
operator&(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm512_and_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
andNot(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm512_andnot_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator|(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm512_or_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator^(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm512_xor_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator+(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm512_add_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator-(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm512_sub_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator*(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm512_mullo_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator==(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm512_cmp_epi32_mask(a.simdInternal_, b.simdInternal_, _MM_CMPINT_EQ)
    };
}

static inline SimdFIBool gmx_simdcall
testBits(SimdFInt32 a)
{
    return {
               _mm512_test_epi32_mask( a.simdInternal_, a.simdInternal_ )
    };
}

static inline SimdFIBool gmx_simdcall
operator<(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm512_cmp_epi32_mask(a.simdInternal_, b.simdInternal_, _MM_CMPINT_LT)
    };
}

static inline SimdFIBool gmx_simdcall
operator&&(SimdFIBool a, SimdFIBool b)
{
    return {
               _mm512_kand(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator||(SimdFIBool a, SimdFIBool b)
{
    return {
               _mm512_kor(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdFIBool a)
{
    return ( avx512Mask2Int(a.simdInternal_) != 0);
}

static inline SimdFInt32 gmx_simdcall
selectByMask(SimdFInt32 a, SimdFIBool m)
{
    return {
               _mm512_mask_mov_epi32(_mm512_setzero_epi32(), m.simdInternal_, a.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
selectByNotMask(SimdFInt32 a, SimdFIBool m)
{
    return {
               _mm512_mask_mov_epi32(a.simdInternal_, m.simdInternal_, _mm512_setzero_epi32())
    };
}

static inline SimdFInt32 gmx_simdcall
blend(SimdFInt32 a, SimdFInt32 b, SimdFIBool sel)
{
    return {
               _mm512_mask_blend_epi32(sel.simdInternal_, a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
cvtR2I(SimdFloat a)
{
    return {
               _mm512_cvtps_epi32(a.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
cvttR2I(SimdFloat a)
{
    return {
               _mm512_cvttps_epi32(a.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
cvtI2R(SimdFInt32 a)
{
    return {
               _mm512_cvtepi32_ps(a.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
cvtB2IB(SimdFBool a)
{
    return {
               a.simdInternal_
    };
}

static inline SimdFBool gmx_simdcall
cvtIB2B(SimdFIBool a)
{
    return {
               a.simdInternal_
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_512_SIMD_FLOAT_H
