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

#ifndef GMX_SIMD_IMPL_X86_MIC_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_X86_MIC_SIMD_DOUBLE_H

#include "config.h"

#include <cassert>
#include <cstdint>

#include <immintrin.h>

#include "gromacs/math/utilities.h"
#include "gromacs/utility/basedefinitions.h"

#include "impl_x86_mic_simd_float.h"

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

        SimdDInt32(std::int32_t i) : simdInternal_(_mm512_set1_epi32(i)) {}

        // Internal utility constructor to simplify return statements
        SimdDInt32(__m512i simd) : simdInternal_(simd) {}

        __m512i  simdInternal_;
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
simdLoad(const double *m, SimdDoubleTag = {})
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
simdLoadU(const double *m, SimdDoubleTag = {})
{
    return {
               _mm512_loadunpackhi_pd(_mm512_loadunpacklo_pd(_mm512_undefined_pd(), m), m+8)
    };
}

static inline void gmx_simdcall
storeU(double *m, SimdDouble a)
{
    _mm512_packstorelo_pd(m, a.simdInternal_);
    _mm512_packstorehi_pd(m+8, a.simdInternal_);

}

static inline SimdDouble gmx_simdcall
setZeroD()
{
    return {
               _mm512_setzero_pd()
    };
}

static inline SimdDInt32 gmx_simdcall
simdLoad(const std::int32_t * m, SimdDInt32Tag)
{
    assert(std::size_t(m) % 32 == 0);
    return {
               _mm512_extload_epi64(m, _MM_UPCONV_EPI64_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE)
    };
}

static inline void gmx_simdcall
store(std::int32_t * m, SimdDInt32 a)
{
    assert(std::size_t(m) % 32 == 0);
    _mm512_mask_packstorelo_epi32(m, _mm512_int2mask(0x00FF), a.simdInternal_);
}

static inline SimdDInt32 gmx_simdcall
simdLoadU(const std::int32_t *m, SimdDInt32Tag)
{
    return {
               _mm512_mask_loadunpackhi_epi32(_mm512_mask_loadunpacklo_epi32(_mm512_undefined_epi32(), _mm512_int2mask(0x00FF), m),
                                              _mm512_int2mask(0x00FF), m+16)
    };
}

static inline void gmx_simdcall
storeU(std::int32_t * m, SimdDInt32 a)
{
    _mm512_mask_packstorelo_epi32(m, _mm512_int2mask(0x00FF), a.simdInternal_);
    _mm512_mask_packstorehi_epi32(m+16, _mm512_int2mask(0x00FF), a.simdInternal_);
}

static inline SimdDInt32 gmx_simdcall
setZeroDI()
{
    return {
               _mm512_setzero_epi32()
    };
}

template<int index>
static inline std::int32_t gmx_simdcall
extract(SimdDInt32 a)
{
    int r;
    _mm512_mask_packstorelo_epi32(&r, _mm512_mask2int(1<<index), a.simdInternal_);
    return r;
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
               _mm512_addn_pd(x.simdInternal_, _mm512_setzero_pd())
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

static inline SimdDouble gmx_simdcall
rsqrt(SimdDouble x)
{
    return {
               _mm512_cvtpslo_pd(_mm512_rsqrt23_ps(_mm512_cvtpd_pslo(x.simdInternal_)))
    };
}

static inline SimdDouble gmx_simdcall
rcp(SimdDouble x)
{
    return {
               _mm512_cvtpslo_pd(_mm512_rcp23_ps(_mm512_cvtpd_pslo(x.simdInternal_)))
    };
}

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
               _mm512_mask_mul_pd(_mm512_setzero_pd(), m.simdInternal_, a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
maskzFma(SimdDouble a, SimdDouble b, SimdDouble c, SimdDBool m)
{
    return {
               _mm512_mask_mov_pd(_mm512_setzero_pd(), m.simdInternal_, _mm512_fmadd_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_))
    };
}

static inline SimdDouble gmx_simdcall
maskzRsqrt(SimdDouble x, SimdDBool m)
{
    return {
               _mm512_cvtpslo_pd(_mm512_mask_rsqrt23_ps(_mm512_setzero_ps(), m.simdInternal_, _mm512_cvtpd_pslo(x.simdInternal_)))
    };
}

static inline SimdDouble gmx_simdcall
maskzRcp(SimdDouble x, SimdDBool m)
{
    return {
               _mm512_cvtpslo_pd(_mm512_mask_rcp23_ps(_mm512_setzero_ps(), m.simdInternal_, _mm512_cvtpd_pslo(x.simdInternal_)))
    };
}

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
               _mm512_gmax_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
min(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_gmin_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
round(SimdDouble x)
{
    return {
               _mm512_roundfxpnt_adjust_pd(x.simdInternal_, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
    };
}

static inline SimdDouble gmx_simdcall
trunc(SimdDouble x)
{
    return {
               _mm512_roundfxpnt_adjust_pd(x.simdInternal_, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
    };
}

static inline SimdDouble
frexp(SimdDouble value, SimdDInt32 * exponent)
{
    __m512d rExponent = _mm512_getexp_pd(value.simdInternal_);
    __m512i iExponent = _mm512_cvtfxpnt_roundpd_epi32lo(rExponent, _MM_FROUND_TO_NEAREST_INT);

    exponent->simdInternal_ = _mm512_add_epi32(iExponent, _mm512_set1_epi32(1));

    return {
               _mm512_getmant_pd(value.simdInternal_, _MM_MANT_NORM_p5_1, _MM_MANT_SIGN_src)
    };
}

template <MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble
ldexp(SimdDouble value, SimdDInt32 exponent)
{
    const __m512i exponentBias = _mm512_set1_epi32(1023);
    __m512i       iExponent    = _mm512_add_epi32(exponent.simdInternal_, exponentBias);

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponent = _mm512_max_epi32(iExponent, _mm512_setzero_epi32());
    }

    iExponent = _mm512_permutevar_epi32(_mm512_set_epi32(7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 0, 0), iExponent);
    iExponent = _mm512_mask_slli_epi32(_mm512_setzero_epi32(), _mm512_int2mask(0xAAAA), iExponent, 20);
    return _mm512_mul_pd(_mm512_castsi512_pd(iExponent), value.simdInternal_);
}

static inline double gmx_simdcall
reduce(SimdDouble a)
{
    return _mm512_reduce_add_pd(a.simdInternal_);
}

// Picky, picky, picky:
// icc-16 complains about "Illegal value of immediate argument to intrinsic"
// unless we use
// 1) Ordered-quiet for ==
// 2) Unordered-quiet for !=
// 3) Ordered-signaling for < and <=

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
               _mm512_cmp_pd_mask(a.simdInternal_, b.simdInternal_, _CMP_NEQ_UQ)
    };
}

static inline SimdDBool gmx_simdcall
operator<(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_cmp_pd_mask(a.simdInternal_, b.simdInternal_, _CMP_LT_OS)
    };
}

static inline SimdDBool gmx_simdcall
operator<=(SimdDouble a, SimdDouble b)
{
    return {
               _mm512_cmp_pd_mask(a.simdInternal_, b.simdInternal_, _CMP_LE_OS)
    };
}

static inline SimdDBool gmx_simdcall
testBits(SimdDouble a)
{
    // This is a bit problematic since Knight's corner does not have any 64-bit integer comparisons,
    // and we cannot use floating-point since values with just a single bit set can evaluate to 0.0.
    // Instead, we do it as
    // 1) Do a logical or of the high/low 32 bits
    // 2) Do a permute so we have the low 32 bits of each value in the low 8 32-bit elements
    // 3) Do an integer comparison, and cast so we just keep the low 8 bits of the mask.
    //
    // By default we will use integers for the masks in the nonbonded kernels, so this shouldn't
    // have any significant performance drawbacks.

    __m512i ia = _mm512_castpd_si512(a.simdInternal_);

    ia = _mm512_or_epi32(ia, _mm512_swizzle_epi32(ia, _MM_SWIZ_REG_CDAB));
    ia = _mm512_permutevar_epi32( _mm512_set_epi32(15, 13, 11, 9, 7, 5, 3, 1, 14, 12, 10, 8, 6, 4, 2, 0), ia);

    return {
               static_cast<__mmask8>(_mm512_cmp_epi32_mask(ia, _mm512_setzero_si512(), _MM_CMPINT_NE))
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
    return _mm512_mask2int(a.simdInternal_) != 0;
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
operator&(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm512_and_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
andNot(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm512_andnot_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator|(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm512_or_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator^(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm512_xor_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator+(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm512_add_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator-(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm512_sub_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
operator*(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm512_mullo_epi32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDIBool gmx_simdcall
operator==(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm512_cmp_epi32_mask(a.simdInternal_, b.simdInternal_, _MM_CMPINT_EQ)
    };
}

static inline SimdDIBool gmx_simdcall
testBits(SimdDInt32 a)
{
    return {
               _mm512_cmp_epi32_mask(a.simdInternal_, _mm512_setzero_si512(), _MM_CMPINT_NE)
    };
}

static inline SimdDIBool gmx_simdcall
operator<(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm512_cmp_epi32_mask(a.simdInternal_, b.simdInternal_, _MM_CMPINT_LT)
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
    return ( _mm512_mask2int(a.simdInternal_) & 0xFF) != 0;
}

static inline SimdDInt32 gmx_simdcall
selectByMask(SimdDInt32 a, SimdDIBool m)
{
    return {
               _mm512_mask_mov_epi32(_mm512_setzero_epi32(), m.simdInternal_, a.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
selectByNotMask(SimdDInt32 a, SimdDIBool m)
{
    return {
               _mm512_mask_mov_epi32(a.simdInternal_, m.simdInternal_, _mm512_setzero_epi32())
    };
}

static inline SimdDInt32 gmx_simdcall
blend(SimdDInt32 a, SimdDInt32 b, SimdDIBool sel)
{
    return {
               _mm512_mask_blend_epi32(sel.simdInternal_, a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
cvtR2I(SimdDouble a)
{
    return {
               _mm512_cvtfxpnt_roundpd_epi32lo(a.simdInternal_, _MM_FROUND_TO_NEAREST_INT)
    };
}

static inline SimdDInt32 gmx_simdcall
cvttR2I(SimdDouble a)
{
    return {
               _mm512_cvtfxpnt_roundpd_epi32lo(a.simdInternal_, _MM_FROUND_TO_ZERO)
    };
}

static inline SimdDouble gmx_simdcall
cvtI2R(SimdDInt32 a)
{
    return {
               _mm512_cvtepi32lo_pd(a.simdInternal_)
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
    __m512i i1 = _mm512_permute4f128_epi32(_mm512_castps_si512(f.simdInternal_), _MM_PERM_DCDC);

    *d0 = _mm512_cvtpslo_pd(f.simdInternal_);
    *d1 = _mm512_cvtpslo_pd(_mm512_castsi512_ps(i1));
}

static inline SimdFloat gmx_simdcall
cvtDD2F(SimdDouble d0, SimdDouble d1)
{
    __m512 f0 = _mm512_cvtpd_pslo(d0.simdInternal_);
    __m512 f1 = _mm512_cvtpd_pslo(d1.simdInternal_);
    return {
               _mm512_mask_permute4f128_ps(f0, _mm512_int2mask(0xFF00), f1, _MM_PERM_BABA)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_MIC_SIMD_DOUBLE_H
