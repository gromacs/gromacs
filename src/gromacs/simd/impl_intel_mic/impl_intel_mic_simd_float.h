/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_INTEL_MIC_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_INTEL_MIC_SIMD_FLOAT_H

#include "config.h"

#include <cmath>
#include <cstdint>

#include <immintrin.h>

#include "impl_intel_mic_common.h"

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdFloat           __m512
#define simdLoadF            _mm512_load_ps
#define simdLoad1F(m)        _mm512_extload_ps(m, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE)
#define simdSet1F            _mm512_set1_ps
#define simdStoreF           _mm512_store_ps
#define simdLoadUF           simdLoadUF_mic
#define simdStoreUF          simdStoreUF_mic
#define simdSetZeroF         _mm512_setzero_ps
#define simdAddF             _mm512_add_ps
#define simdSubF             _mm512_sub_ps
#define simdMulF             _mm512_mul_ps
#define simdFmaddF           _mm512_fmadd_ps
#define simdFmsubF           _mm512_fmsub_ps
#define simdFnmaddF          _mm512_fnmadd_ps
#define simdFnmsubF          _mm512_fnmsub_ps
#define simdAndF(a, b)        _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define simdAndNotF(a, b)     _mm512_castsi512_ps(_mm512_andnot_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define simdOrF(a, b)         _mm512_castsi512_ps(_mm512_or_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define simdXorF(a, b)        _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define simdRsqrtF           _mm512_rsqrt23_ps
#define simdRcpF             _mm512_rcp23_ps
#define simdAbsF(x)         simdAndNotF(_mm512_set1_ps(GMX_FLOAT_NEGZERO), x)
#define simdNegF(x)         _mm512_addn_ps(x, _mm512_setzero_ps())
#define simdMaxF             _mm512_gmax_ps
#define simdMinF             _mm512_gmin_ps
#define simdRoundF(x)        _mm512_round_ps(x, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define simdTruncF(x)        _mm512_round_ps(x, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define simdFractionF(x)     _mm512_sub_ps(x, simdTruncF(x))
#define simdGetExponentF(x) _mm512_getexp_ps(x)
#define simdGetMantissaF(x) _mm512_getmant_ps(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero)
#define simdSetExponentF(x) simdSetExponentF_mic(x)
/* integer datatype corresponding to float: SimdFInt32 */
#define SimdFInt32          __m512i
#define simdLoadFI           _mm512_load_epi32
#define simdSet1FI           _mm512_set1_epi32
#define simdStoreFI          _mm512_store_epi32
#define simdLoadUFI          simdLoadUFI_mic
#define simdStoreUFI         simdStoreUFI_mic
#define simdExtractFI        simdExtractFI_mic
#define simdSetZeroFI        _mm512_setzero_epi32
#define simdCvtF2I(a)        _mm512_cvtfxpnt_round_adjustps_epi32(a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define simdCvttF2I(a)       _mm512_cvtfxpnt_round_adjustps_epi32(a, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define simdCvtI2F(a)        _mm512_cvtfxpnt_round_adjustepi32_ps(a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
/* Integer logical ops on SimdFInt32 */
#define simdSlliFI           _mm512_slli_epi32
#define simdSrliFI           _mm512_srli_epi32
#define simdAndFI            _mm512_and_epi32
#define simdAndNotFI         _mm512_andnot_epi32
#define simdOrFI             _mm512_or_epi32
#define simdXorFI            _mm512_xor_epi32
/* Integer arithmetic ops on SimdFInt32 */
#define simdAddFI            _mm512_add_epi32
#define simdSubFI            _mm512_sub_epi32
#define simdMulFI            _mm512_mullo_epi32
/* Boolean & comparison operations on SimdFloat */
#define SimdFBool           __mmask16
#define simdCmpEqF(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_EQ_OQ)
#define simdCmpLtF(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_LT_OS)
#define simdCmpLeF(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_LE_OS)
#define simdAndFB            _mm512_kand
#define simdAndNotFB(a, b)   _mm512_knot(_mm512_kor(a, b))
#define simdOrFB             _mm512_kor
#define simdAnyTrueFB        _mm512_mask2int
#define simdMaskF(a, sel)    _mm512_mask_mov_ps(_mm512_setzero_ps(), sel, a)
#define simdMaskNotF(a, sel) _mm512_mask_mov_ps(_mm512_setzero_ps(), _mm512_knot(sel), a)
#define simdBlendF(a, b, sel)    _mm512_mask_blend_ps(sel, a, b)
#define simdReduceF(a)       _mm512_reduce_add_ps(a)
/* Boolean & comparison operations on SimdFInt32 */
#define SimdFIBool          __mmask16
#define simdCmpEqFI(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_EQ)
#define simdCmpLtFI(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_LT)
#define simdAndFIB           _mm512_kand
#define simdOrFIB            _mm512_kor
#define simdAnyTrueFIB       _mm512_mask2int
#define simdMaskFI(a, sel)    _mm512_mask_mov_epi32(_mm512_setzero_epi32(), sel, a)
#define simdMaskNotFI(a, sel) _mm512_mask_mov_epi32(_mm512_setzero_epi32(), _mm512_knot(sel), a)
#define simdBlendFI(a, b, sel)    _mm512_mask_blend_epi32(sel, a, b)
/* Conversions between different booleans */
#define simdCvtFB2FIB(x)     (x)
#define simdCvtFIB2FB(x)     (x)

/* MIC provides full single precision of some neat functions: */
/* 1/sqrt(x) and 1/x work fine in simd_math.h, and won't use extra iterations */
#define simdExp2F            simdExp2F_mic
#define simdExpF             simdExpF_mic
#define simdLogF             simdLogF_mic

/* load store float */
static inline __m512 gmx_simdcall
simdLoadUF_mic(const float * m)
{
    return _mm512_loadunpackhi_ps(_mm512_loadunpacklo_ps(_mm512_undefined_ps(), m), m+16);
}

static inline void gmx_simdcall
simdStoreUF_mic(float * m, __m512 s)
{
    _mm512_packstorelo_ps(m, s);
    _mm512_packstorehi_ps(m+16, s);
}

/* load store fint32 */
static inline __m512i gmx_simdcall
simdLoadUFI_mic(const std::int32_t * m)
{
    return _mm512_loadunpackhi_epi32(_mm512_loadunpacklo_epi32(_mm512_undefined_epi32(), m), m+16);
}

static inline void gmx_simdcall
simdStoreUFI_mic(std::int32_t * m, __m512i s)
{
    _mm512_packstorelo_epi32(m, s);
    _mm512_packstorehi_epi32(m+16, s);
}

/* extract */
static inline std::int32_t gmx_simdcall
simdExtractFI_mic(SimdFInt32 a, int index)
{
    int r;
    _mm512_mask_packstorelo_epi32(&r, _mm512_mask2int(1<<index), a);
    return r;
}

/* This is likely faster than the built in scale operation (lat 8, t-put 3)
 * since we only work on the integer part and use shifts. TODO: check. given that scale also only does integer
 */
static inline __m512 gmx_simdcall
simdSetExponentF_mic(__m512 a)
{
    __m512i       iexp         = simdCvtF2I(a);

    const __m512i expbias      = _mm512_set1_epi32(127);
    iexp = _mm512_slli_epi32(_mm512_add_epi32(iexp, expbias), 23);
    return _mm512_castsi512_ps(iexp);

    /* scale alternative:
       return _mm512_scale_ps(_mm512_set1_ps(1), iexp);
     */
}

static inline __m512 gmx_simdcall
simdExp2F_mic(__m512 x)
{
    return _mm512_exp223_ps(_mm512_cvtfxpnt_round_adjustps_epi32(x, _MM_ROUND_MODE_NEAREST, _MM_EXPADJ_24));
}

static inline __m512 gmx_simdcall
simdExpF_mic(__m512 x)
{
    const SimdFloat         argscale    = simdSet1F(1.44269504088896341f);
    const SimdFloat         invargscale = simdSet1F(-0.69314718055994528623f);
    __m512                  xscaled     = _mm512_mul_ps(x, argscale);
    __m512                  r           = simdExp2F_mic(xscaled);

    /* simdExp2F_mic() provides 23 bits of accuracy, but we ruin some of that
     * with the argument scaling due to single-precision rounding, where the
     * rounding error is amplified exponentially. To correct this, we find the
     * difference between the scaled argument and the true one (extended precision
     * arithmetics does not appear to be necessary to fulfill our accuracy requirements)
     * and then multiply by the exponent of this correction since exp(a+b)=exp(a)*exp(b).
     * Note that this only adds two instructions (and maybe some constant loads).
     */
    x         = simdFmaddF(invargscale, xscaled, x);
    /* x will now be a _very_ small number, so approximate exp(x)=1+x.
     * We should thus apply the correction as r'=r*(1+x)=r+r*x
     */
    r         = simdFmaddF(r, x, r);
    return r;
}

static inline __m512 gmx_simdcall
simdLogF_mic(__m512 x)
{
    return _mm512_mul_ps(_mm512_set1_ps(0.693147180559945286226764), _mm512_log2ae23_ps(x));
}

#endif /* GMX_SIMD_IMPL_INTEL_MIC_SIMD_FLOAT_H */
