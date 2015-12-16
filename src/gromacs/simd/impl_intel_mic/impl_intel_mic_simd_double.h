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

#ifndef GMX_SIMD_IMPL_INTEL_MIC_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_INTEL_MIC_SIMD_DOUBLE_H

#include "config.h"

#include <assert.h>

#include <cmath>
#include <cstdint>

#include <immintrin.h>

#include "impl_intel_mic_common.h"

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdDouble          __m512d
#define simdLoadD            _mm512_load_pd
#define simdLoad1D(m)        _mm512_extload_pd(m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE)
#define simdSet1D            _mm512_set1_pd
#define simdStoreD           _mm512_store_pd
#define simdLoadUD           simdLoadUD_mic
#define simdStoreUD          simdStoreUD_mic
#define simdSetZeroD         _mm512_setzero_pd
#define simdAddD             _mm512_add_pd
#define simdSubD             _mm512_sub_pd
#define simdMulD             _mm512_mul_pd
#define simdFmaddD           _mm512_fmadd_pd
#define simdFmsubD           _mm512_fmsub_pd
#define simdFnmaddD          _mm512_fnmadd_pd
#define simdFnmsubD          _mm512_fnmsub_pd
#define simdAndD(a, b)       _mm512_castsi512_pd(_mm512_and_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define simdAndNotD(a, b)    _mm512_castsi512_pd(_mm512_andnot_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define simdOrD(a, b)        _mm512_castsi512_pd(_mm512_or_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define simdXorD(a, b)       _mm512_castsi512_pd(_mm512_xor_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define simdRsqrtD(x)        _mm512_cvtpslo_pd(_mm512_rsqrt23_ps(_mm512_cvtpd_pslo(x)))
#define simdRcpD(x)          _mm512_cvtpslo_pd(_mm512_rcp23_ps(_mm512_cvtpd_pslo(x)))
#define simdAbsD(x)         simdAndNotD(_mm512_set1_pd(GMX_DOUBLE_NEGZERO), x)
#define simdNegD(x)         _mm512_addn_pd(x, _mm512_setzero_pd())
#define simdMaxD             _mm512_gmax_pd
#define simdMinD             _mm512_gmin_pd
#define simdRoundD(a)        _mm512_roundfxpnt_adjust_pd(a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define simdTruncD(a)        _mm512_roundfxpnt_adjust_pd(a, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define simdFractionD(x)     _mm512_sub_pd(x, simdTruncD(x))
#define simdGetExponentD(x) _mm512_getexp_pd(x)
#define simdGetMantissaD(x) _mm512_getmant_pd(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero)
#define simdSetExponentD(x) simdSetExponentD_mic(x)
/* integer datatype corresponding to float: SimdFInt32
   Doesn't use mask other than where required. No side effect expected for operating on the (unused) upper 8.
 */
#define SimdDInt32          __m512i
#define simdLoadDI(m)        _mm512_extload_epi64(m, _MM_UPCONV_EPI64_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE)
#define simdSet1DI           _mm512_set1_epi32
#define simdStoreDI          simdStoreDI_mic
#define simdLoadUDI          simdLoadUDI_mic
#define simdStoreUDI         simdStoreUDI_mic
#define simdExtractDI        simdExtractDI_mic
#define simdSetZeroDI        _mm512_setzero_epi32
#define simdCvtD2I(a)        _mm512_cvtfxpnt_roundpd_epi32lo(a, _MM_FROUND_TO_NEAREST_INT)
#define simdCvttD2I(a)       _mm512_cvtfxpnt_roundpd_epi32lo(a, _MM_FROUND_TO_ZERO)
#define simdCvtI2D           _mm512_cvtepi32lo_pd
/* Integer logical ops on SimdFInt32 */
#define simdSlliDI           _mm512_slli_epi32
#define simdSrliDI           _mm512_srli_epi32
#define simdAndDI            _mm512_and_epi32
#define simdAndNotDI         _mm512_andnot_epi32
#define simdOrDI             _mm512_or_epi32
#define simdXorDI            _mm512_xor_epi32
/* Integer arithmetic ops on SimdFInt32 */
#define simdAddDI            _mm512_add_epi32
#define simdSubDI            _mm512_sub_epi32
#define simdMulDI            _mm512_mullo_epi32
/* Boolean & comparison operations on SimdFloat */
#define SimdDBool           __mmask8
#define simdCmpEqD(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_EQ_OQ)
#define simdCmpLtD(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_LT_OS)
#define simdCmpLeD(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_LE_OS)
#define simdAndDB            _mm512_kand
#define simdOrDB             _mm512_kor
#define simdAnyTrueDB(x)     _mm512_mask2int(x)
#define simdMaskD(a, sel)    _mm512_mask_mov_pd(_mm512_setzero_pd(), sel, a)
#define simdMaskNotD(a, sel) _mm512_mask_mov_pd(_mm512_setzero_pd(), _mm512_knot(sel), a)
#define simdBlendD(a, b, sel)    _mm512_mask_blend_pd(sel, a, b)
#define simdReduceD(a)       _mm512_reduce_add_pd(a)
/* Boolean & comparison operations on SimdFInt32 */
#define SimdDIBool          __mmask16
#define simdCmpEqDI(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_EQ)
#define simdCmpLtDI(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_LT)
#define simdAndDIB           _mm512_kand
#define simdOrDIB            _mm512_kor
#define simdAnyTrueDIB(x)    (_mm512_mask2int(x)&0xFF)
#define simdMaskDI(a, sel)    _mm512_mask_mov_epi32(_mm512_setzero_epi32(), sel, a)
#define simdMaskNotDI(a, sel) _mm512_mask_mov_epi32(_mm512_setzero_epi32(), _mm512_knot(sel), a)
#define simdBlendDI(a, b, sel)    _mm512_mask_blend_epi32(sel, a, b)
/* Conversions between booleans. Double & dint stuff is stored in low bits */
#define simdCvtDB2DIB(x)     (x)
#define simdCvtDIB2DB(x)     (x)

/* Float/double conversion */
#define simdCvtF2DD          simdCvtF2DD_mic
#define simdCvtDD2F          simdCvtDD2F_mic


#define PERM_LOW2HIGH _MM_PERM_BABA
#define PERM_HIGH2LOW _MM_PERM_DCDC

#define mask_loh _mm512_int2mask(0x00FF) /* would be better a constant - but can't initialize with a function call. */
#define mask_hih _mm512_int2mask(0xFF00)

/* load store double */
static inline __m512d gmx_simdcall
simdLoadUD_mic(const double * m)
{
    return _mm512_loadunpackhi_pd(_mm512_loadunpacklo_pd(_mm512_undefined_pd(), m), m+8);
}

static inline void gmx_simdcall
simdStoreUD_mic(double * m, __m512d s)
{
    _mm512_packstorelo_pd(m, s);
    _mm512_packstorehi_pd(m+8, s);
}

/* load store dint32 */
static inline void gmx_simdcall
simdStoreDI_mic(std::int32_t * m, __m512i s)
{
    assert((size_t)m%32 == 0);
    _mm512_mask_packstorelo_epi32(m, mask_loh, s);
}

static inline __m512i gmx_simdcall
simdLoadUDI_mic(const std::int32_t * m)
{
    return _mm512_mask_loadunpackhi_epi32(_mm512_mask_loadunpacklo_epi32(_mm512_undefined_epi32(), mask_loh, m), mask_loh, m+16);
}

static inline void gmx_simdcall
simdStoreUDI_mic(std::int32_t * m, __m512i s)
{
    _mm512_mask_packstorelo_epi32(m, mask_loh, s);
    _mm512_mask_packstorehi_epi32(m+16, mask_loh, s);
}

static inline std::int32_t gmx_simdcall
simdExtractDI_mic(SimdDInt32 a, int index)
{
    int r;
    _mm512_mask_packstorelo_epi32(&r, _mm512_mask2int(1<<index), a);
    return r;
}

/* This is likely faster than the built in scale operation (lat 8, t-put 3)
 * since we only work on the integer part and use shifts. TODO: check. given that scale also only does integer
 */
static inline __m512d gmx_simdcall
simdSetExponentD_mic(__m512d a)
{
    const __m512i expbias      = _mm512_set1_epi32(1023);
    __m512i       iexp         = _mm512_cvtfxpnt_roundpd_epi32lo(a, _MM_FROUND_TO_NEAREST_INT);
    iexp = _mm512_permutevar_epi32(_mm512_set_epi32(7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 0, 0), iexp);
    iexp = _mm512_mask_slli_epi32(_mm512_setzero_epi32(), _mm512_int2mask(0xAAAA), _mm512_add_epi32(iexp, expbias), 20);
    return _mm512_castsi512_pd(iexp);
}

static inline void gmx_simdcall
simdCvtF2DD_mic(__m512 f, __m512d * d0, __m512d * d1)
{
    __m512i i1 = _mm512_permute4f128_epi32(_mm512_castps_si512(f), PERM_HIGH2LOW);

    *d0 = _mm512_cvtpslo_pd(f);
    *d1 = _mm512_cvtpslo_pd(_mm512_castsi512_ps(i1));
}

static inline __m512 gmx_simdcall
simdCvtDD2F_mic(__m512d d0, __m512d d1)
{
    __m512 f0 = _mm512_cvtpd_pslo(d0);
    __m512 f1 = _mm512_cvtpd_pslo(d1);
    return _mm512_mask_permute4f128_ps(f0, mask_hih, f1, PERM_LOW2HIGH);
}

#endif /* GMX_SIMD_IMPL_INTEL_MIC_SIMD_DOUBLE_H */
