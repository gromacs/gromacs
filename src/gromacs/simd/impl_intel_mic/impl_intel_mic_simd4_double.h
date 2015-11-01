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

#ifndef GMX_SIMD_IMPL_INTEL_MIC_SIMD4_DOUBLE_H
#define GMX_SIMD_IMPL_INTEL_MIC_SIMD4_DOUBLE_H

#include "config.h"

#include <assert.h>
#include <math.h>

#include <immintrin.h>

#include "impl_intel_mic_common.h"
#include "impl_intel_mic_simd_double.h"

/****************************************************
 *      DOUBLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
#define Simd4Double          __m512d
#define simd4Mask              _mm512_int2mask(0xF)
#define simd4LoadD(m)         _mm512_mask_extload_pd(_mm512_undefined_pd(), simd4Mask, m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE)
#define simd4Load1D(m)        _mm512_mask_extload_pd(_mm512_undefined_pd(), simd4Mask, m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE)
#define simd4Set1D            _mm512_set1_pd
#define simd4StoreD           simd4StoreD_mic
#define simd4LoadUD           simd4LoadUD_mic
#define simd4StoreUD          simd4StoreUD_mic
#define simd4SetZeroD         _mm512_setzero_pd
#define simd4AddD(a, b)       _mm512_mask_add_pd(_mm512_undefined_pd(), simd4Mask, a, b)
#define simd4SubD(a, b)       _mm512_mask_sub_pd(_mm512_undefined_pd(), simd4Mask, a, b)
#define simd4MulD(a, b)       _mm512_mask_mul_pd(_mm512_undefined_pd(), simd4Mask, a, b)
#define simd4FmaddD(a, b, c)  _mm512_mask_fmadd_pd(a, simd4Mask, b, c)
#define simd4FmsubD(a, b, c)  _mm512_mask_fmsub_pd(a, simd4Mask, b, c)
#define simd4FnmaddD(a, b, c) _mm512_mask_fnmadd_pd(a, simd4Mask, b, c)
#define simd4FnmsubD(a, b, c) _mm512_mask_fnmsub_pd(a, simd4Mask, b, c)
#define simd4AndD(a, b)       _mm512_castsi512_pd(_mm512_mask_and_epi32(_mm512_undefined_epi32(), mask_loh, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define simd4AndNotD(a, b)    _mm512_castsi512_pd(_mm512_mask_andnot_epi32(_mm512_undefined_epi32(), mask_loh, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define simd4OrD(a, b)        _mm512_castsi512_pd(_mm512_mask_or_epi32(_mm512_undefined_epi32(), mask_loh, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define simd4XorD(a, b)       _mm512_castsi512_pd(_mm512_mask_xor_epi32(_mm512_undefined_epi32(), mask_loh, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define simd4RsqrtD(a)        _mm512_mask_cvtpslo_pd(_mm512_undefined_pd(), simd4Mask, _mm512_mask_rsqrt23_ps(_mm512_undefined_ps(), simd4Mask, _mm512_mask_cvtpd_pslo(_mm512_undefined_ps(), simd4Mask, x)))
#define simd4AbsD(x)         simd4AndNotD(_mm512_set1_pd(GMX_DOUBLE_NEGZERO), x)
#define simd4NegD(x)         _mm512_mask_addn_pd(_mm512_undefined_pd(), simd4Mask, x, _mm512_setzero_pd())
#define simd4MaxD(a, b)       _mm512_mask_gmax_pd(_mm512_undefined_pd(), simd4Mask, a, b)
#define simd4MinD(a, b)       _mm512_mask_gmin_pd(_mm512_undefined_pd(), simd4Mask, a, b)
#define simd4RoundD(a)        _mm512_mask_roundfxpnt_adjust_pd(_mm512_undefined_pd(), simd4Mask, a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define simd4TruncD(a)        _mm512_mask_roundfxpnt_adjust_pd(_mm512_undefined_pd(), simd4Mask, a, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define simd4DotProductD(a, b) _mm512_mask_reduce_add_pd(_mm512_int2mask(7), _mm512_mask_mul_pd(_mm512_undefined_pd(), _mm512_int2mask(7), a, b))
#define Simd4DBool           __mmask16
#define simd4CmpEqD(a, b)     _mm512_mask_cmp_pd_mask(simd4Mask, a, b, _CMP_EQ_OQ)
#define simd4CmpLtD(a, b)     _mm512_mask_cmp_pd_mask(simd4Mask, a, b, _CMP_LT_OS)
#define simd4CmpLeD(a, b)     _mm512_mask_cmp_pd_mask(simd4Mask, a, b, _CMP_LE_OS)
#define simd4AndDB            _mm512_kand
#define simd4OrDB             _mm512_kor
#define simd4AnyTrueDB(x)     (_mm512_mask2int(x)&0xF)
#define simd4MaskD(a, sel)    _mm512_mask_mov_pd(_mm512_setzero_pd(), sel, a)
#define simd4MaskNotD(a, sel) _mm512_mask_mov_pd(_mm512_setzero_pd(), _mm512_knot(sel), a)
#define simd4BlendD(a, b, sel)    _mm512_mask_blend_pd(sel, a, b)
#define simd4ReduceD(x)       _mm512_mask_reduce_add_pd(_mm512_int2mask(0xF), x)

/* Implementation helpers */
static inline void gmx_simdcall
simd4StoreD_mic(double * m, __m512d s)
{
    assert((size_t)m%32 == 0);
    _mm512_mask_packstorelo_pd(m, simd4Mask, s);
}

static inline __m512d gmx_simdcall
simd4LoadUD_mic(const double * m)
{
    return _mm512_mask_loadunpackhi_pd(_mm512_mask_loadunpacklo_pd(_mm512_undefined_pd(), simd4Mask, m), simd4Mask, m+8);
}

static inline void gmx_simdcall
simd4StoreUD_mic(double * m, __m512d s)
{
    _mm512_mask_packstorelo_pd(m, simd4Mask, s);
    _mm512_mask_packstorehi_pd(m+8, simd4Mask, s);
}

#endif /* GMX_SIMD_IMPL_INTEL_MIC_SIMD4_DOUBLE_H */
