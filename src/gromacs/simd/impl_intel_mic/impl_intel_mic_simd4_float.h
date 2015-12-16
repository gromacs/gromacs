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

#ifndef GMX_SIMD_IMPL_INTEL_MIC_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_INTEL_MIC_SIMD4_FLOAT_H

#include "config.h"

#include <assert.h>
#include <math.h>

#include <immintrin.h>

#include "impl_intel_mic_common.h"
#include "impl_intel_mic_simd_float.h"

/****************************************************
 *      SINGLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
/* Load and store are guranteed to only access the 4 floats. All arithmetic operations
   only operate on the 4 elements (to avoid floating excpetions). But other operations
   are not gurateed to not modify the other 12 elements. E.g. setzero or blendzero
   set the upper 12 to zero. */
#define Simd4Float           __m512
#define simd4Mask              _mm512_int2mask(0xF)
#define simd4LoadF(m)         _mm512_mask_extload_ps(_mm512_undefined_ps(), simd4Mask, m, _MM_UPCONV_PS_NONE, _MM_BROADCAST_4X16, _MM_HINT_NONE)
#define simd4Load1F(m)        _mm512_mask_extload_ps(_mm512_undefined_ps(), simd4Mask, m, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE)
#define simd4Set1F            _mm512_set1_ps
#define simd4StoreF           simd4StoreF_mic
#define simd4LoadUF           simd4LoadUF_mic
#define simd4StoreUF          simd4StoreUF_mic
#define simd4SetZeroF         _mm512_setzero_ps
#define simd4AddF(a, b)       _mm512_mask_add_ps(_mm512_undefined_ps(), simd4Mask, a, b)
#define simd4SubF(a, b)       _mm512_mask_sub_ps(_mm512_undefined_ps(), simd4Mask, a, b)
#define simd4MulF(a, b)       _mm512_mask_mul_ps(_mm512_undefined_ps(), simd4Mask, a, b)
#define simd4FmaddF(a, b, c)  _mm512_mask_fmadd_ps(a, simd4Mask, b, c)
#define simd4FmsubF(a, b, c)  _mm512_mask_fmsub_ps(a, simd4Mask, b, c)
#define simd4FnmaddF(a, b, c) _mm512_mask_fnmadd_ps(a, simd4Mask, b, c)
#define simd4FnmsubF(a, b, c) _mm512_mask_fnmsub_ps(a, simd4Mask, b, c)
#define simd4AndF(a, b)       _mm512_castsi512_ps(_mm512_mask_and_epi32(_mm512_undefined_epi32(), simd4Mask, _mm512_castps_si512(a), _mm512_castps_si512(b)))
#define simd4AndNotF(a, b)    _mm512_castsi512_ps(_mm512_mask_andnot_epi32(_mm512_undefined_epi32(), simd4Mask, _mm512_castps_si512(a), _mm512_castps_si512(b)))
#define simd4OrF(a, b)        _mm512_castsi512_ps(_mm512_mask_or_epi32(_mm512_undefined_epi32(), simd4Mask, _mm512_castps_si512(a), _mm512_castps_si512(b)))
#define simd4XorF(a, b)       _mm512_castsi512_ps(_mm512_mask_xor_epi32(_mm512_undefined_epi32(), simd4Mask, _mm512_castps_si512(a), _mm512_castps_si512(b)))
#define simd4RsqrtF(a)        _mm512_mask_rsqrt23_ps(_mm512_undefined_ps(), simd4Mask, a)
#define simd4AbsF(x)         simd4AndNotF(_mm512_set1_ps(GMX_FLOAT_NEGZERO), x)
#define simd4NegF(x)         _mm512_mask_addn_ps(_mm512_undefined_ps(), simd4Mask, x, _mm512_setzero_ps())
#define simd4MaxF(a, b)       _mm512_mask_gmax_ps(_mm512_undefined_ps(), simd4Mask, a, b)
#define simd4MinF(a, b)       _mm512_mask_gmin_ps(_mm512_undefined_ps(), simd4Mask, a, b)
#define simd4RoundF(x)        _mm512_mask_round_ps(_mm512_undefined_ps(), simd4Mask, x, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define simd4TruncF(x)        _mm512_mask_round_ps(_mm512_undefined_ps(), simd4Mask, x, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define simd4DotProductF(a, b) _mm512_mask_reduce_add_ps(_mm512_int2mask(7), _mm512_mask_mul_ps(_mm512_undefined_ps(), _mm512_int2mask(7), a, b))
#define Simd4FBool           __mmask16
#define simd4CmpEqF(a, b)     _mm512_mask_cmp_ps_mask(simd4Mask, a, b, _CMP_EQ_OQ)
#define simd4CmpLtF(a, b)     _mm512_mask_cmp_ps_mask(simd4Mask, a, b, _CMP_LT_OS)
#define simd4CmpLeF(a, b)     _mm512_mask_cmp_ps_mask(simd4Mask, a, b, _CMP_LE_OS)
#define simd4AndFB            _mm512_kand
#define simd4OrFB             _mm512_kor
#define simd4AnyTrueFB(x)     (_mm512_mask2int(x)&0xF)
#define simd4MaskF(a, sel)    _mm512_mask_mov_ps(_mm512_setzero_ps(), sel, a)
#define simd4MaskNotF(a, sel) _mm512_mask_mov_ps(_mm512_setzero_ps(), _mm512_knot(sel), a)
#define simd4BlendF(a, b, sel)    _mm512_mask_blend_ps(sel, a, b)
#define simd4ReduceF(x)       _mm512_mask_reduce_add_ps(_mm512_int2mask(0xF), x)

/* Implementation helpers */

/* load store simd4 */
static inline void gmx_simdcall
simd4StoreF_mic(float * m, __m512 s)
{
    assert((size_t)m%16 == 0);
    _mm512_mask_packstorelo_ps(m, simd4Mask, s);
}

static inline __m512 gmx_simdcall
simd4LoadUF_mic(const float * m)
{
    return _mm512_mask_loadunpackhi_ps(_mm512_mask_loadunpacklo_ps(_mm512_undefined_ps(), simd4Mask, m), simd4Mask, m+16);
}

static inline void gmx_simdcall
simd4StoreUF_mic(float * m, __m512 s)
{
    _mm512_mask_packstorelo_ps(m, simd4Mask, s);
    _mm512_mask_packstorehi_ps(m+16, simd4Mask, s);
}

#endif /* GMX_SIMD_IMPL_INTEL_MIC_SIMD4_FLOAT_H */
