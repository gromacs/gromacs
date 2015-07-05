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
#define gmx_simd4_float_t           __m512
#define gmx_simd4_mask              _mm512_int2mask(0xF)
#define gmx_simd4_load_f(m)         _mm512_mask_extload_ps(_mm512_undefined_ps(), gmx_simd4_mask, m, _MM_UPCONV_PS_NONE, _MM_BROADCAST_4X16, _MM_HINT_NONE)
#define gmx_simd4_load1_f(m)        _mm512_mask_extload_ps(_mm512_undefined_ps(), gmx_simd4_mask, m, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE)
#define gmx_simd4_set1_f            _mm512_set1_ps
#define gmx_simd4_store_f           gmx_simd4_store_f_mic
#define gmx_simd4_loadu_f           gmx_simd4_loadu_f_mic
#define gmx_simd4_storeu_f          gmx_simd4_storeu_f_mic
#define gmx_simd4_setzero_f         _mm512_setzero_ps
#define gmx_simd4_add_f(a, b)       _mm512_mask_add_ps(_mm512_undefined_ps(), gmx_simd4_mask, a, b)
#define gmx_simd4_sub_f(a, b)       _mm512_mask_sub_ps(_mm512_undefined_ps(), gmx_simd4_mask, a, b)
#define gmx_simd4_mul_f(a, b)       _mm512_mask_mul_ps(_mm512_undefined_ps(), gmx_simd4_mask, a, b)
#define gmx_simd4_fmadd_f(a, b, c)  _mm512_mask_fmadd_ps(a, gmx_simd4_mask, b, c)
#define gmx_simd4_fmsub_f(a, b, c)  _mm512_mask_fmsub_ps(a, gmx_simd4_mask, b, c)
#define gmx_simd4_fnmadd_f(a, b, c) _mm512_mask_fnmadd_ps(a, gmx_simd4_mask, b, c)
#define gmx_simd4_fnmsub_f(a, b, c) _mm512_mask_fnmsub_ps(a, gmx_simd4_mask, b, c)
#define gmx_simd4_and_f(a, b)       _mm512_castsi512_ps(_mm512_mask_and_epi32(_mm512_undefined_epi32(), gmx_simd4_mask, _mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd4_andnot_f(a, b)    _mm512_castsi512_ps(_mm512_mask_andnot_epi32(_mm512_undefined_epi32(), gmx_simd4_mask, _mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd4_or_f(a, b)        _mm512_castsi512_ps(_mm512_mask_or_epi32(_mm512_undefined_epi32(), gmx_simd4_mask, _mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd4_xor_f(a, b)       _mm512_castsi512_ps(_mm512_mask_xor_epi32(_mm512_undefined_epi32(), gmx_simd4_mask, _mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd4_rsqrt_f(a)        _mm512_mask_rsqrt23_ps(_mm512_undefined_ps(), gmx_simd4_mask, a)
#define gmx_simd4_fabs_f(x)         gmx_simd4_andnot_f(_mm512_set1_ps(GMX_FLOAT_NEGZERO), x)
#define gmx_simd4_fneg_f(x)         _mm512_mask_addn_ps(_mm512_undefined_ps(), gmx_simd4_mask, x, _mm512_setzero_ps())
#define gmx_simd4_max_f(a, b)       _mm512_mask_gmax_ps(_mm512_undefined_ps(), gmx_simd4_mask, a, b)
#define gmx_simd4_min_f(a, b)       _mm512_mask_gmin_ps(_mm512_undefined_ps(), gmx_simd4_mask, a, b)
#define gmx_simd4_round_f(x)        _mm512_mask_round_ps(_mm512_undefined_ps(), gmx_simd4_mask, x, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define gmx_simd4_trunc_f(x)        _mm512_mask_round_ps(_mm512_undefined_ps(), gmx_simd4_mask, x, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define gmx_simd4_dotproduct3_f(a, b) _mm512_mask_reduce_add_ps(_mm512_int2mask(7), _mm512_mask_mul_ps(_mm512_undefined_ps(), _mm512_int2mask(7), a, b))
#define gmx_simd4_fbool_t           __mmask16
#define gmx_simd4_cmpeq_f(a, b)     _mm512_mask_cmp_ps_mask(gmx_simd4_mask, a, b, _CMP_EQ_OQ)
#define gmx_simd4_cmplt_f(a, b)     _mm512_mask_cmp_ps_mask(gmx_simd4_mask, a, b, _CMP_LT_OS)
#define gmx_simd4_cmple_f(a, b)     _mm512_mask_cmp_ps_mask(gmx_simd4_mask, a, b, _CMP_LE_OS)
#define gmx_simd4_and_fb            _mm512_kand
#define gmx_simd4_or_fb             _mm512_kor
#define gmx_simd4_anytrue_fb(x)     (_mm512_mask2int(x)&0xF)
#define gmx_simd4_blendzero_f(a, sel)    _mm512_mask_mov_ps(_mm512_setzero_ps(), sel, a)
#define gmx_simd4_blendnotzero_f(a, sel) _mm512_mask_mov_ps(_mm512_setzero_ps(), _mm512_knot(sel), a)
#define gmx_simd4_blendv_f(a, b, sel)    _mm512_mask_blend_ps(sel, a, b)
#define gmx_simd4_reduce_f(x)       _mm512_mask_reduce_add_ps(_mm512_int2mask(0xF), x)

/* Implementation helpers */

/* load store simd4 */
static gmx_inline void gmx_simdcall
gmx_simd4_store_f_mic(float * m, __m512 s)
{
    assert((size_t)m%16 == 0);
    _mm512_mask_packstorelo_ps(m, gmx_simd4_mask, s);
}

static gmx_inline __m512 gmx_simdcall
gmx_simd4_loadu_f_mic(const float * m)
{
    return _mm512_mask_loadunpackhi_ps(_mm512_mask_loadunpacklo_ps(_mm512_undefined_ps(), gmx_simd4_mask, m), gmx_simd4_mask, m+16);
}

static gmx_inline void gmx_simdcall
gmx_simd4_storeu_f_mic(float * m, __m512 s)
{
    _mm512_mask_packstorelo_ps(m, gmx_simd4_mask, s);
    _mm512_mask_packstorehi_ps(m+16, gmx_simd4_mask, s);
}

#endif /* GMX_SIMD_IMPL_INTEL_MIC_SIMD4_FLOAT_H */
