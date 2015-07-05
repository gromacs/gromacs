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
#define gmx_simd4_double_t          __m512d
#define gmx_simd4_mask              _mm512_int2mask(0xF)
#define gmx_simd4_load_d(m)         _mm512_mask_extload_pd(_mm512_undefined_pd(), gmx_simd4_mask, m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE)
#define gmx_simd4_load1_d(m)        _mm512_mask_extload_pd(_mm512_undefined_pd(), gmx_simd4_mask, m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE)
#define gmx_simd4_set1_d            _mm512_set1_pd
#define gmx_simd4_store_d           gmx_simd4_store_d_mic
#define gmx_simd4_loadu_d           gmx_simd4_loadu_d_mic
#define gmx_simd4_storeu_d          gmx_simd4_storeu_d_mic
#define gmx_simd4_setzero_d         _mm512_setzero_pd
#define gmx_simd4_add_d(a, b)       _mm512_mask_add_pd(_mm512_undefined_pd(), gmx_simd4_mask, a, b)
#define gmx_simd4_sub_d(a, b)       _mm512_mask_sub_pd(_mm512_undefined_pd(), gmx_simd4_mask, a, b)
#define gmx_simd4_mul_d(a, b)       _mm512_mask_mul_pd(_mm512_undefined_pd(), gmx_simd4_mask, a, b)
#define gmx_simd4_fmadd_d(a, b, c)  _mm512_mask_fmadd_pd(a, gmx_simd4_mask, b, c)
#define gmx_simd4_fmsub_d(a, b, c)  _mm512_mask_fmsub_pd(a, gmx_simd4_mask, b, c)
#define gmx_simd4_fnmadd_d(a, b, c) _mm512_mask_fnmadd_pd(a, gmx_simd4_mask, b, c)
#define gmx_simd4_fnmsub_d(a, b, c) _mm512_mask_fnmsub_pd(a, gmx_simd4_mask, b, c)
#define gmx_simd4_and_d(a, b)       _mm512_castsi512_pd(_mm512_mask_and_epi32(_mm512_undefined_epi32(), mask_loh, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd4_andnot_d(a, b)    _mm512_castsi512_pd(_mm512_mask_andnot_epi32(_mm512_undefined_epi32(), mask_loh, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd4_or_d(a, b)        _mm512_castsi512_pd(_mm512_mask_or_epi32(_mm512_undefined_epi32(), mask_loh, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd4_xor_d(a, b)       _mm512_castsi512_pd(_mm512_mask_xor_epi32(_mm512_undefined_epi32(), mask_loh, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd4_rsqrt_d(a)        _mm512_mask_cvtpslo_pd(_mm512_undefined_pd(), gmx_simd4_mask, _mm512_mask_rsqrt23_ps(_mm512_undefined_ps(), gmx_simd4_mask, _mm512_mask_cvtpd_pslo(_mm512_undefined_ps(), gmx_simd4_mask, x)))
#define gmx_simd4_fabs_d(x)         gmx_simd4_andnot_d(_mm512_set1_pd(GMX_DOUBLE_NEGZERO), x)
#define gmx_simd4_fneg_d(x)         _mm512_mask_addn_pd(_mm512_undefined_pd(), gmx_simd4_mask, x, _mm512_setzero_pd())
#define gmx_simd4_max_d(a, b)       _mm512_mask_gmax_pd(_mm512_undefined_pd(), gmx_simd4_mask, a, b)
#define gmx_simd4_min_d(a, b)       _mm512_mask_gmin_pd(_mm512_undefined_pd(), gmx_simd4_mask, a, b)
#define gmx_simd4_round_d(a)        _mm512_mask_roundfxpnt_adjust_pd(_mm512_undefined_pd(), gmx_simd4_mask, a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define gmx_simd4_trunc_d(a)        _mm512_mask_roundfxpnt_adjust_pd(_mm512_undefined_pd(), gmx_simd4_mask, a, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define gmx_simd4_dotproduct3_d(a, b) _mm512_mask_reduce_add_pd(_mm512_int2mask(7), _mm512_mask_mul_pd(_mm512_undefined_pd(), _mm512_int2mask(7), a, b))
#define gmx_simd4_dbool_t           __mmask16
#define gmx_simd4_cmpeq_d(a, b)     _mm512_mask_cmp_pd_mask(gmx_simd4_mask, a, b, _CMP_EQ_OQ)
#define gmx_simd4_cmplt_d(a, b)     _mm512_mask_cmp_pd_mask(gmx_simd4_mask, a, b, _CMP_LT_OS)
#define gmx_simd4_cmple_d(a, b)     _mm512_mask_cmp_pd_mask(gmx_simd4_mask, a, b, _CMP_LE_OS)
#define gmx_simd4_and_db            _mm512_kand
#define gmx_simd4_or_db             _mm512_kor
#define gmx_simd4_anytrue_db(x)     (_mm512_mask2int(x)&0xF)
#define gmx_simd4_blendzero_d(a, sel)    _mm512_mask_mov_pd(_mm512_setzero_pd(), sel, a)
#define gmx_simd4_blendnotzero_d(a, sel) _mm512_mask_mov_pd(_mm512_setzero_pd(), _mm512_knot(sel), a)
#define gmx_simd4_blendv_d(a, b, sel)    _mm512_mask_blend_pd(sel, a, b)
#define gmx_simd4_reduce_d(x)       _mm512_mask_reduce_add_pd(_mm512_int2mask(0xF), x)

/* Implementation helpers */
static gmx_inline void gmx_simdcall
gmx_simd4_store_d_mic(double * m, __m512d s)
{
    assert((size_t)m%32 == 0);
    _mm512_mask_packstorelo_pd(m, gmx_simd4_mask, s);
}

static gmx_inline __m512d gmx_simdcall
gmx_simd4_loadu_d_mic(const double * m)
{
    return _mm512_mask_loadunpackhi_pd(_mm512_mask_loadunpacklo_pd(_mm512_undefined_pd(), gmx_simd4_mask, m), gmx_simd4_mask, m+8);
}

static gmx_inline void gmx_simdcall
gmx_simd4_storeu_d_mic(double * m, __m512d s)
{
    _mm512_mask_packstorelo_pd(m, gmx_simd4_mask, s);
    _mm512_mask_packstorehi_pd(m+8, gmx_simd4_mask, s);
}

#endif /* GMX_SIMD_IMPL_INTEL_MIC_SIMD4_DOUBLE_H */
