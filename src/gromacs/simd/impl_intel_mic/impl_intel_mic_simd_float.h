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

#include <math.h>

#include <immintrin.h>

#include "impl_intel_mic_common.h"

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_float_t           __m512
#define gmx_simd_load_f            _mm512_load_ps
#define gmx_simd_load1_f(m)        _mm512_extload_ps(m, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE)
#define gmx_simd_set1_f            _mm512_set1_ps
#define gmx_simd_store_f           _mm512_store_ps
#define gmx_simd_loadu_f           gmx_simd_loadu_f_mic
#define gmx_simd_storeu_f          gmx_simd_storeu_f_mic
#define gmx_simd_setzero_f         _mm512_setzero_ps
#define gmx_simd_add_f             _mm512_add_ps
#define gmx_simd_sub_f             _mm512_sub_ps
#define gmx_simd_mul_f             _mm512_mul_ps
#define gmx_simd_fmadd_f           _mm512_fmadd_ps
#define gmx_simd_fmsub_f           _mm512_fmsub_ps
#define gmx_simd_fnmadd_f          _mm512_fnmadd_ps
#define gmx_simd_fnmsub_f          _mm512_fnmsub_ps
#define gmx_simd_and_f(a, b)        _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd_andnot_f(a, b)     _mm512_castsi512_ps(_mm512_andnot_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd_or_f(a, b)         _mm512_castsi512_ps(_mm512_or_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd_xor_f(a, b)        _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd_rsqrt_f           _mm512_rsqrt23_ps
#define gmx_simd_rcp_f             _mm512_rcp23_ps
#define gmx_simd_fabs_f(x)         gmx_simd_andnot_f(_mm512_set1_ps(GMX_FLOAT_NEGZERO), x)
#define gmx_simd_fneg_f(x)         _mm512_addn_ps(x, _mm512_setzero_ps())
#define gmx_simd_max_f             _mm512_gmax_ps
#define gmx_simd_min_f             _mm512_gmin_ps
#define gmx_simd_round_f(x)        _mm512_round_ps(x, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define gmx_simd_trunc_f(x)        _mm512_round_ps(x, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define gmx_simd_fraction_f(x)     _mm512_sub_ps(x, gmx_simd_trunc_f(x))
#define gmx_simd_get_exponent_f(x) _mm512_getexp_ps(x)
#define gmx_simd_get_mantissa_f(x) _mm512_getmant_ps(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero)
#define gmx_simd_set_exponent_f(x) gmx_simd_set_exponent_f_mic(x)
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_fint32_t          __m512i
#define gmx_simd_load_fi           _mm512_load_epi32
#define gmx_simd_set1_fi           _mm512_set1_epi32
#define gmx_simd_store_fi          _mm512_store_epi32
#define gmx_simd_loadu_fi          gmx_simd_loadu_fi_mic
#define gmx_simd_storeu_fi         gmx_simd_storeu_fi_mic
#define gmx_simd_extract_fi        gmx_simd_extract_fi_mic
#define gmx_simd_setzero_fi        _mm512_setzero_epi32
#define gmx_simd_cvt_f2i(a)        _mm512_cvtfxpnt_round_adjustps_epi32(a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define gmx_simd_cvtt_f2i(a)       _mm512_cvtfxpnt_round_adjustps_epi32(a, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define gmx_simd_cvt_i2f(a)        _mm512_cvtfxpnt_round_adjustepi32_ps(a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
/* Integer logical ops on gmx_simd_fint32_t */
#define gmx_simd_slli_fi           _mm512_slli_epi32
#define gmx_simd_srli_fi           _mm512_srli_epi32
#define gmx_simd_and_fi            _mm512_and_epi32
#define gmx_simd_andnot_fi         _mm512_andnot_epi32
#define gmx_simd_or_fi             _mm512_or_epi32
#define gmx_simd_xor_fi            _mm512_xor_epi32
/* Integer arithmetic ops on gmx_simd_fint32_t */
#define gmx_simd_add_fi            _mm512_add_epi32
#define gmx_simd_sub_fi            _mm512_sub_epi32
#define gmx_simd_mul_fi            _mm512_mullo_epi32
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_fbool_t           __mmask16
#define gmx_simd_cmpeq_f(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_EQ_OQ)
#define gmx_simd_cmplt_f(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_LT_OS)
#define gmx_simd_cmple_f(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_LE_OS)
#define gmx_simd_and_fb            _mm512_kand
#define gmx_simd_andnot_fb(a, b)   _mm512_knot(_mm512_kor(a, b))
#define gmx_simd_or_fb             _mm512_kor
#define gmx_simd_anytrue_fb        _mm512_mask2int
#define gmx_simd_blendzero_f(a, sel)    _mm512_mask_mov_ps(_mm512_setzero_ps(), sel, a)
#define gmx_simd_blendnotzero_f(a, sel) _mm512_mask_mov_ps(_mm512_setzero_ps(), _mm512_knot(sel), a)
#define gmx_simd_blendv_f(a, b, sel)    _mm512_mask_blend_ps(sel, a, b)
#define gmx_simd_reduce_f(a)       _mm512_reduce_add_ps(a)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_fibool_t          __mmask16
#define gmx_simd_cmpeq_fi(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_EQ)
#define gmx_simd_cmplt_fi(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_LT)
#define gmx_simd_and_fib           _mm512_kand
#define gmx_simd_or_fib            _mm512_kor
#define gmx_simd_anytrue_fib       _mm512_mask2int
#define gmx_simd_blendzero_fi(a, sel)    _mm512_mask_mov_epi32(_mm512_setzero_epi32(), sel, a)
#define gmx_simd_blendnotzero_fi(a, sel) _mm512_mask_mov_epi32(_mm512_setzero_epi32(), _mm512_knot(sel), a)
#define gmx_simd_blendv_fi(a, b, sel)    _mm512_mask_blend_epi32(sel, a, b)
/* Conversions between different booleans */
#define gmx_simd_cvt_fb2fib(x)     (x)
#define gmx_simd_cvt_fib2fb(x)     (x)

/* MIC provides full single precision of some neat functions: */
/* 1/sqrt(x) and 1/x work fine in simd_math.h, and won't use extra iterations */
#define gmx_simd_exp2_f            gmx_simd_exp2_f_mic
#define gmx_simd_exp_f             gmx_simd_exp_f_mic
#define gmx_simd_log_f             gmx_simd_log_f_mic

/* load store float */
static gmx_inline __m512 gmx_simdcall
gmx_simd_loadu_f_mic(const float * m)
{
    return _mm512_loadunpackhi_ps(_mm512_loadunpacklo_ps(_mm512_undefined_ps(), m), m+16);
}

static gmx_inline void gmx_simdcall
gmx_simd_storeu_f_mic(float * m, __m512 s)
{
    _mm512_packstorelo_ps(m, s);
    _mm512_packstorehi_ps(m+16, s);
}

/* load store fint32 */
static gmx_inline __m512i gmx_simdcall
gmx_simd_loadu_fi_mic(const gmx_int32_t * m)
{
    return _mm512_loadunpackhi_epi32(_mm512_loadunpacklo_epi32(_mm512_undefined_epi32(), m), m+16);
}

static gmx_inline void gmx_simdcall
gmx_simd_storeu_fi_mic(gmx_int32_t * m, __m512i s)
{
    _mm512_packstorelo_epi32(m, s);
    _mm512_packstorehi_epi32(m+16, s);
}

/* extract */
static gmx_inline gmx_int32_t gmx_simdcall
gmx_simd_extract_fi_mic(gmx_simd_fint32_t a, int index)
{
    int r;
    _mm512_mask_packstorelo_epi32(&r, _mm512_mask2int(1<<index), a);
    return r;
}

/* This is likely faster than the built in scale operation (lat 8, t-put 3)
 * since we only work on the integer part and use shifts. TODO: check. given that scale also only does integer
 */
static gmx_inline __m512 gmx_simdcall
gmx_simd_set_exponent_f_mic(__m512 a)
{
    __m512i       iexp         = gmx_simd_cvt_f2i(a);

    const __m512i expbias      = _mm512_set1_epi32(127);
    iexp = _mm512_slli_epi32(_mm512_add_epi32(iexp, expbias), 23);
    return _mm512_castsi512_ps(iexp);

    /* scale alternative:
       return _mm512_scale_ps(_mm512_set1_ps(1), iexp);
     */
}

static gmx_inline __m512 gmx_simdcall
gmx_simd_exp2_f_mic(__m512 x)
{
    return _mm512_exp223_ps(_mm512_cvtfxpnt_round_adjustps_epi32(x, _MM_ROUND_MODE_NEAREST, _MM_EXPADJ_24));
}

static gmx_inline __m512 gmx_simdcall
gmx_simd_exp_f_mic(__m512 x)
{
    const gmx_simd_float_t  argscale    = gmx_simd_set1_f(1.44269504088896341f);
    const gmx_simd_float_t  invargscale = gmx_simd_set1_f(-0.69314718055994528623f);
    __m512                  xscaled     = _mm512_mul_ps(x, argscale);
    __m512                  r           = gmx_simd_exp2_f_mic(xscaled);

    /* gmx_simd_exp2_f_mic() provides 23 bits of accuracy, but we ruin some of that
     * with the argument scaling due to single-precision rounding, where the
     * rounding error is amplified exponentially. To correct this, we find the
     * difference between the scaled argument and the true one (extended precision
     * arithmetics does not appear to be necessary to fulfill our accuracy requirements)
     * and then multiply by the exponent of this correction since exp(a+b)=exp(a)*exp(b).
     * Note that this only adds two instructions (and maybe some constant loads).
     */
    x         = gmx_simd_fmadd_f(invargscale, xscaled, x);
    /* x will now be a _very_ small number, so approximate exp(x)=1+x.
     * We should thus apply the correction as r'=r*(1+x)=r+r*x
     */
    r         = gmx_simd_fmadd_f(r, x, r);
    return r;
}

static gmx_inline __m512 gmx_simdcall
gmx_simd_log_f_mic(__m512 x)
{
    return _mm512_mul_ps(_mm512_set1_ps(0.693147180559945286226764), _mm512_log2ae23_ps(x));
}

#endif /* GMX_SIMD_IMPL_INTEL_MIC_SIMD_FLOAT_H */
