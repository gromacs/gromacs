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

#ifndef GMX_SIMD_IMPL_X86_AVX_512F_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_512F_SIMD_FLOAT_H

#include <math.h>

#include <immintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_avx_512f_common.h"

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_float_t           __m512
#define gmx_simd_load_f            _mm512_load_ps
/* Avoid using _mm512_extload_ps() since it is not available on gcc-4.9 */
#define gmx_simd_load1_f(m)        _mm512_broadcastss_ps(_mm_broadcast_ss(m))
#define gmx_simd_set1_f            _mm512_set1_ps
#define gmx_simd_store_f           _mm512_store_ps
#define gmx_simd_loadu_f           _mm512_loadu_ps
#define gmx_simd_storeu_f          _mm512_storeu_ps
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
#define gmx_simd_rsqrt_f           _mm512_rsqrt14_ps
#define gmx_simd_rcp_f             _mm512_rcp14_ps
#define gmx_simd_fabs_f(x)         _mm512_abs_ps(x)
#define gmx_simd_fneg_f(x)         gmx_simd_xor_f(x, _mm512_set1_ps(GMX_FLOAT_NEGZERO))
#define gmx_simd_max_f             _mm512_max_ps
#define gmx_simd_min_f             _mm512_min_ps
#define gmx_simd_round_f(x)        _mm512_roundscale_ps(x, 0)
#define gmx_simd_trunc_f(x)        _mm512_trunc_ps(x)
#define gmx_simd_fraction_f(x)     _mm512_sub_ps(x, gmx_simd_trunc_f(x))
#define gmx_simd_get_exponent_f(x) _mm512_getexp_ps(x)
#define gmx_simd_get_mantissa_f(x) _mm512_getmant_ps(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero)
#define gmx_simd_set_exponent_f(x) gmx_simd_set_exponent_f_x86_avx_512f(x)
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_fint32_t          __m512i
#define gmx_simd_load_fi           _mm512_load_si512
#define gmx_simd_set1_fi           _mm512_set1_epi32
#define gmx_simd_store_fi          _mm512_store_si512
#define gmx_simd_loadu_fi          _mm512_loadu_si512
#define gmx_simd_storeu_fi         _mm512_storeu_si512
#undef  gmx_simd_extract_fi
#define gmx_simd_setzero_fi        _mm512_setzero_epi32
#define gmx_simd_cvt_f2i           _mm512_cvtps_epi32
#define gmx_simd_cvtt_f2i          _mm512_cvttps_epi32
#define gmx_simd_cvt_i2f           _mm512_cvtepi32_ps
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
#define gmx_simd_andnot_fb(a, b)   _mm512_kandn(a, b)
#define gmx_simd_or_fb             _mm512_kor
#define gmx_simd_anytrue_fb        _mm512_mask2int
#define gmx_simd_blendzero_f(a, sel)    _mm512_mask_mov_ps(_mm512_setzero_ps(), sel, a)
#define gmx_simd_blendnotzero_f(a, sel) _mm512_mask_mov_ps(a, sel, _mm512_setzero_ps())
#define gmx_simd_blendv_f(a, b, sel)    _mm512_mask_blend_ps(sel, a, b)
#define gmx_simd_reduce_f(a)       gmx_simd_reduce_f_x86_avx_512f(a)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_fibool_t          __mmask16
#define gmx_simd_cmpeq_fi(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_EQ)
#define gmx_simd_cmplt_fi(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_LT)
#define gmx_simd_and_fib           _mm512_kand
#define gmx_simd_or_fib            _mm512_kor
#define gmx_simd_anytrue_fib       _mm512_mask2int
#define gmx_simd_blendzero_fi(a, sel)    _mm512_mask_mov_epi32(_mm512_setzero_epi32(), sel, a)
#define gmx_simd_blendnotzero_fi(a, sel) _mm512_mask_mov_epi32(a, sel, _mm512_setzero_epi32())
#define gmx_simd_blendv_fi(a, b, sel)    _mm512_mask_blend_epi32(sel, a, b)
/* Conversions between different booleans */
#define gmx_simd_cvt_fb2fib(x)     (x)
#define gmx_simd_cvt_fib2fb(x)     (x)


/* Implementation helper functions */

/* This is likely faster than the built in scale operation (lat 8, t-put 3)
 * since we only work on the integer part and use shifts, meaning it will run
 * on the integer ports that are typically less utilized in our kernels.
 */
static gmx_inline __m512
gmx_simd_set_exponent_f_x86_avx_512f(__m512 a)
{
    const __m512i expbias      = _mm512_set1_epi32(127);
    __m512i       iexp         = gmx_simd_cvt_f2i(a);

    iexp = _mm512_slli_epi32(_mm512_add_epi32(iexp, expbias), 23);
    return _mm512_castsi512_ps(iexp);
}

static gmx_inline float
gmx_simd_reduce_f_x86_avx_512f(__m512 a)
{
    __m128 b;
    a = _mm512_add_ps(a, _mm512_shuffle_f32x4(a, a, _MM_PERM_DCDC));
    a = _mm512_add_ps(a, _mm512_shuffle_f32x4(a, a, _MM_PERM_ABAB));
    b = _mm512_castps512_ps128(a);
    b = _mm_hadd_ps(b, b);
    b = _mm_hadd_ps(b, b);
    return _mm_cvtss_f32(b);
}

#endif /* GMX_SIMD_IMPL_X86_AVX_512F_SIMD_FLOAT_H */
