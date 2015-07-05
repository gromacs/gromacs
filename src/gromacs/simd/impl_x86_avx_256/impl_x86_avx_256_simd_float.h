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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_256_SIMD_FLOAT_H

#include "config.h"

#include <immintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_avx_256_common.h"

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_float_t           __m256
#define gmx_simd_load_f            _mm256_load_ps
#define gmx_simd_load1_f           _mm256_broadcast_ss
#define gmx_simd_set1_f            _mm256_set1_ps
#define gmx_simd_store_f           _mm256_store_ps
#define gmx_simd_loadu_f           _mm256_loadu_ps
#define gmx_simd_storeu_f          _mm256_storeu_ps
#define gmx_simd_setzero_f         _mm256_setzero_ps
#define gmx_simd_add_f             _mm256_add_ps
#define gmx_simd_sub_f             _mm256_sub_ps
#define gmx_simd_mul_f             _mm256_mul_ps
#define gmx_simd_fmadd_f(a, b, c)    _mm256_add_ps(_mm256_mul_ps(a, b), c)
#define gmx_simd_fmsub_f(a, b, c)    _mm256_sub_ps(_mm256_mul_ps(a, b), c)
#define gmx_simd_fnmadd_f(a, b, c)   _mm256_sub_ps(c, _mm256_mul_ps(a, b))
#define gmx_simd_fnmsub_f(a, b, c)   _mm256_sub_ps(_mm256_setzero_ps(), gmx_simd_fmadd_f(a, b, c))
#define gmx_simd_and_f             _mm256_and_ps
#define gmx_simd_andnot_f          _mm256_andnot_ps
#define gmx_simd_or_f              _mm256_or_ps
#define gmx_simd_xor_f             _mm256_xor_ps
#define gmx_simd_rsqrt_f           _mm256_rsqrt_ps
#define gmx_simd_rcp_f             _mm256_rcp_ps
#define gmx_simd_fabs_f(x)         _mm256_andnot_ps(_mm256_set1_ps(GMX_FLOAT_NEGZERO), x)
#define gmx_simd_fneg_f(x)         _mm256_xor_ps(x, _mm256_set1_ps(GMX_FLOAT_NEGZERO))
#define gmx_simd_max_f             _mm256_max_ps
#define gmx_simd_min_f             _mm256_min_ps
#define gmx_simd_round_f(x)        _mm256_round_ps(x, _MM_FROUND_NINT)
#define gmx_simd_trunc_f(x)        _mm256_round_ps(x, _MM_FROUND_TRUNC)
#define gmx_simd_fraction_f(x)     _mm256_sub_ps(x, gmx_simd_trunc_f(x))
#define gmx_simd_get_exponent_f    gmx_simd_get_exponent_f_avx_256
#define gmx_simd_get_mantissa_f    gmx_simd_get_mantissa_f_avx_256
#define gmx_simd_set_exponent_f    gmx_simd_set_exponent_f_avx_256
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_fint32_t          __m256i
#define gmx_simd_load_fi(m)        _mm256_load_si256((__m256i const*)(m))
#define gmx_simd_set1_fi           _mm256_set1_epi32
#define gmx_simd_store_fi(m, x)    _mm256_store_si256((__m256i *)(m), x)
#define gmx_simd_loadu_fi(m)       _mm256_loadu_si256((__m256i const*)(m))
#define gmx_simd_storeu_fi(m, x)   _mm256_storeu_si256((__m256i *)(m), x)
#define gmx_simd_setzero_fi        _mm256_setzero_si256
#define gmx_simd_cvt_f2i           _mm256_cvtps_epi32
#define gmx_simd_cvtt_f2i          _mm256_cvttps_epi32
#define gmx_simd_cvt_i2f           _mm256_cvtepi32_ps
#define gmx_simd_extract_fi(x, i)   _mm_extract_epi32(_mm256_extractf128_si256(x, (i)>>2), (i)&0x3)
/* Integer logical ops on gmx_simd_fint32_t */
/* gmx_simd_add_fi not supported     */
/* gmx_simd_sub_fi not supported     */
/* gmx_simd_mul_fi not supported     */
/* gmx_simd_slli_fi not supported    */
/* gmx_simd_srli_fi not supported    */
/* gmx_simd_and_fi not supported     */
/* gmx_simd_andnot_fi not supported  */
/* gmx_simd_or_fi not supported      */
/* gmx_simd_xor_fi not supported     */
/* Integer arithmetic ops on gmx_simd_fint32_t */
/* gmx_simd_add_fi not supported     */
/* gmx_simd_sub_fi not supported     */
/* gmx_simd_mul_fi not supported     */
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_fbool_t           __m256
#define gmx_simd_cmpeq_f(a, b)      _mm256_cmp_ps(a, b, _CMP_EQ_OQ)
#define gmx_simd_cmplt_f(a, b)      _mm256_cmp_ps(a, b, _CMP_LT_OQ)
#define gmx_simd_cmple_f(a, b)      _mm256_cmp_ps(a, b, _CMP_LE_OQ)
#define gmx_simd_and_fb            _mm256_and_ps
#define gmx_simd_or_fb             _mm256_or_ps
#define gmx_simd_anytrue_fb        _mm256_movemask_ps
#define gmx_simd_blendzero_f       _mm256_and_ps
#define gmx_simd_blendnotzero_f(a, sel)  _mm256_andnot_ps(sel, a)
#define gmx_simd_blendv_f          _mm256_blendv_ps
#define gmx_simd_reduce_f          gmx_simd_reduce_f_avx_256
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_fibool_t          __m256i
/* gmx_simd_cmpeq_fi not supported        */
/* gmx_simd_cmplt_fi not supported        */
/* gmx_simd_and_fib not supported         */
/* gmx_simd_or_fib not supported          */
/* gmx_simd_anytrue_fib not supported     */
/* gmx_simd_blendzero_fi not supported    */
/* gmx_simd_blendnotzero_fi not supported    */
/* gmx_simd_blendv_fi not supported       */
/* Conversions between different booleans */
#define gmx_simd_cvt_fb2fib        _mm256_castps_si256
#define gmx_simd_cvt_fib2fb        _mm256_castsi256_ps

/*********************************************************
 * SIMD SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 *********************************************************/
static gmx_inline __m256 gmx_simdcall
gmx_simd_get_exponent_f_avx_256(__m256 x)
{
    const __m256  expmask      = _mm256_castsi256_ps(_mm256_set1_epi32(0x7F800000));
    const __m128i expbias      = _mm_set1_epi32(127);
    __m256i       iexp256;
    __m128i       iexp128a, iexp128b;

    iexp256   = _mm256_castps_si256(_mm256_and_ps(x, expmask));
    iexp128b  = _mm256_extractf128_si256(iexp256, 0x1);
    iexp128a  = _mm256_castsi256_si128(iexp256);
    iexp128a  = _mm_srli_epi32(iexp128a, 23);
    iexp128b  = _mm_srli_epi32(iexp128b, 23);
    iexp128a  = _mm_sub_epi32(iexp128a, expbias);
    iexp128b  = _mm_sub_epi32(iexp128b, expbias);
    iexp256   = _mm256_castsi128_si256(iexp128a);
    iexp256   = _mm256_insertf128_si256(iexp256, iexp128b, 0x1);
    return _mm256_cvtepi32_ps(iexp256);
}

static gmx_inline __m256 gmx_simdcall
gmx_simd_get_mantissa_f_avx_256(__m256 x)
{
    const __m256 mantmask   = _mm256_castsi256_ps(_mm256_set1_epi32(0x007FFFFF));
    const __m256 one        = _mm256_set1_ps(1.0);

    x = _mm256_and_ps(x, mantmask);
    return _mm256_or_ps(x, one);
}

static gmx_inline __m256 gmx_simdcall
gmx_simd_set_exponent_f_avx_256(__m256 x)
{
    const __m128i expbias      = _mm_set1_epi32(127);
    __m256i       iexp256;
    __m128i       iexp128a, iexp128b;

    iexp256   = _mm256_cvtps_epi32(x);
    iexp128b  = _mm256_extractf128_si256(iexp256, 0x1);
    iexp128a  = _mm256_castsi256_si128(iexp256);
    iexp128a  = _mm_slli_epi32(_mm_add_epi32(iexp128a, expbias), 23);
    iexp128b  = _mm_slli_epi32(_mm_add_epi32(iexp128b, expbias), 23);
    iexp256   = _mm256_castsi128_si256(iexp128a);
    iexp256   = _mm256_insertf128_si256(iexp256, iexp128b, 0x1);
    return _mm256_castsi256_ps(iexp256);
}

static gmx_inline float gmx_simdcall
gmx_simd_reduce_f_avx_256(__m256 a)
{
    float  f;

    __m128 a0, a1;
    a  = _mm256_hadd_ps(a, a);
    a  = _mm256_hadd_ps(a, a);
    a0 = _mm256_castps256_ps128(a);
    a1 = _mm256_extractf128_ps(a, 0x1);
    a0 = _mm_add_ss(a0, a1);
    _mm_store_ss(&f, a0);
    return f;
}

#endif /* GMX_SIMD_IMPL_X86_AVX_256_SIMD_FLOAT_H */
