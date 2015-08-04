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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_H
#define GMX_SIMD_IMPL_X86_AVX_256_H

#include "config.h"

#include <immintrin.h>

/* It is cleaner to start the AVX implementation from scratch rather than
 * first inheriting from SSE4.1, which in turn inherits from SSE2. However,
 * the capabilities still form a superset.
 */
#define GMX_SIMD_X86_SSE2_OR_HIGHER
#define GMX_SIMD_X86_SSE4_1_OR_HIGHER
#define GMX_SIMD_X86_AVX_256_OR_HIGHER


/* x86 256-bit AVX SIMD instruction wrappers
 *
 * Please see documentation in gromacs/simd/simd.h for defines.
 */

/* Capability definitions for 256-bit AVX - no inheritance from SSE */
#define GMX_SIMD_HAVE_FLOAT
#define GMX_SIMD_HAVE_DOUBLE
#define GMX_SIMD_HAVE_SIMD_HARDWARE
#define GMX_SIMD_HAVE_LOADU
#define GMX_SIMD_HAVE_STOREU
#define GMX_SIMD_HAVE_LOGICAL
#undef  GMX_SIMD_HAVE_FMA
#undef  GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD_HAVE_FINT32
#define GMX_SIMD_HAVE_FINT32_EXTRACT     /* Emulated */
#undef  GMX_SIMD_HAVE_FINT32_LOGICAL     /* AVX1 cannot do 256-bit int shifts */
#undef  GMX_SIMD_HAVE_FINT32_ARITHMETICS /* AVX1 cannot do 256-bit int +,-,*  */
#define GMX_SIMD_HAVE_DINT32
#define GMX_SIMD_HAVE_DINT32_EXTRACT     /* Native, dint uses 128-bit SIMD    */
#define GMX_SIMD_HAVE_DINT32_LOGICAL
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS
#define GMX_SIMD4_HAVE_FLOAT
#define GMX_SIMD4_HAVE_DOUBLE

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH         8
#define GMX_SIMD_DOUBLE_WIDTH        4
#define GMX_SIMD_FINT32_WIDTH        8
#define GMX_SIMD_DINT32_WIDTH        4
#define GMX_SIMD_RSQRT_BITS         11
#define GMX_SIMD_RCP_BITS           11

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

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_double_t          __m256d
#define gmx_simd_load_d            _mm256_load_pd
#define gmx_simd_load1_d           _mm256_broadcast_sd
#define gmx_simd_set1_d            _mm256_set1_pd
#define gmx_simd_store_d           _mm256_store_pd
#define gmx_simd_loadu_d           _mm256_loadu_pd
#define gmx_simd_storeu_d          _mm256_storeu_pd
#define gmx_simd_setzero_d         _mm256_setzero_pd
#define gmx_simd_add_d             _mm256_add_pd
#define gmx_simd_sub_d             _mm256_sub_pd
#define gmx_simd_mul_d             _mm256_mul_pd
#define gmx_simd_fmadd_d(a, b, c)    _mm256_add_pd(_mm256_mul_pd(a, b), c)
#define gmx_simd_fmsub_d(a, b, c)    _mm256_sub_pd(_mm256_mul_pd(a, b), c)
#define gmx_simd_fnmadd_d(a, b, c)   _mm256_sub_pd(c, _mm256_mul_pd(a, b))
#define gmx_simd_fnmsub_d(a, b, c)   _mm256_sub_pd(_mm256_setzero_pd(), gmx_simd_fmadd_d(a, b, c))
#define gmx_simd_and_d             _mm256_and_pd
#define gmx_simd_andnot_d          _mm256_andnot_pd
#define gmx_simd_or_d              _mm256_or_pd
#define gmx_simd_xor_d             _mm256_xor_pd
#define gmx_simd_rsqrt_d(x)        _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x)))
#define gmx_simd_rcp_d(x)          _mm256_cvtps_pd(_mm_rcp_ps(_mm256_cvtpd_ps(x)))
#define gmx_simd_fabs_d(x)         _mm256_andnot_pd(_mm256_set1_pd(-0.0), x)
#define gmx_simd_fneg_d(x)         _mm256_xor_pd(x, _mm256_set1_pd(-0.0))
#define gmx_simd_max_d             _mm256_max_pd
#define gmx_simd_min_d             _mm256_min_pd
#define gmx_simd_round_d(x)        _mm256_round_pd(x, _MM_FROUND_NINT)
#define gmx_simd_trunc_d(x)        _mm256_round_pd(x, _MM_FROUND_TRUNC)
#define gmx_simd_fraction_d(x)     _mm256_sub_pd(x, gmx_simd_trunc_d(x))
#define gmx_simd_get_exponent_d    gmx_simd_get_exponent_d_avx_256
#define gmx_simd_get_mantissa_d    gmx_simd_get_mantissa_d_avx_256
#define gmx_simd_set_exponent_d    gmx_simd_set_exponent_d_avx_256
/* integer datatype corresponding to double: gmx_simd_dint32_t */
#define gmx_simd_dint32_t          __m128i
#define gmx_simd_load_di(m)        _mm_load_si128((const __m128i *)(m))
#define gmx_simd_set1_di           _mm_set1_epi32
#define gmx_simd_store_di(m, x)     _mm_store_si128((__m128i *)(m), x)
#define gmx_simd_loadu_di(m)       _mm_loadu_si128((const __m128i *)(m))
#define gmx_simd_storeu_di(m, x)    _mm_storeu_si128((__m128i *)(m), x)
#define gmx_simd_setzero_di        _mm_setzero_si128
#define gmx_simd_cvt_d2i           _mm256_cvtpd_epi32
#define gmx_simd_cvtt_d2i          _mm256_cvttpd_epi32
#define gmx_simd_cvt_i2d           _mm256_cvtepi32_pd
#define gmx_simd_extract_di        _mm_extract_epi32
/* Integer logical ops on gmx_simd_dint32_t */
#define gmx_simd_slli_di           _mm_slli_epi32
#define gmx_simd_srli_di           _mm_srli_epi32
#define gmx_simd_and_di            _mm_and_si128
#define gmx_simd_andnot_di         _mm_andnot_si128
#define gmx_simd_or_di             _mm_or_si128
#define gmx_simd_xor_di            _mm_xor_si128
/* Integer arithmetic ops on integer datatype corresponding to double */
#define gmx_simd_add_di            _mm_add_epi32
#define gmx_simd_sub_di            _mm_sub_epi32
#define gmx_simd_mul_di            _mm_mullo_epi32
/* Boolean & comparison operations on gmx_simd_double_t */
#define gmx_simd_dbool_t           __m256d
#define gmx_simd_cmpeq_d(a, b)      _mm256_cmp_pd(a, b, _CMP_EQ_OQ)
#define gmx_simd_cmplt_d(a, b)      _mm256_cmp_pd(a, b, _CMP_LT_OQ)
#define gmx_simd_cmple_d(a, b)      _mm256_cmp_pd(a, b, _CMP_LE_OQ)
#define gmx_simd_and_db            _mm256_and_pd
#define gmx_simd_or_db             _mm256_or_pd
#define gmx_simd_anytrue_db        _mm256_movemask_pd
#define gmx_simd_blendzero_d       _mm256_and_pd
#define gmx_simd_blendnotzero_d(a, sel)  _mm256_andnot_pd(sel, a)
#define gmx_simd_blendv_d          _mm256_blendv_pd
#define gmx_simd_reduce_d          gmx_simd_reduce_d_avx_256
/* Boolean & comparison operations on gmx_simd_dint32_t */
#define gmx_simd_dibool_t          __m128i
#define gmx_simd_cmpeq_di          _mm_cmpeq_epi32
#define gmx_simd_cmplt_di          _mm_cmplt_epi32
#define gmx_simd_and_dib           _mm_and_si128
#define gmx_simd_or_dib            _mm_or_si128
#define gmx_simd_anytrue_dib       _mm_movemask_epi8
#define gmx_simd_blendzero_di      _mm_and_si128
#define gmx_simd_blendnotzero_di(a, sel)  _mm_andnot_si128(sel, a)
#define gmx_simd_blendv_di         _mm_blendv_epi8
/* Conversions between different booleans */
#define gmx_simd_cvt_db2dib        gmx_simd_cvt_db2dib_avx_256
#define gmx_simd_cvt_dib2db        gmx_simd_cvt_dib2db_avx_256
/* Float/double conversion */
#define gmx_simd_cvt_f2dd          gmx_simd_cvt_f2dd_avx_256
#define gmx_simd_cvt_dd2f          gmx_simd_cvt_dd2f_avx_256

/****************************************************
 *      SINGLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
#define gmx_simd4_float_t          __m128
#define gmx_simd4_load_f           _mm_load_ps
#define gmx_simd4_load1_f          _mm_broadcast_ss
#define gmx_simd4_set1_f           _mm_set1_ps
#define gmx_simd4_store_f          _mm_store_ps
#define gmx_simd4_loadu_f          _mm_loadu_ps
#define gmx_simd4_storeu_f         _mm_storeu_ps
#define gmx_simd4_setzero_f        _mm_setzero_ps
#define gmx_simd4_add_f            _mm_add_ps
#define gmx_simd4_sub_f            _mm_sub_ps
#define gmx_simd4_mul_f            _mm_mul_ps
#define gmx_simd4_fmadd_f(a, b, c)   _mm_add_ps(_mm_mul_ps(a, b), c)
#define gmx_simd4_fmsub_f(a, b, c)   _mm_sub_ps(_mm_mul_ps(a, b), c)
#define gmx_simd4_fnmadd_f(a, b, c)  _mm_sub_ps(c, _mm_mul_ps(a, b))
#define gmx_simd4_fnmsub_f(a, b, c)  _mm_sub_ps(_mm_setzero_ps(), gmx_simd4_fmadd_f(a, b, c))
#define gmx_simd4_and_f            _mm_and_ps
#define gmx_simd4_andnot_f         _mm_andnot_ps
#define gmx_simd4_or_f             _mm_or_ps
#define gmx_simd4_xor_f            _mm_xor_ps
#define gmx_simd4_rsqrt_f          _mm_rsqrt_ps
#define gmx_simd4_fabs_f(x)        _mm_andnot_ps(_mm_set1_ps(-0.0), x)
#define gmx_simd4_fneg_f(x)        _mm_xor_ps(x, _mm_set1_ps(-0.0))
#define gmx_simd4_max_f            _mm_max_ps
#define gmx_simd4_min_f            _mm_min_ps
#define gmx_simd4_round_f(x)       _mm_round_ps(x, _MM_FROUND_NINT)
#define gmx_simd4_trunc_f(x)       _mm_round_ps(x, _MM_FROUND_TRUNC)
#define gmx_simd4_dotproduct3_f    gmx_simd4_dotproduct3_f_avx_256
#define gmx_simd4_fbool_t          __m128
#define gmx_simd4_cmpeq_f          _mm_cmpeq_ps
#define gmx_simd4_cmplt_f          _mm_cmplt_ps
#define gmx_simd4_cmple_f          _mm_cmple_ps
#define gmx_simd4_and_fb           _mm_and_ps
#define gmx_simd4_or_fb            _mm_or_ps
#define gmx_simd4_anytrue_fb       _mm_movemask_ps
#define gmx_simd4_blendzero_f      _mm_and_ps
#define gmx_simd4_blendnotzero_f(a, sel)  _mm_andnot_ps(sel, a)
#define gmx_simd4_blendv_f         _mm_blendv_ps
#define gmx_simd4_reduce_f         gmx_simd4_reduce_f_avx_256

/****************************************************
 *      DOUBLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
#define gmx_simd4_double_t          gmx_simd_double_t
#define gmx_simd4_load_d            gmx_simd_load_d
#define gmx_simd4_load1_d           gmx_simd_load1_d
#define gmx_simd4_set1_d            gmx_simd_set1_d
#define gmx_simd4_store_d           gmx_simd_store_d
#define gmx_simd4_loadu_d           gmx_simd_loadu_d
#define gmx_simd4_storeu_d          gmx_simd_storeu_d
#define gmx_simd4_setzero_d         gmx_simd_setzero_d
#define gmx_simd4_add_d             gmx_simd_add_d
#define gmx_simd4_sub_d             gmx_simd_sub_d
#define gmx_simd4_mul_d             gmx_simd_mul_d
#define gmx_simd4_fmadd_d           gmx_simd_fmadd_d
#define gmx_simd4_fmsub_d           gmx_simd_fmsub_d
#define gmx_simd4_fnmadd_d          gmx_simd_fnmadd_d
#define gmx_simd4_fnmsub_d          gmx_simd_fnmsub_d
#define gmx_simd4_and_d             gmx_simd_and_d
#define gmx_simd4_andnot_d          gmx_simd_andnot_d
#define gmx_simd4_or_d              gmx_simd_or_d
#define gmx_simd4_xor_d             gmx_simd_xor_d
#define gmx_simd4_rsqrt_d           gmx_simd_rsqrt_d
#define gmx_simd4_fabs_d            gmx_simd_fabs_d
#define gmx_simd4_fneg_d            gmx_simd_fneg_d
#define gmx_simd4_max_d             gmx_simd_max_d
#define gmx_simd4_min_d             gmx_simd_min_d
#define gmx_simd4_round_d           gmx_simd_round_d
#define gmx_simd4_trunc_d           gmx_simd_trunc_d
#define gmx_simd4_dotproduct3_d     gmx_simd4_dotproduct3_d_avx_256
#define gmx_simd4_dbool_t           gmx_simd_dbool_t
#define gmx_simd4_cmpeq_d           gmx_simd_cmpeq_d
#define gmx_simd4_cmplt_d           gmx_simd_cmplt_d
#define gmx_simd4_cmple_d           gmx_simd_cmple_d
#define gmx_simd4_and_db            gmx_simd_and_db
#define gmx_simd4_or_db             gmx_simd_or_db
#define gmx_simd4_anytrue_db        gmx_simd_anytrue_db
#define gmx_simd4_blendzero_d       gmx_simd_blendzero_d
#define gmx_simd4_blendnotzero_d    gmx_simd_blendnotzero_d
#define gmx_simd4_blendv_d          gmx_simd_blendv_d
#define gmx_simd4_reduce_d          gmx_simd_reduce_d
/* SIMD4 float/double conversion */
#define gmx_simd4_cvt_f2d           _mm256_cvtps_pd
#define gmx_simd4_cvt_d2f           _mm256_cvtpd_ps

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

/*********************************************************
 * SIMD DOUBLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 *********************************************************/
static gmx_inline __m256d gmx_simdcall
gmx_simd_get_exponent_d_avx_256(__m256d x)
{
    const __m256d expmask      = _mm256_castsi256_pd( _mm256_set1_epi64x(0x7FF0000000000000LL));
    const __m128i expbias      = _mm_set1_epi32(1023);
    __m256i       iexp256;
    __m128i       iexp128a, iexp128b;

    iexp256   = _mm256_castpd_si256(_mm256_and_pd(x, expmask));
    iexp128b  = _mm256_extractf128_si256(iexp256, 0x1);
    iexp128a  = _mm256_castsi256_si128(iexp256);
    iexp128a  = _mm_srli_epi64(iexp128a, 52);
    iexp128b  = _mm_srli_epi64(iexp128b, 52);
    iexp128a  = _mm_shuffle_epi32(iexp128a, _MM_SHUFFLE(1, 1, 2, 0));
    iexp128b  = _mm_shuffle_epi32(iexp128b, _MM_SHUFFLE(2, 0, 1, 1));
    iexp128a  = _mm_or_si128(iexp128a, iexp128b);
    iexp128a  = _mm_sub_epi32(iexp128a, expbias);
    return _mm256_cvtepi32_pd(iexp128a);
}

static gmx_inline __m256d gmx_simdcall
gmx_simd_get_mantissa_d_avx_256(__m256d x)
{
    const __m256d mantmask   = _mm256_castsi256_pd(_mm256_set1_epi64x(0x000FFFFFFFFFFFFFLL));
    const __m256d one        = _mm256_set1_pd(1.0);

    x = _mm256_and_pd(x, mantmask);
    return _mm256_or_pd(x, one);
}

static gmx_inline __m256d gmx_simdcall
gmx_simd_set_exponent_d_avx_256(__m256d x)
{
    const __m128i expbias      = _mm_set1_epi32(1023);
    __m128i       iexp128a, iexp128b;

    iexp128a = _mm256_cvtpd_epi32(x);
    iexp128a = _mm_add_epi32(iexp128a, expbias);
    iexp128b = _mm_shuffle_epi32(iexp128a, _MM_SHUFFLE(3, 3, 2, 2));
    iexp128a = _mm_shuffle_epi32(iexp128a, _MM_SHUFFLE(1, 1, 0, 0));
    iexp128b = _mm_slli_epi64(iexp128b, 52);
    iexp128a = _mm_slli_epi64(iexp128a, 52);
    return _mm256_castsi256_pd(_mm256_insertf128_si256(_mm256_castsi128_si256(iexp128a), iexp128b, 0x1));
}

static gmx_inline double gmx_simdcall
gmx_simd_reduce_d_avx_256(__m256d a)
{
    double  f;
    __m128d a0, a1;
    a  = _mm256_hadd_pd(a, a);
    a0 = _mm256_castpd256_pd128(a);
    a1 = _mm256_extractf128_pd(a, 0x1);
    a0 = _mm_add_sd(a0, a1);
    _mm_store_sd(&f, a0);
    return f;
}

static gmx_inline gmx_simd_dibool_t gmx_simdcall
gmx_simd_cvt_db2dib_avx_256(gmx_simd_dbool_t a)
{
    __m128i a1 = _mm256_extractf128_si256(_mm256_castpd_si256(a), 0x1);
    __m128i a0 = _mm256_castsi256_si128(_mm256_castpd_si256(a));
    a0 = _mm_shuffle_epi32(a0, _MM_SHUFFLE(2, 0, 2, 0));
    a1 = _mm_shuffle_epi32(a1, _MM_SHUFFLE(2, 0, 2, 0));
    return _mm_blend_epi16(a0, a1, 0xF0);
}

static gmx_inline gmx_simd_dbool_t gmx_simdcall
gmx_simd_cvt_dib2db_avx_256(gmx_simd_dibool_t a)
{
    __m128i a1 = _mm_shuffle_epi32(a, _MM_SHUFFLE(3, 3, 2, 2));
    __m128i a0 = _mm_shuffle_epi32(a, _MM_SHUFFLE(1, 1, 0, 0));
    return _mm256_castsi256_pd(_mm256_insertf128_si256(_mm256_castsi128_si256(a0), a1, 0x1));
}

static gmx_inline void gmx_simdcall
gmx_simd_cvt_f2dd_avx_256(__m256 f, __m256d *d0, __m256d *d1)
{
    *d0 = _mm256_cvtps_pd(_mm256_castps256_ps128(f));
    *d1 = _mm256_cvtps_pd(_mm256_extractf128_ps(f, 0x1));
}

static gmx_inline __m256 gmx_simdcall
gmx_simd_cvt_dd2f_avx_256(__m256d d0, __m256d d1)
{
    __m128 f0 = _mm256_cvtpd_ps(d0);
    __m128 f1 = _mm256_cvtpd_ps(d1);
    return _mm256_insertf128_ps(_mm256_castps128_ps256(f0), f1, 0x1);
}

/* SIMD4 reduce helper */
static gmx_inline float gmx_simdcall
gmx_simd4_reduce_f_avx_256(__m128 a)
{
    float f;
    a = _mm_hadd_ps(a, a);
    a = _mm_hadd_ps(a, a);
    _mm_store_ss(&f, a);
    return f;
}

/* SIMD4 Dotproduct helper function */
static gmx_inline float gmx_simdcall
gmx_simd4_dotproduct3_f_avx_256(__m128 a, __m128 b)
{
    float  f;
    __m128 c;
    a = _mm_mul_ps(a, b);
    c = _mm_add_ps(a, _mm_permute_ps(a, _MM_SHUFFLE(0, 3, 2, 1)));
    c = _mm_add_ps(c, _mm_permute_ps(a, _MM_SHUFFLE(1, 0, 3, 2)));
    _mm_store_ss(&f, c);
    return f;
}

static gmx_inline double gmx_simdcall
gmx_simd4_dotproduct3_d_avx_256(__m256d a, __m256d b)
{
    double  d;
    __m128d tmp1, tmp2;
    a    = _mm256_mul_pd(a, b);
    tmp1 = _mm256_castpd256_pd128(a);
    tmp2 = _mm256_extractf128_pd(a, 0x1);

    tmp1 = _mm_add_pd(tmp1, _mm_permute_pd(tmp1, _MM_SHUFFLE2(0, 1)));
    tmp1 = _mm_add_pd(tmp1, tmp2);
    _mm_store_sd(&d, tmp1);
    return d;
}


#endif /* GMX_SIMD_IMPL_X86_AVX_256_H */
