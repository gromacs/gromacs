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

#include <assert.h>
#include <math.h>

#include <immintrin.h>

#define GMX_SIMD_V2

#ifdef GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG
#    define gmx_mm_maskload_ps(mem, mask)        _mm_maskload_ps((mem), _mm_castsi128_ps(mask))
#    define gmx_mm_maskstore_ps(mem, mask, x)    _mm_maskstore_ps((mem), _mm_castsi128_ps(mask), (x))
#else
#    define gmx_mm_maskload_ps(mem, mask)        _mm_maskload_ps((mem), (mask))
#    define gmx_mm_maskstore_ps(mem, mask, x)    _mm_maskstore_ps((mem), (mask), (x))
#endif

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
#define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT
#define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_DOUBLE
#define GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT
#undef  GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE /* No need for half-simd, width is 4 */
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
#define gmx_simd_mul_mask_f(a, b, m)         _mm256_and_ps(_mm256_mul_ps(a, b), m)
#define gmx_simd_fmadd_mask_f(a, b, c, m)    _mm256_and_ps(gmx_simd_fmadd_f(a, b, c), m)
#ifdef NDEBUG
#    define gmx_simd_rcp_mask_f(a, m)        _mm256_and_ps(_mm256_rcp_ps(a), m)
#    define gmx_simd_rsqrt_mask_f(a, m)      _mm256_and_ps(_mm256_rsqrt_ps(a), m)
#else
/* For masked rcp/rsqrt we need to make sure we do not use the masked-out arguments if FP exceptions are enabled */
#    define gmx_simd_rcp_mask_f(a, m)        _mm256_and_ps(_mm256_rcp_ps(_mm256_blendv_ps(_mm256_set1_ps(1.0f), a, m)), m)
#    define gmx_simd_rsqrt_mask_f(a, m)      _mm256_and_ps(_mm256_rsqrt_ps(_mm256_blendv_ps(_mm256_set1_ps(1.0f), a, m)), m)
#endif
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
#define gmx_simd_cmpeq_f(a, b)     _mm256_cmp_ps(a, b, _CMP_EQ_OQ)
#define gmx_simd_cmpnz_f(a)        _mm256_cmp_ps(a, _mm256_setzero_ps(), _CMP_NEQ_OQ)
#define gmx_simd_cmplt_f(a, b)     _mm256_cmp_ps(a, b, _CMP_LT_OQ)
#define gmx_simd_cmple_f(a, b)     _mm256_cmp_ps(a, b, _CMP_LE_OQ)
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
/* gmx_simd_blendnotzero_fi not supported */
/* gmx_simd_blendv_fi not supported       */
/* Conversions between different booleans */
#define gmx_simd_cvt_fb2fib        _mm256_castps_si256
#define gmx_simd_cvt_fib2fb        _mm256_castsi256_ps
/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_f             gmx_simd_gather_load_transpose_f_avx_256
#define gmx_simd_best_pair_alignment_f               gmx_simd_best_pair_alignment_f_avx_256
#define gmx_simd_gather_loadu_transpose_f            gmx_simd_gather_loadu_transpose_f_avx_256
#define gmx_simd_transpose_scatter_storeu_f          gmx_simd_transpose_scatter_storeu_f_avx_256
#define gmx_simd_transpose_scatter_incru_f           gmx_simd_transpose_scatter_incru_f_avx_256
#define gmx_simd_transpose_scatter_decru_f           gmx_simd_transpose_scatter_decru_f_avx_256
#define gmx_simd_expand_scalars_to_triplets_f        gmx_simd_expand_scalars_to_triplets_f_avx_256
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_avx_256
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_avx_256
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_avx_256
/* Half-simd-width utilities (only needed in single, where width is 8) */
#define gmx_simd_load_dual_hsimd_f                   gmx_simd_load_dual_hsimd_f_avx_256
#define gmx_simd_loaddup_hsimd_f                     gmx_simd_loaddup_hsimd_f_avx_256
#define gmx_simd_load1_dual_hsimd_f                  gmx_simd_load1_dual_hsimd_f_avx_256
#define gmx_simd_store_dual_hsimd_f                  gmx_simd_store_dual_hsimd_f_avx_256
#define gmx_simd_decr_hsimd_f                        gmx_simd_decr_hsimd_f_avx_256
#define gmx_simd_gather_load_transpose_hsimd_f       gmx_simd_gather_load_transpose_hsimd_f_avx_256
#define gmx_simd_reduce_incr_4_return_sum_hsimd_f    gmx_simd_reduce_incr_4_return_sum_hsimd_f_avx_256

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
#define gmx_simd_mul_mask_d(a, b, m)       _mm256_and_pd(_mm256_mul_pd(a, b), m)
#define gmx_simd_fmadd_mask_d(a, b, c, m)  _mm256_and_pd(gmx_simd_fmadd_d(a, b, c), m)
#ifdef NDEBUG
#    define gmx_simd_rcp_mask_d(a, m)        _mm256_and_pd(gmx_simd_rcp_d(a), m)
#    define gmx_simd_rsqrt_mask_d(a, m)      _mm256_and_pd(gmx_simd_rsqrt_d(a), m)
#else
/* For masked rcp/rsqrt we need to make sure we do not use the masked-out arguments if FP exceptions are enabled */
#    define gmx_simd_rcp_mask_d(a, m)        _mm256_and_pd(gmx_simd_rcp_d(_mm256_blendv_pd(_mm256_set1_pd(1.0), a, m)), m)
#    define gmx_simd_rsqrt_mask_d(a, m)      _mm256_and_pd(gmx_simd_rsqrt_d(_mm256_blendv_pd(_mm256_set1_pd(1.0), a, m)), m)
#endif
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
#define gmx_simd_store_di(m, x)    _mm_store_si128((__m128i *)(m), x)
#define gmx_simd_loadu_di(m)       _mm_loadu_si128((const __m128i *)(m))
#define gmx_simd_storeu_di(m, x)   _mm_storeu_si128((__m128i *)(m), x)
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
#define gmx_simd_cmpeq_d(a, b)     _mm256_cmp_pd(a, b, _CMP_EQ_OQ)
#define gmx_simd_cmpnz_d(a)        _mm256_cmp_pd(a, _mm256_setzero_pd(), _CMP_NEQ_OQ)
#define gmx_simd_cmplt_d(a, b)     _mm256_cmp_pd(a, b, _CMP_LT_OQ)
#define gmx_simd_cmple_d(a, b)     _mm256_cmp_pd(a, b, _CMP_LE_OQ)
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
/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_d             gmx_simd_gather_load_transpose_d_avx_256
#define gmx_simd_best_pair_alignment_d               gmx_simd_best_pair_alignment_d_avx_256
#define gmx_simd_gather_loadu_transpose_d            gmx_simd_gather_loadu_transpose_d_avx_256
#define gmx_simd_transpose_scatter_storeu_d          gmx_simd_transpose_scatter_storeu_d_avx_256
#define gmx_simd_transpose_scatter_incru_d           gmx_simd_transpose_scatter_incru_d_avx_256
#define gmx_simd_transpose_scatter_decru_d           gmx_simd_transpose_scatter_decru_d_avx_256
#define gmx_simd_expand_scalars_to_triplets_d        gmx_simd_expand_scalars_to_triplets_d_avx_256
#define gmx_simd_gather_load_bysimdint_transpose_d   gmx_simd_gather_load_bysimdint_transpose_d_avx_256
#define gmx_simd_gather_loadu_bysimdint_transpose_d  gmx_simd_gather_loadu_bysimdint_transpose_d_avx_256
#define gmx_simd_reduce_incr_4_return_sum_d          gmx_simd_reduce_incr_4_return_sum_d_avx_256

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
#define gmx_simd4_transpose_f      gmx_simd4_transpose_f_avx_256
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
#define gmx_simd4_transpose_d       gmx_simd4_transpose_d_avx_256
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

    __m128 t0;
    t0 = _mm_add_ps(_mm256_castps256_ps128(a), _mm256_extractf128_ps(a, 0x1));
    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    _MM_EXTRACT_FLOAT(f, t0, 0);
    return f;
}

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

/*
 * This is an internal helper macro used by the three functions storing, incrementing, or
 * decrementing data. Do NOT use it outside this file.
 *
 * Input v0: [x0 x1 x2 x3 x4 x5 x6 x7]
 * Input v1: [y0 y1 y2 y3 y4 y5 y6 y7]
 * Input v2: [z0 z1 z2 z3 z4 z5 z6 z7]
 * Input v3: Unused
 *
 * Output v0: [x0 y0 z0 -  x4 y4 z4 - ]
 * Output v1: [x1 y1 z1 -  x5 y5 z5 - ]
 * Output v2: [x2 y2 z2 -  x6 y6 z6 - ]
 * Output v3: [x3 y3 z3 -  x7 y7 z7 - ]
 *
 * Here, - means undefined. Note that such values will not be zero!
 */
#define GMX_MM_TRANSPOSE_3TO4X2_FLOAT(v0, v1, v2, v3)                  \
    {                                                                   \
        __m256 gmx_mm_t1 = _mm256_unpacklo_ps(v0, v1);                  \
        __m256 gmx_mm_t2 = _mm256_unpackhi_ps(v0, v1);                  \
        v0 = _mm256_shuffle_ps(gmx_mm_t1, v2, _MM_SHUFFLE(0, 0, 1, 0)); \
        v1 = _mm256_shuffle_ps(gmx_mm_t1, v2, _MM_SHUFFLE(0, 1, 3, 2)); \
        v3 = _mm256_shuffle_ps(gmx_mm_t2, v2, _MM_SHUFFLE(0, 3, 3, 2)); \
        v2 = _mm256_shuffle_ps(gmx_mm_t2, v2, _MM_SHUFFLE(0, 2, 1, 0)); \
    }

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_avx_256(const float *        base,
                                         const gmx_int32_t    offset[],
                                         gmx_simd_float_t    &v0,
                                         gmx_simd_float_t    &v1,
                                         gmx_simd_float_t    &v2,
                                         gmx_simd_float_t    &v3)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;

    assert((size_t)offset % 32 == 0);
    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);

    t1  = _mm_load_ps( base + align * offset[0] );
    t2  = _mm_load_ps( base + align * offset[1] );
    t3  = _mm_load_ps( base + align * offset[2] );
    t4  = _mm_load_ps( base + align * offset[3] );
    t5  = _mm_load_ps( base + align * offset[4] );
    t6  = _mm_load_ps( base + align * offset[5] );
    t7  = _mm_load_ps( base + align * offset[6] );
    t8  = _mm_load_ps( base + align * offset[7] );

    v0  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    v1  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    v2  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    v3  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(v0, v1);
    tB  = _mm256_unpacklo_ps(v2, v3);
    tC  = _mm256_unpackhi_ps(v0, v1);
    tD  = _mm256_unpackhi_ps(v2, v3);

    v0  = _mm256_shuffle_ps(tA, tB, _MM_SHUFFLE(1, 0, 1, 0));
    v1  = _mm256_shuffle_ps(tA, tB, _MM_SHUFFLE(3, 2, 3, 2));
    v2  = _mm256_shuffle_ps(tC, tD, _MM_SHUFFLE(1, 0, 1, 0));
    v3  = _mm256_shuffle_ps(tC, tD, _MM_SHUFFLE(3, 2, 3, 2));
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_avx_256(const float *        base,
                                         const gmx_int32_t    offset[],
                                         gmx_simd_float_t    &v0,
                                         gmx_simd_float_t    &v1)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;

    assert((size_t)offset % 32 == 0);
    assert((size_t)base % 8 == 0);
    assert(align % 2 == 0);

    t1  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[0] ) );
    t2  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[1] ) );
    t3  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[2] ) );
    t4  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[3] ) );
    t5  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[4] ) );
    t6  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[5] ) );
    t7  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[6] ) );
    t8  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[7] ) );

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(tA, tC);
    tB  = _mm256_unpacklo_ps(tB, tD);
    v0  = _mm256_unpacklo_ps(tA, tB);
    v1  = _mm256_unpackhi_ps(tA, tB);
}

static const int gmx_simd_best_pair_alignment_f_avx_256 = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_f_avx_256(const float *        base,
                                          const gmx_int32_t    offset[],
                                          gmx_simd_float_t    &v0,
                                          gmx_simd_float_t    &v1,
                                          gmx_simd_float_t    &v2)
{
    __m256  t1, t2, t3, t4, t5, t6, t7, t8;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    assert((size_t)offset % 32 == 0);

    if ( (align & 0x3) == 0)
    {
        /* With alignment 4 or better we can read a byte beyond triplets and use _mm_loadu_ps() */
        t1  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps( base + align * offset[0] )),
                                   _mm_loadu_ps( base + align * offset[4] ), 0x1);
        t2  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps(base + align * offset[1] )),
                                   _mm_loadu_ps( base + align * offset[5] ), 0x1);
        t3  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps(base + align * offset[2] )),
                                   _mm_loadu_ps( base + align * offset[6] ), 0x1);
        t4  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps(base + align * offset[3] )),
                                   _mm_loadu_ps( base + align * offset[7] ), 0x1);
    }
    else
    {
        /* Arbitrary alignment */
        t1  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[0], mask)),
                                   gmx_mm_maskload_ps( base + align * offset[4], mask), 0x1);
        t2  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[1], mask)),
                                   gmx_mm_maskload_ps( base + align * offset[5], mask), 0x1);
        t3  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[2], mask)),
                                   gmx_mm_maskload_ps( base + align * offset[6], mask), 0x1);
        t4  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[3], mask)),
                                   gmx_mm_maskload_ps( base + align * offset[7], mask), 0x1);
    }
    t5  = _mm256_unpacklo_ps(t1, t2);
    t6  = _mm256_unpacklo_ps(t3, t4);
    t7  = _mm256_unpackhi_ps(t1, t2);
    t8  = _mm256_unpackhi_ps(t3, t4);
    v0  = _mm256_shuffle_ps(t5, t6, _MM_SHUFFLE(1, 0, 1, 0));
    v1  = _mm256_shuffle_ps(t5, t6, _MM_SHUFFLE(3, 2, 3, 2));
    v2  = _mm256_shuffle_ps(t7, t8, _MM_SHUFFLE(1, 0, 1, 0));
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_f_avx_256(float *              base,
                                            const gmx_int32_t    offset[],
                                            gmx_simd_float_t     v0,
                                            gmx_simd_float_t     v1,
                                            gmx_simd_float_t     v2)
{
    __m256  v3;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    assert((size_t)offset % 32 == 0);

    GMX_MM_TRANSPOSE_3TO4X2_FLOAT(v0, v1, v2, v3);
    gmx_mm_maskstore_ps( base + align * offset[0], mask, _mm256_castps256_ps128(v0));
    gmx_mm_maskstore_ps( base + align * offset[1], mask, _mm256_castps256_ps128(v1));
    gmx_mm_maskstore_ps( base + align * offset[2], mask, _mm256_castps256_ps128(v2));
    gmx_mm_maskstore_ps( base + align * offset[3], mask, _mm256_castps256_ps128(v3));
    gmx_mm_maskstore_ps( base + align * offset[4], mask, _mm256_extractf128_ps(v0, 0x1));
    gmx_mm_maskstore_ps( base + align * offset[5], mask, _mm256_extractf128_ps(v1, 0x1));
    gmx_mm_maskstore_ps( base + align * offset[6], mask, _mm256_extractf128_ps(v2, 0x1));
    gmx_mm_maskstore_ps( base + align * offset[7], mask, _mm256_extractf128_ps(v3, 0x1));
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_f_avx_256(float *              base,
                                           const gmx_int32_t    offset[],
                                           gmx_simd_float_t     v0,
                                           gmx_simd_float_t     v1,
                                           gmx_simd_float_t     v2)
{
    __m256  v3, t0, t1, t2, t3;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    assert((size_t)offset % 32 == 0);

    t0  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[0], mask)),
                               gmx_mm_maskload_ps( base + align * offset[4], mask), 0x1);
    t1  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[1], mask)),
                               gmx_mm_maskload_ps( base + align * offset[5], mask), 0x1);
    t2  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[2], mask)),
                               gmx_mm_maskload_ps( base + align * offset[6], mask), 0x1);
    t3  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[3], mask)),
                               gmx_mm_maskload_ps( base + align * offset[7], mask), 0x1);

    GMX_MM_TRANSPOSE_3TO4X2_FLOAT(v0, v1, v2, v3);
    t0          = _mm256_add_ps(t0, v0);
    t1          = _mm256_add_ps(t1, v1);
    t2          = _mm256_add_ps(t2, v2);
    t3          = _mm256_add_ps(t3, v3);

    gmx_mm_maskstore_ps(base + align * offset[0], mask, _mm256_castps256_ps128(t0));
    gmx_mm_maskstore_ps(base + align * offset[1], mask, _mm256_castps256_ps128(t1));
    gmx_mm_maskstore_ps(base + align * offset[2], mask, _mm256_castps256_ps128(t2));
    gmx_mm_maskstore_ps(base + align * offset[3], mask, _mm256_castps256_ps128(t3));
    gmx_mm_maskstore_ps(base + align * offset[4], mask, _mm256_extractf128_ps(t0, 0x1));
    gmx_mm_maskstore_ps(base + align * offset[5], mask, _mm256_extractf128_ps(t1, 0x1));
    gmx_mm_maskstore_ps(base + align * offset[6], mask, _mm256_extractf128_ps(t2, 0x1));
    gmx_mm_maskstore_ps(base + align * offset[7], mask, _mm256_extractf128_ps(t3, 0x1));
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_f_avx_256(float *              base,
                                           const gmx_int32_t    offset[],
                                           gmx_simd_float_t     v0,
                                           gmx_simd_float_t     v1,
                                           gmx_simd_float_t     v2)
{
    __m256  v3, t0, t1, t2, t3;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    assert((size_t)offset % 32 == 0);

    t0  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[0], mask)),
                               gmx_mm_maskload_ps( base + align * offset[4], mask), 0x1);
    t1  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[1], mask)),
                               gmx_mm_maskload_ps( base + align * offset[5], mask), 0x1);
    t2  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[2], mask)),
                               gmx_mm_maskload_ps( base + align * offset[6], mask), 0x1);
    t3  = _mm256_insertf128_ps(_mm256_castps128_ps256(gmx_mm_maskload_ps( base + align * offset[3], mask)),
                               gmx_mm_maskload_ps( base + align * offset[7], mask), 0x1);

    GMX_MM_TRANSPOSE_3TO4X2_FLOAT(v0, v1, v2, v3);
    t0          = _mm256_sub_ps(t0, v0);
    t1          = _mm256_sub_ps(t1, v1);
    t2          = _mm256_sub_ps(t2, v2);
    t3          = _mm256_sub_ps(t3, v3);

    gmx_mm_maskstore_ps(base + align * offset[0], mask, _mm256_castps256_ps128(t0));
    gmx_mm_maskstore_ps(base + align * offset[1], mask, _mm256_castps256_ps128(t1));
    gmx_mm_maskstore_ps(base + align * offset[2], mask, _mm256_castps256_ps128(t2));
    gmx_mm_maskstore_ps(base + align * offset[3], mask, _mm256_castps256_ps128(t3));
    gmx_mm_maskstore_ps(base + align * offset[4], mask, _mm256_extractf128_ps(t0, 0x1));
    gmx_mm_maskstore_ps(base + align * offset[5], mask, _mm256_extractf128_ps(t1, 0x1));
    gmx_mm_maskstore_ps(base + align * offset[6], mask, _mm256_extractf128_ps(t2, 0x1));
    gmx_mm_maskstore_ps(base + align * offset[7], mask, _mm256_extractf128_ps(t3, 0x1));
}

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_f_avx_256(gmx_simd_float_t   scalar,
                                              gmx_simd_float_t  &triplets0,
                                              gmx_simd_float_t  &triplets1,
                                              gmx_simd_float_t  &triplets2)
{
    __m256 t0 = _mm256_permute2f128_ps(scalar, scalar, 0x21);
    __m256 t1 = _mm256_permute_ps(scalar, _MM_SHUFFLE(1, 0, 0, 0));
    __m256 t2 = _mm256_permute_ps(t0, _MM_SHUFFLE(2, 2, 1, 1));
    __m256 t3 = _mm256_permute_ps(scalar, _MM_SHUFFLE(3, 3, 3, 2));
    triplets0 = _mm256_blend_ps(t1, t2, 0xF0);
    triplets1 = _mm256_blend_ps(t3, t1, 0xF0);
    triplets2 = _mm256_blend_ps(t2, t3, 0xF0);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_avx_256(const float *       base,
                                                   gmx_simd_fint32_t   simdoffset,
                                                   gmx_simd_float_t   &v0,
                                                   gmx_simd_float_t   &v1,
                                                   gmx_simd_float_t   &v2,
                                                   gmx_simd_float_t   &v3)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;
    /* All x86 compilers we know that support SIMD also support GCC alignment attributes */
    int __attribute__ ((aligned (32))) offset[GMX_SIMD_FLOAT_WIDTH];

    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);

    /* Plain AVX does not support 256-bit integer shift operations,
     * so we always store to memory and multiply in the integer units instead.
     */
    _mm256_store_si256( (__m256i *)offset, simdoffset);

    t1  = _mm_load_ps( base + align * offset[0] );
    t2  = _mm_load_ps( base + align * offset[1] );
    t3  = _mm_load_ps( base + align * offset[2] );
    t4  = _mm_load_ps( base + align * offset[3] );
    t5  = _mm_load_ps( base + align * offset[4] );
    t6  = _mm_load_ps( base + align * offset[5] );
    t7  = _mm_load_ps( base + align * offset[6] );
    t8  = _mm_load_ps( base + align * offset[7] );

    v0  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    v1  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    v2  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    v3  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(v0, v1);
    tB  = _mm256_unpacklo_ps(v2, v3);
    tC  = _mm256_unpackhi_ps(v0, v1);
    tD  = _mm256_unpackhi_ps(v2, v3);

    v0  = _mm256_shuffle_ps(tA, tB, _MM_SHUFFLE(1, 0, 1, 0));
    v1  = _mm256_shuffle_ps(tA, tB, _MM_SHUFFLE(3, 2, 3, 2));
    v2  = _mm256_shuffle_ps(tC, tD, _MM_SHUFFLE(1, 0, 1, 0));
    v3  = _mm256_shuffle_ps(tC, tD, _MM_SHUFFLE(3, 2, 3, 2));
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_avx_256(const float *       base,
                                                   gmx_simd_fint32_t   simdoffset,
                                                   gmx_simd_float_t   &v0,
                                                   gmx_simd_float_t   &v1)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;
    /* All x86 compilers we know that support SIMD also support GCC alignment attributes */
    int __attribute__ ((aligned (32))) offset[GMX_SIMD_FLOAT_WIDTH];

    assert((size_t)base % 8 == 0);
    assert(align % 2 == 0);

    /* Plain AVX does not support 256-bit integer shift operations,
     * so we always store to memory and multiply in the integer units instead.
     */
    _mm256_store_si256( (__m256i *)offset, simdoffset);

    t1  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[0] ) );
    t2  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[1] ) );
    t3  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[2] ) );
    t4  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[3] ) );
    t5  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[4] ) );
    t6  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[5] ) );
    t7  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[6] ) );
    t8  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[7] ) );

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(tA, tC);
    tB  = _mm256_unpacklo_ps(tB, tD);
    v0  = _mm256_unpacklo_ps(tA, tB);
    v1  = _mm256_unpackhi_ps(tA, tB);
}

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_f_avx_256(const float *       base,
                                                    gmx_simd_fint32_t   simdoffset,
                                                    gmx_simd_float_t   &v0,
                                                    gmx_simd_float_t   &v1)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    __m256 tA, tB, tC, tD;
    /* All x86 compilers we know that support SIMD also support GCC alignment attributes */
    int __attribute__ ((aligned (32))) offset[GMX_SIMD_FLOAT_WIDTH];

    /* Plain AVX does not support 256-bit integer shift operations,
     * so we always store to memory and multiply in the integer units instead.
     */
    _mm256_store_si256( (__m256i *)offset, simdoffset);

    t1  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[0] ) );
    t2  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[1] ) );
    t3  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[2] ) );
    t4  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[3] ) );
    t5  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[4] ) );
    t6  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[5] ) );
    t7  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[6] ) );
    t8  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)( base + align * offset[7] ) );

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t4), t8, 0x1);

    tA  = _mm256_unpacklo_ps(tA, tC);
    tB  = _mm256_unpacklo_ps(tB, tD);
    v0  = _mm256_unpacklo_ps(tA, tB);
    v1  = _mm256_unpackhi_ps(tA, tB);
}

static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_f_avx_256(float *           m,
                                            gmx_simd_float_t  v0,
                                            gmx_simd_float_t  v1,
                                            gmx_simd_float_t  v2,
                                            gmx_simd_float_t  v3)
{
    float  f;
    __m128 t0, t2;

    assert((size_t)m % 16 == 0);

    v0 = _mm256_hadd_ps(v0, v1);
    v2 = _mm256_hadd_ps(v2, v3);
    v0 = _mm256_hadd_ps(v0, v2);
    t0 = _mm_add_ps(_mm256_castps256_ps128(v0), _mm256_extractf128_ps(v0, 0x1));

    t2 = _mm_add_ps(t0, _mm_load_ps(m));
    _mm_store_ps(m, t2);

    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    _MM_EXTRACT_FLOAT(f, t0, 0);
    return f;
}

/******************************************************
 * Single precision half-simd-width utility functions *
 ******************************************************/
static gmx_inline gmx_simd_float_t
gmx_simd_load_dual_hsimd_f_avx_256(const float * m0,
                                   const float * m1)
{
    assert((size_t)m0 % 16 == 0);
    assert((size_t)m1 % 16 == 0);
    return _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_load_ps(m0)), _mm_load_ps(m1), 0x1);
}

static gmx_inline gmx_simd_float_t
gmx_simd_loaddup_hsimd_f_avx_256(const float * m)
{
    assert((size_t)m % 16 == 0);
    return _mm256_broadcast_ps((const __m128 *)m);
}

static gmx_inline gmx_simd_float_t
gmx_simd_load1_dual_hsimd_f_avx_256(const float * m)
{
    __m128 t0, t1;
    t0 = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)m);
    t1 = _mm_permute_ps(t0, _MM_SHUFFLE(1, 1, 1, 1));
    t0 = _mm_permute_ps(t0, _MM_SHUFFLE(0, 0, 0, 0));
    return _mm256_insertf128_ps(_mm256_castps128_ps256(t0), t1, 0x1);
}


static gmx_inline void
gmx_simd_store_dual_hsimd_f_avx_256(float *           m0,
                                    float *           m1,
                                    gmx_simd_float_t  a)
{
    assert((size_t)m0 % 16 == 0);
    assert((size_t)m1 % 16 == 0);
    _mm_store_ps(m0, _mm256_castps256_ps128(a));
    _mm_store_ps(m1, _mm256_extractf128_ps(a, 0x1));
}

static gmx_inline void
gmx_simd_decr_hsimd_f_avx_256(float *          m,
                              gmx_simd_float_t a)
{
    assert((size_t)m % 16 == 0);
    __m128 asum = _mm_add_ps(_mm256_castps256_ps128(a), _mm256_extractf128_ps(a, 0x1));
    _mm_store_ps(m, _mm_sub_ps(_mm_load_ps(m), asum));
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_hsimd_f_avx_256(const float *       base0,
                                               const float *       base1,
                                               gmx_int32_t         offset[],
                                               gmx_simd_float_t   &v0,
                                               gmx_simd_float_t   &v1)
{
    __m128 t0, t1, t2, t3, t4, t5, t6, t7;
    __m256 tA, tB, tC, tD;

    assert((size_t)offset % 16 == 0);
    assert((size_t)base0 % 8 == 0);
    assert((size_t)base1 % 8 == 0);
    assert(align % 2 == 0);

    t0  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base0 + align * offset[0]));
    t1  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base0 + align * offset[1]));
    t2  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base0 + align * offset[2]));
    t3  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base0 + align * offset[3]));
    t4  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base1 + align * offset[0]));
    t5  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base1 + align * offset[1]));
    t6  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base1 + align * offset[2]));
    t7  = _mm_loadl_pi(_mm_setzero_ps(), (__m64 *)(base1 + align * offset[3]));

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t0), t4, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);

    tA  = _mm256_unpacklo_ps(tA, tC);
    tB  = _mm256_unpacklo_ps(tB, tD);
    v0  = _mm256_unpacklo_ps(tA, tB);
    v1  = _mm256_unpackhi_ps(tA, tB);
}


static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_hsimd_f_avx_256(float *            m,
                                                  gmx_simd_float_t   v0,
                                                  gmx_simd_float_t   v1)
{
    __m128 t0, t1;
    float  f;

    assert((size_t)m % 16 == 0);

    v0  = _mm256_hadd_ps(v0, v1);
    t0  = _mm256_extractf128_ps(v0, 0x1);
    t0  = _mm_hadd_ps(_mm256_castps256_ps128(v0), t0);
    t0  = _mm_permute_ps(t0, _MM_SHUFFLE(3, 1, 2, 0));

    t1  = _mm_add_ps(t0, _mm_load_ps(m));
    _mm_store_ps(m, t1);

    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    _MM_EXTRACT_FLOAT(f, t0, 0);
    return f;
}
#endif

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
    a  = _mm256_add_pd(a, _mm256_permute_pd(a, 0b0101 ));
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

/****************************************************
 * Double precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

/* Internal macro: Full 4x4 transpose of __m256d */
#define GMX_MM_TRANSPOSE4_DOUBLE(v0, v1, v2, v3)                \
    {                                                            \
        __m256d gmx_mm_t1, gmx_mm_t2, gmx_mm_t3, gmx_mm_t4;      \
        gmx_mm_t1 = _mm256_unpacklo_pd(v0, v1);                  \
        gmx_mm_t2 = _mm256_unpackhi_pd(v0, v1);                  \
        gmx_mm_t3 = _mm256_unpacklo_pd(v2, v3);                  \
        gmx_mm_t4 = _mm256_unpackhi_pd(v2, v3);                  \
        v0        = _mm256_permute2f128_pd(gmx_mm_t1, gmx_mm_t3, 0x20); \
        v1        = _mm256_permute2f128_pd(gmx_mm_t2, gmx_mm_t4, 0x20); \
        v2        = _mm256_permute2f128_pd(gmx_mm_t1, gmx_mm_t3, 0x31); \
        v3        = _mm256_permute2f128_pd(gmx_mm_t2, gmx_mm_t4, 0x31); \
    }

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d_avx_256(const double *        base,
                                         const gmx_int32_t     offset[],
                                         gmx_simd_double_t    &v0,
                                         gmx_simd_double_t    &v1,
                                         gmx_simd_double_t    &v2,
                                         gmx_simd_double_t    &v3)
{
    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 32 == 0);
    assert(align % 4 == 0);

    v0  = _mm256_load_pd( base + align * offset[0] );
    v1  = _mm256_load_pd( base + align * offset[1] );
    v2  = _mm256_load_pd( base + align * offset[2] );
    v3  = _mm256_load_pd( base + align * offset[3] );
    GMX_MM_TRANSPOSE4_DOUBLE(v0, v1, v2, v3);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d_avx_256(const double *        base,
                                         const gmx_int32_t     offset[],
                                         gmx_simd_double_t    &v0,
                                         gmx_simd_double_t    &v1)
{
    __m128d t1, t2, t3, t4;
    __m256d tA, tB;

    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 16 == 0);
    assert(align % 2 == 0);

    t1   = _mm_load_pd( base + align * offset[0] );
    t2   = _mm_load_pd( base + align * offset[1] );
    t3   = _mm_load_pd( base + align * offset[2] );
    t4   = _mm_load_pd( base + align * offset[3] );
    tA   = _mm256_insertf128_pd(_mm256_castpd128_pd256(t1), t3, 0x1);
    tB   = _mm256_insertf128_pd(_mm256_castpd128_pd256(t2), t4, 0x1);

    v0  = _mm256_unpacklo_pd(tA, tB);
    v1  = _mm256_unpackhi_pd(tA, tB);
}

static const int gmx_simd_best_pair_alignment_d_avx_256 = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_d_avx_256(const double *        base,
                                          const gmx_int32_t     offset[],
                                          gmx_simd_double_t    &v0,
                                          gmx_simd_double_t    &v1,
                                          gmx_simd_double_t    &v2)
{
    assert((size_t)offset % 16 == 0);

    __m256i mask = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_set_epi32(-1, -1, -1, -1)), _mm_set_epi32(0, 0, -1, -1), 0x1);
    __m256d t1, t2, t3, t4, t5, t6, t7, t8;
    t1  = _mm256_maskload_pd(base + align * offset[0], mask);
    t2  = _mm256_maskload_pd(base + align * offset[1], mask);
    t3  = _mm256_maskload_pd(base + align * offset[2], mask);
    t4  = _mm256_maskload_pd(base + align * offset[3], mask);
    t5  = _mm256_unpacklo_pd(t1, t2);
    t6  = _mm256_unpackhi_pd(t1, t2);
    t7  = _mm256_unpacklo_pd(t3, t4);
    t8  = _mm256_unpackhi_pd(t3, t4);
    v0  = _mm256_permute2f128_pd(t5, t7, 0x20);
    v1  = _mm256_permute2f128_pd(t6, t8, 0x20);
    v2  = _mm256_permute2f128_pd(t5, t7, 0x31);
}



template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_d_avx_256(double *              base,
                                            const gmx_int32_t     offset[],
                                            gmx_simd_double_t     v0,
                                            gmx_simd_double_t     v1,
                                            gmx_simd_double_t     v2)
{
    __m256i mask = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_set_epi32(-1, -1, -1, -1)), _mm_set_epi32(0, 0, -1, -1), 0x1);
    __m256d v3   = _mm256_setzero_pd();

    assert((size_t)offset % 16 == 0);

    GMX_MM_TRANSPOSE4_DOUBLE(v0, v1, v2, v3);

    _mm256_maskstore_pd(base + align * offset[0], mask, v0);
    _mm256_maskstore_pd(base + align * offset[1], mask, v1);
    _mm256_maskstore_pd(base + align * offset[2], mask, v2);
    _mm256_maskstore_pd(base + align * offset[3], mask, v3);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_d_avx_256(double *              base,
                                           const gmx_int32_t     offset[],
                                           gmx_simd_double_t     v0,
                                           gmx_simd_double_t     v1,
                                           gmx_simd_double_t     v2)
{
    __m256i mask = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_set_epi32(-1, -1, -1, -1)), _mm_set_epi32(0, 0, -1, -1), 0x1);
    __m256d v3   = _mm256_setzero_pd();
    __m256d t0, t1, t2, t3;

    assert((size_t)offset % 16 == 0);

    t0  = _mm256_maskload_pd(base + align * offset[0], mask);
    t1  = _mm256_maskload_pd(base + align * offset[1], mask);
    t2  = _mm256_maskload_pd(base + align * offset[2], mask);
    t3  = _mm256_maskload_pd(base + align * offset[3], mask);

    GMX_MM_TRANSPOSE4_DOUBLE(v0, v1, v2, v3);

    t0  = _mm256_add_pd(t0, v0);
    t1  = _mm256_add_pd(t1, v1);
    t2  = _mm256_add_pd(t2, v2);
    t3  = _mm256_add_pd(t3, v3);

    _mm256_maskstore_pd(base + align * offset[0], mask, t0);
    _mm256_maskstore_pd(base + align * offset[1], mask, t1);
    _mm256_maskstore_pd(base + align * offset[2], mask, t2);
    _mm256_maskstore_pd(base + align * offset[3], mask, t3);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_d_avx_256(double *              base,
                                           const gmx_int32_t     offset[],
                                           gmx_simd_double_t     v0,
                                           gmx_simd_double_t     v1,
                                           gmx_simd_double_t     v2)
{
    __m256i mask = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_set_epi32(-1, -1, -1, -1)), _mm_set_epi32(0, 0, -1, -1), 0x1);
    __m256d v3   = _mm256_setzero_pd();
    __m256d t0, t1, t2, t3;

    assert((size_t)offset % 16 == 0);

    t0  = _mm256_maskload_pd(base + align * offset[0], mask);
    t1  = _mm256_maskload_pd(base + align * offset[1], mask);
    t2  = _mm256_maskload_pd(base + align * offset[2], mask);
    t3  = _mm256_maskload_pd(base + align * offset[3], mask);

    GMX_MM_TRANSPOSE4_DOUBLE(v0, v1, v2, v3);

    t0  = _mm256_sub_pd(t0, v0);
    t1  = _mm256_sub_pd(t1, v1);
    t2  = _mm256_sub_pd(t2, v2);
    t3  = _mm256_sub_pd(t3, v3);

    _mm256_maskstore_pd(base + align * offset[0], mask, t0);
    _mm256_maskstore_pd(base + align * offset[1], mask, t1);
    _mm256_maskstore_pd(base + align * offset[2], mask, t2);
    _mm256_maskstore_pd(base + align * offset[3], mask, t3);
}

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_d_avx_256(gmx_simd_double_t   scalar,
                                              gmx_simd_double_t  &triplets0,
                                              gmx_simd_double_t  &triplets1,
                                              gmx_simd_double_t  &triplets2)
{
    __m256d t0 = _mm256_permute2f128_pd(scalar, scalar, 0x21);
    __m256d t1 = _mm256_permute_pd(scalar, 0b0000);
    __m256d t2 = _mm256_permute_pd(scalar, 0b1111);
    triplets0 = _mm256_blend_pd(t1, t0, 0b1100);
    triplets1 = _mm256_blend_pd(t2, t1, 0b1100);
    triplets2 = _mm256_blend_pd(t0, t2, 0b1100);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_avx_256(const double *       base,
                                                   gmx_simd_dint32_t    offset,
                                                   gmx_simd_double_t   &v0,
                                                   gmx_simd_double_t   &v1,
                                                   gmx_simd_double_t   &v2,
                                                   gmx_simd_double_t   &v3)
{
    assert((size_t)base % 32 == 0);
    assert(align % 4 == 0);

    /* Use optimized bit-shift multiply for the most common alignments */
    if (align == 4)
    {
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 8)
    {
        offset = _mm_slli_epi32(offset, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset = _mm_add_epi32(offset, _mm_slli_epi32(offset, 1));
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 16)
    {
        offset = _mm_slli_epi32(offset, 4);
    }

    if (align == 4 || align == 8 || align == 12 || align == 16)
    {
        v0  = _mm256_load_pd(base + _mm_extract_epi32(offset, 0));
        v1  = _mm256_load_pd(base + _mm_extract_epi32(offset, 1));
        v2  = _mm256_load_pd(base + _mm_extract_epi32(offset, 2));
        v3  = _mm256_load_pd(base + _mm_extract_epi32(offset, 3));
    }
    else
    {
        v0  = _mm256_load_pd(base + align * _mm_extract_epi32(offset, 0));
        v1  = _mm256_load_pd(base + align * _mm_extract_epi32(offset, 1));
        v2  = _mm256_load_pd(base + align * _mm_extract_epi32(offset, 2));
        v3  = _mm256_load_pd(base + align * _mm_extract_epi32(offset, 3));
    }
    GMX_MM_TRANSPOSE4_DOUBLE(v0, v1, v2, v3);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_avx_256(const double *       base,
                                                   gmx_simd_dint32_t    offset,
                                                   gmx_simd_double_t   &v0,
                                                   gmx_simd_double_t   &v1)
{
    __m128d t1, t2, t3, t4;
    __m256d tA, tB;

    assert((size_t)base % 16 == 0);
    assert(align % 2 == 0);

    /* Use optimized bit-shift multiply for the most common alignments */
    if (align == 2)
    {
        offset = _mm_slli_epi32(offset, 1);
    }
    else if (align == 4)
    {
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 6)
    {
        /* multiply by 3, then by 2 */
        offset = _mm_add_epi32(offset, _mm_slli_epi32(offset, 1));
        offset = _mm_slli_epi32(offset, 1);
    }
    else if (align == 8)
    {
        offset = _mm_slli_epi32(offset, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset = _mm_add_epi32(offset, _mm_slli_epi32(offset, 1));
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 16)
    {
        offset = _mm_slli_epi32(offset, 4);
    }

    if (align == 2 || align == 4 || align == 6 ||
        align == 8 || align == 12 || align == 16)
    {
        t1  = _mm_load_pd(base + _mm_extract_epi32(offset, 0));
        t2  = _mm_load_pd(base + _mm_extract_epi32(offset, 1));
        t3  = _mm_load_pd(base + _mm_extract_epi32(offset, 2));
        t4  = _mm_load_pd(base + _mm_extract_epi32(offset, 3));
    }
    else
    {
        t1  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 0));
        t2  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 1));
        t3  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 2));
        t4  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 3));
    }

    tA  = _mm256_insertf128_pd(_mm256_castpd128_pd256(t1), t3, 0x1);
    tB  = _mm256_insertf128_pd(_mm256_castpd128_pd256(t2), t4, 0x1);
    v0  = _mm256_unpacklo_pd(tA, tB);
    v1  = _mm256_unpackhi_pd(tA, tB);
}


template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_d_avx_256(const double *       base,
                                                    gmx_simd_dint32_t    offset,
                                                    gmx_simd_double_t   &v0,
                                                    gmx_simd_double_t   &v1)
{
    __m128d t1, t2, t3, t4;
    __m256d tA, tB;

    /* Use optimized bit-shift multiply for the most common alignments. */
    /* Do nothing for align == 1 */
    if (align == 2)
    {
        offset = _mm_slli_epi32(offset, 1);
    }
    else if (align == 4)
    {
        offset = _mm_slli_epi32(offset, 2);
    }

    if (align == 1 || align == 2 || align == 4)
    {
        t1   = _mm_loadu_pd(base + _mm_extract_epi32(offset, 0));
        t2   = _mm_loadu_pd(base + _mm_extract_epi32(offset, 1));
        t3   = _mm_loadu_pd(base + _mm_extract_epi32(offset, 2));
        t4   = _mm_loadu_pd(base + _mm_extract_epi32(offset, 3));
    }
    else
    {
        t1   = _mm_loadu_pd(base + align * _mm_extract_epi32(offset, 0));
        t2   = _mm_loadu_pd(base + align * _mm_extract_epi32(offset, 1));
        t3   = _mm_loadu_pd(base + align * _mm_extract_epi32(offset, 2));
        t4   = _mm_loadu_pd(base + align * _mm_extract_epi32(offset, 3));
    }
    tA  = _mm256_insertf128_pd(_mm256_castpd128_pd256(t1), t3, 0x1);
    tB  = _mm256_insertf128_pd(_mm256_castpd128_pd256(t2), t4, 0x1);

    v0  = _mm256_unpacklo_pd(tA, tB);
    v1  = _mm256_unpackhi_pd(tA, tB);
}


static gmx_inline double gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_d_avx_256(double * m,
                                            gmx_simd_double_t v0, gmx_simd_double_t v1,
                                            gmx_simd_double_t v2, gmx_simd_double_t v3)
{
    __m256d t0, t1, t2;

    assert((size_t)m % 32 == 0);

    t0 = _mm256_hadd_pd(v0, v1);
    t1 = _mm256_hadd_pd(v2, v3);
    t2 = _mm256_permute2f128_pd(t0, t1, 0x21);
    t0 = _mm256_add_pd(t0, t2);
    t1 = _mm256_add_pd(t1, t2);
    t0 = _mm256_blend_pd(t0, t1, 0b1100);

    t1 = _mm256_add_pd(t0, _mm256_load_pd(m));
    _mm256_store_pd(m, t1);

    return gmx_simd_reduce_d_avx_256(t0);
}
#endif


/**************************
 * SIMD4 helper functions *
 **************************/

/* SIMD4 reduce helper */
static gmx_inline float gmx_simdcall
gmx_simd4_reduce_f_avx_256(__m128 a)
{
    float f;
    a = _mm_add_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(1, 0, 3, 2)));
    a = _mm_add_ss(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(0, 3, 2, 1)));
    _MM_EXTRACT_FLOAT(f, a, 0);
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

#ifdef __cplusplus
static gmx_inline void gmx_simdcall
gmx_simd4_transpose_f_avx_256(gmx_simd4_float_t &v0, gmx_simd4_float_t &v1,
                              gmx_simd4_float_t &v2, gmx_simd4_float_t &v3)
{
    _MM_TRANSPOSE4_PS(v0, v1, v2, v3);
}
#endif

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

#ifdef __cplusplus
static gmx_inline void gmx_simdcall
gmx_simd4_transpose_d_avx_256(gmx_simd4_double_t &v0, gmx_simd4_double_t &v1,
                              gmx_simd4_double_t &v2, gmx_simd4_double_t &v3)
{
    GMX_MM_TRANSPOSE4_DOUBLE(v0, v1, v2, v3);
}
#endif


/* Function to check whether SIMD operations have resulted in overflow */
static int
gmx_simd_check_and_reset_overflow(void)
{
    int MXCSR;
    int sse_overflow;

    MXCSR = _mm_getcsr();
    /* The overflow flag is bit 3 in the register */
    if (MXCSR & 0x0008)
    {
        sse_overflow = 1;
        /* Set the overflow flag to zero */
        MXCSR = MXCSR & 0xFFF7;
        _mm_setcsr(MXCSR);
    }
    else
    {
        sse_overflow = 0;
    }
    return sse_overflow;
}


#endif /* GMX_SIMD_IMPL_X86_AVX_256_H */
