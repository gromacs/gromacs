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

#ifndef GMX_SIMD_IMPL_X86_AVX_512F_H
#define GMX_SIMD_IMPL_X86_AVX_512F_H

#include <math.h>

#include <immintrin.h>

/* Intel AVX-512F */

/* This implementation works on the Intel compiler (version 14.0.2), but gcc-4.9
 * is still missing the cast operations needed to reinterpret between
 * 512-bit SIMD variables.
 * In addition, the truncate intrinsic is not yet available on gcc, and if
 * we try to implement it through _mm512_roundscale_round_ps() we get strange
 * errors about invalid immediate values when compiling Gromacs, but not a
 * simple test program. Hopefully these glitches will have been fixed in gcc
 * by the time AVX-512F-capable hardware appears.
 *
 * A couple of additional notes:
 *
 * - We avoid using the reduce calls provided by the intel compiler, since they
 *   are not instructions but software routines, and not implemented on gcc.
 * - gcc-4.9 does not yet provide int2mask and mask2int. This will be trivial
 *   to accomplish with plain casts, though.
 * - gcc-4.9 can do conversions between float/integer/double with simple casts
 *   like (__m512) as long as the variables are of the same size. However, it
 *   does not work to cast e.g. __m512 to __m128 this way. Since the sharing
 *   of xmm/ymm/zmm registers is more or less the whole idea with AVX-512
 *   (compared to MIC), this will have to be fixed in gcc.
 */

/* Capability definitions for AVX-512 SIMD. */
#define GMX_SIMD_HAVE_FLOAT
#define GMX_SIMD_HAVE_DOUBLE
#define GMX_SIMD_HAVE_SIMD_HARDWARE
#define GMX_SIMD_HAVE_LOADU
#define GMX_SIMD_HAVE_STOREU
#define GMX_SIMD_HAVE_LOGICAL
#define GMX_SIMD_HAVE_FMA
#undef  GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD_HAVE_FINT32
/* Technically it is straightforward to emulate extract on AVX-512F through
 * memory operations, but when applied to 16 elements as part of a table lookup
 * it will be faster to just store the entire vector once, so we avoid setting it.
 */
#undef  GMX_SIMD_HAVE_FINT32_EXTRACT
#define GMX_SIMD_HAVE_FINT32_LOGICAL
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS
#define GMX_SIMD_HAVE_DINT32
#undef  GMX_SIMD_HAVE_DINT32_EXTRACT
#define GMX_SIMD_HAVE_DINT32_LOGICAL
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS
#define GMX_SIMD4_HAVE_FLOAT
#define GMX_SIMD4_HAVE_DOUBLE

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH        16
#define GMX_SIMD_DOUBLE_WIDTH        8
#define GMX_SIMD_FINT32_WIDTH       16
#define GMX_SIMD_DINT32_WIDTH        8
#define GMX_SIMD_RSQRT_BITS         14
#define GMX_SIMD_RCP_BITS           14

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
#define gmx_simd_fneg_f(x)         gmx_simd_xor_f(x, _mm512_set1_ps(-0.0f))
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




/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_double_t          __m512d
#define gmx_simd_load_d            _mm512_load_pd
/* Avoid using _mm512_extload_pd() since it is not available on gcc-4.9 */
#define gmx_simd_load1_d(m)        _mm512_broadcastsd_pd(_mm_load1_pd(m))
#define gmx_simd_set1_d            _mm512_set1_pd
#define gmx_simd_store_d           _mm512_store_pd
#define gmx_simd_loadu_d           _mm512_loadu_pd
#define gmx_simd_storeu_d          _mm512_storeu_pd
#define gmx_simd_setzero_d         _mm512_setzero_pd
#define gmx_simd_add_d             _mm512_add_pd
#define gmx_simd_sub_d             _mm512_sub_pd
#define gmx_simd_mul_d             _mm512_mul_pd
#define gmx_simd_fmadd_d           _mm512_fmadd_pd
#define gmx_simd_fmsub_d           _mm512_fmsub_pd
#define gmx_simd_fnmadd_d          _mm512_fnmadd_pd
#define gmx_simd_fnmsub_d          _mm512_fnmsub_pd
#define gmx_simd_and_d(a, b)       _mm512_castsi512_pd(_mm512_and_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_andnot_d(a, b)    _mm512_castsi512_pd(_mm512_andnot_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_or_d(a, b)        _mm512_castsi512_pd(_mm512_or_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_xor_d(a, b)       _mm512_castsi512_pd(_mm512_xor_epi32(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd_rsqrt_d           _mm512_rsqrt14_pd
#define gmx_simd_rcp_d             _mm512_rcp14_pd
#define gmx_simd_fabs_d(x)         _mm512_abs_pd(x)
#define gmx_simd_fneg_d(x)         gmx_simd_xor_d(x, _mm512_set1_pd(-0.0))
#define gmx_simd_max_d             _mm512_max_pd
#define gmx_simd_min_d             _mm512_min_pd
#define gmx_simd_round_d(x)        _mm512_roundscale_pd(x, 0)
#define gmx_simd_trunc_d(x)        _mm512_trunc_pd(x)
#define gmx_simd_fraction_d(x)     _mm512_sub_pd(x, gmx_simd_trunc_d(x))
#define gmx_simd_get_exponent_d(x) _mm512_getexp_pd(x)
#define gmx_simd_get_mantissa_d(x) _mm512_getmant_pd(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero)
#define gmx_simd_set_exponent_d(x) gmx_simd_set_exponent_d_x86_avx_512f(x)
/* integer datatype corresponding to double: gmx_simd_dint32_t
   Doesn't use mask other than where required. No side effect expected for operating on the (unused) upper 8.
 */
#define gmx_simd_dint32_t          __m256i
#define gmx_simd_load_di(m)        _mm256_load_si256((const __m256i *)m)
#define gmx_simd_set1_di           _mm256_set1_epi32
#define gmx_simd_store_di(m, x)     _mm256_store_si256((__m256i *)m, x)
#define gmx_simd_loadu_di(m)       _mm256_loadu_si256((const __m256i *)m)
#define gmx_simd_storeu_di(m, x)    _mm256_storeu_si256((__m256i *)m, x)
#undef  gmx_simd_extract_di
#define gmx_simd_setzero_di        _mm256_setzero_si256
#define gmx_simd_cvt_d2i           _mm512_cvtpd_epi32
#define gmx_simd_cvtt_d2i          _mm512_cvttpd_epi32
#define gmx_simd_cvt_i2d           _mm512_cvtepi32_pd
/* Integer logical ops on gmx_simd_dint32_t */
#define gmx_simd_slli_di           _mm256_slli_epi32
#define gmx_simd_srli_di           _mm256_srli_epi32
#define gmx_simd_and_di            _mm256_and_si256
#define gmx_simd_andnot_di         _mm256_andnot_si256
#define gmx_simd_or_di             _mm256_or_si256
#define gmx_simd_xor_di            _mm256_xor_si256
/* Integer arithmetic ops on gmx_simd_dint32_t */
#define gmx_simd_add_di            _mm256_add_epi32
#define gmx_simd_sub_di            _mm256_sub_epi32
#define gmx_simd_mul_di            _mm256_mullo_epi32
/* Boolean & comparison operations on gmx_simd_double_t */
#define gmx_simd_dbool_t           __mmask8
#define gmx_simd_cmpeq_d(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_EQ_OQ)
#define gmx_simd_cmplt_d(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_LT_OS)
#define gmx_simd_cmple_d(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_LE_OS)
#define gmx_simd_and_db            _mm512_kand
#define gmx_simd_or_db             _mm512_kor
#define gmx_simd_anytrue_db(x)     _mm512_mask2int(x)
#define gmx_simd_blendzero_d(a, sel)    _mm512_mask_mov_pd(_mm512_setzero_pd(), sel, a)
#define gmx_simd_blendnotzero_d(a, sel) _mm512_mask_mov_pd(a, sel, _mm512_setzero_pd())
#define gmx_simd_blendv_d(a, b, sel)    _mm512_mask_blend_pd(sel, a, b)
#define gmx_simd_reduce_d(a)       gmx_simd_reduce_d_x86_avx_512f(a)
/* Boolean & comparison operations on gmx_simd_dint32_t */
#define gmx_simd_dibool_t          __mmask16
#define gmx_simd_cmpeq_di(a, b)    _mm512_mask_cmp_epi32_mask(_mm512_int2mask(0xFF), _mm512_castsi256_si512(a), _mm512_castsi256_si512(b), _MM_CMPINT_EQ)
#define gmx_simd_cmplt_di(a, b)    _mm512_mask_cmp_epi32_mask(_mm512_int2mask(0xFF), _mm512_castsi256_si512(a), _mm512_castsi256_si512(b), _MM_CMPINT_LT)
#define gmx_simd_and_dib           _mm512_kand
#define gmx_simd_or_dib            _mm512_kor
#define gmx_simd_anytrue_dib(x)    (_mm512_mask2int(x)&0xFF)
#define gmx_simd_blendzero_di(a, sel)    _mm512_castsi512_si256(_mm512_mask_mov_epi32(_mm512_setzero_si512(), sel, _mm512_castsi256_si512(a)))
#define gmx_simd_blendnotzero_di(a, sel) _mm512_castsi512_si256(_mm512_mask_mov_epi32(_mm512_castsi256_si512(a), sel, _mm512_setzero_si512()))
#define gmx_simd_blendv_di(a, b, sel)    _mm512_castsi512_si256(_mm512_mask_blend_epi32(sel, _mm512_castsi256_si512(a), _mm512_castsi256_si512(b)))
/* Conversions between booleans. Double & dint stuff is stored in low bits */
#define gmx_simd_cvt_db2dib(x)     (x)
#define gmx_simd_cvt_dib2db(x)     (x)

/* Float/double conversion */
#define gmx_simd_cvt_f2dd          gmx_simd_cvt_f2dd_x86_avx_512f
#define gmx_simd_cvt_dd2f          gmx_simd_cvt_dd2f_x86_avx_512f




/****************************************************
 *      SINGLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
/* We use __m128 to only access part of the registers, but for a few operations
 * we cast to the full register width when those operations are cheaper. We
 * also save some register space by using mask registers for the booleans.
 */
#define gmx_simd4_float_t           __m128
#define gmx_simd4_load_f            _mm_load_ps
#define gmx_simd4_load1_f           _mm_load1_ps
#define gmx_simd4_set1_f            _mm_set1_ps
#define gmx_simd4_store_f           _mm_store_ps
#define gmx_simd4_loadu_f           _mm_loadu_ps
#define gmx_simd4_storeu_f          _mm_storeu_ps
#define gmx_simd4_setzero_f         _mm_setzero_ps
#define gmx_simd4_add_f             _mm_add_ps
#define gmx_simd4_sub_f             _mm_sub_ps
#define gmx_simd4_mul_f             _mm_mul_ps
#define gmx_simd4_fmadd_f           _mm_fmadd_ps
#define gmx_simd4_fmsub_f           _mm_fmsub_ps
#define gmx_simd4_fnmadd_f          _mm_fnmadd_ps
#define gmx_simd4_fnmsub_f          _mm_fnmsub_ps
#define gmx_simd4_and_f             _mm_and_ps
#define gmx_simd4_andnot_f          _mm_andnot_ps
#define gmx_simd4_or_f              _mm_or_ps
#define gmx_simd4_xor_f             _mm_xor_ps
/* We need to use the new table lookup instructions since we have specified
 * 14 bits of accuracy for AVX-512F.
 */
#define gmx_simd4_rsqrt_f(x)        _mm512_castps512_ps128(_mm512_rsqrt14_ps(_mm512_castps128_ps512(x)))
/* abs/neg cannot cause FP exceptions, so we can operate on entire register */
#define gmx_simd4_fabs_f(x)         _mm512_castps512_ps128(_mm512_abs_ps(_mm512_castps128_ps512(x)))
#define gmx_simd4_fneg_f(x)         _mm_xor_ps(x, _mm_set1_ps(-0.0f))
#define gmx_simd4_max_f             _mm_max_ps
#define gmx_simd4_min_f             _mm_min_ps
#define gmx_simd4_round_f(x)        _mm_round_ps(x, _MM_FROUND_NINT)
#define gmx_simd4_trunc_f(x)        _mm_round_ps(x, _MM_FROUND_TRUNC)
#define gmx_simd4_dotproduct3_f(a, b) gmx_simd4_dotproduct3_f_x86_avx_512f(a, b)
#define gmx_simd4_fbool_t           __mmask16
#define gmx_simd4_cmpeq_f(a, b)     _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), _mm512_castps128_ps512(a), _mm512_castps128_ps512(b), _CMP_EQ_OQ)
#define gmx_simd4_cmplt_f(a, b)     _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), _mm512_castps128_ps512(a), _mm512_castps128_ps512(b), _CMP_LT_OS)
#define gmx_simd4_cmple_f(a, b)     _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), _mm512_castps128_ps512(a), _mm512_castps128_ps512(b), _CMP_LE_OS)
#define gmx_simd4_and_fb            _mm512_kand
#define gmx_simd4_or_fb             _mm512_kor
#define gmx_simd4_anytrue_fb(x)     (_mm512_mask2int(x)&0xF)
#define gmx_simd4_blendzero_f(a, sel)    _mm512_castps512_ps128(_mm512_mask_mov_ps(_mm512_setzero_ps(), sel, _mm512_castps128_ps512(a)))
#define gmx_simd4_blendnotzero_f(a, sel) _mm512_castps512_ps128(_mm512_mask_mov_ps(_mm512_setzero_ps(), _mm512_knot(sel), _mm512_castps128_ps512(a)))
#define gmx_simd4_blendv_f(a, b, sel)    _mm512_castps512_ps128(_mm512_mask_blend_ps(sel, _mm512_castps128_ps512(a), _mm512_castps128_ps512(b)))
#define gmx_simd4_reduce_f(x)       gmx_simd4_reduce_f_x86_avx_512f(x)

/****************************************************
 *      DOUBLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
/* We use __m256d to only access part of the registers, but for a few operations
 * we cast to the full register width when those operations are cheaper. We
 * also save some register space by using mask registers for the booleans.
 */
#define gmx_simd4_double_t          __m256d
#define gmx_simd4_load_d            _mm256_load_pd
#define gmx_simd4_load1_d           _mm256_broadcast_sd
#define gmx_simd4_set1_d            _mm256_set1_pd
#define gmx_simd4_store_d           _mm256_store_pd
#define gmx_simd4_loadu_d           _mm256_loadu_pd
#define gmx_simd4_storeu_d          _mm256_storeu_pd
#define gmx_simd4_setzero_d         _mm256_setzero_pd
#define gmx_simd4_add_d             _mm256_add_pd
#define gmx_simd4_sub_d             _mm256_sub_pd
#define gmx_simd4_mul_d             _mm256_mul_pd
#define gmx_simd4_fmadd_d           _mm256_fmadd_pd
#define gmx_simd4_fmsub_d           _mm256_fmsub_pd
#define gmx_simd4_fnmadd_d          _mm256_fnmadd_pd
#define gmx_simd4_fnmsub_d          _mm256_fnmsub_pd
#define gmx_simd4_and_d             _mm256_and_pd
#define gmx_simd4_andnot_d          _mm256_andnot_pd
#define gmx_simd4_or_d              _mm256_or_pd
#define gmx_simd4_xor_d             _mm256_xor_pd
/* We need to use the new table lookup instructions since we have specified
 * 14 bits of accuracy for AVX-512F.
 */
#define gmx_simd4_rsqrt_d(x)        _mm512_castpd512_pd256(_mm512_rsqrt14_pd(_mm512_castpd256_pd512(x)))
/* abs/neg cannot cause FP exceptions, so we can operate on entire register */
#define gmx_simd4_fabs_d(x)         _mm512_castpd512_pd256(_mm512_abs_pd(_mm512_castpd256_pd512(x)))
#define gmx_simd4_fneg_d(x)         _mm256_xor_pd(x, _mm256_set1_pd(-0.0))
#define gmx_simd4_max_d(a, b)       _mm256_max_pd(a, b)
#define gmx_simd4_min_d(a, b)       _mm256_min_pd(a, b)
#define gmx_simd4_round_d(a)        _mm256_round_pd(a, _MM_FROUND_NINT)
#define gmx_simd4_trunc_d(a)        _mm256_round_pd(a, _MM_FROUND_TRUNC)
#define gmx_simd4_dotproduct3_d(a, b) gmx_simd4_dotproduct3_d_x86_avx_512f(a, b)
#define gmx_simd4_dbool_t           __mmask16
#define gmx_simd4_cmpeq_d(a, b)     _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), _mm512_castpd256_pd512(a), _mm512_castpd256_pd512(b), _CMP_EQ_OQ)
#define gmx_simd4_cmplt_d(a, b)     _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), _mm512_castpd256_pd512(a), _mm512_castpd256_pd512(b), _CMP_LT_OS)
#define gmx_simd4_cmple_d(a, b)     _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), _mm512_castpd256_pd512(a), _mm512_castpd256_pd512(b), _CMP_LE_OS)
#define gmx_simd4_and_db            _mm512_kand
#define gmx_simd4_or_db             _mm512_kor
#define gmx_simd4_anytrue_db(x)     (_mm512_mask2int(x)&0xF)
#define gmx_simd4_blendzero_d(a, sel)    _mm512_castpd512_pd256(_mm512_mask_mov_pd(_mm512_setzero_pd(), sel, _mm512_castpd256_pd512(a)))
#define gmx_simd4_blendnotzero_d(a, sel) _mm512_castpd512_pd256(_mm512_mask_mov_pd(_mm512_setzero_pd(), _mm512_knot(sel), _mm512_castpd256_pd512(a)))
#define gmx_simd4_blendv_d(a, b, sel)    _mm512_castpd512_pd256(_mm512_mask_blend_pd(sel, _mm512_castpd256_pd512(a), _mm512_castpd256_pd512(b)))
#define gmx_simd4_reduce_d(x)       gmx_simd4_reduce_d_x86_avx_512f(x)


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

static gmx_inline __m512d
gmx_simd_set_exponent_d_x86_avx_512f(__m512d a)
{
    const __m512i expbias      = _mm512_set1_epi32(1023);
    __m512i       iexp         = _mm512_castsi256_si512(gmx_simd_cvt_d2i(a));

    iexp = _mm512_permutevar_epi32(_mm512_set_epi32(7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 0, 0), iexp);
    iexp = _mm512_mask_slli_epi32(_mm512_setzero_epi32(), _mm512_int2mask(0xAAAA), _mm512_add_epi32(iexp, expbias), 20);
    return _mm512_castsi512_pd(iexp);
}

static gmx_inline double
gmx_simd_reduce_d_x86_avx_512f(__m512d a)
{
    __m128d b;
    a = _mm512_add_pd(a, _mm512_shuffle_f64x2(a, a, _MM_PERM_DCDC));
    a = _mm512_add_pd(a, _mm512_shuffle_f64x2(a, a, _MM_PERM_ABAB));
    b = _mm512_castpd512_pd128(a);
    b = _mm_hadd_pd(b, b);
    return _mm_cvtsd_f64(b);
}

static gmx_inline void
gmx_simd_cvt_f2dd_x86_avx_512f(__m512 f, __m512d * d0, __m512d * d1)
{
    *d0 = _mm512_cvtpslo_pd(f);
    *d1 = _mm512_cvtpslo_pd(_mm512_shuffle_f32x4(f, f, _MM_PERM_DCDC));
}

static gmx_inline __m512
gmx_simd_cvt_dd2f_x86_avx_512f(__m512d d0, __m512d d1)
{
    __m512 f0 = _mm512_cvtpd_pslo(d0);
    __m512 f1 = _mm512_cvtpd_pslo(d1);
    return _mm512_shuffle_f32x4(f0, f1, _MM_PERM_BABA);
}

static gmx_inline float
gmx_simd4_reduce_f_x86_avx_512f(__m128 a)
{
    float f;
    a = _mm_hadd_ps(a, a);
    a = _mm_hadd_ps(a, a);
    _mm_store_ss(&f, a);
    return f;
}

static gmx_inline float
gmx_simd4_dotproduct3_f_x86_avx_512f(__m128 a, __m128 b)
{
    float  f;
    __m128 c;
    a = _mm_mul_ps(a, b);
    c = _mm_add_ps(a, _mm_permute_ps(a, _MM_SHUFFLE(0, 3, 2, 1)));
    c = _mm_add_ps(c, _mm_permute_ps(a, _MM_SHUFFLE(1, 0, 3, 2)));
    _mm_store_ss(&f, c);
    return f;
}

static gmx_inline double
gmx_simd4_reduce_d_x86_avx_512f(__m256d a)
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

static gmx_inline double
gmx_simd4_dotproduct3_d_x86_avx_512f(__m256d a, __m256d b)
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

#endif /* GMX_SIMD_IMPL_X86_AVX_512F_H */
