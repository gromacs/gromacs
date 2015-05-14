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

#include "config.h"

#include <assert.h>
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

#define GMX_SIMD_V2

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
#define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT
#define GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT

#define GMX_SIMD_HAVE_DINT32
#undef  GMX_SIMD_HAVE_DINT32_EXTRACT
#define GMX_SIMD_HAVE_DINT32_LOGICAL
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS
#define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_DOUBLE
#define GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE
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
#define gmx_simd_load1_f(m)        _mm512_set1_ps(*m)
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
#define gmx_simd_mul_mask_f(a, b, m)       _mm512_maskz_mul_ps(m, a, b)
#define gmx_simd_fmadd_mask_f(a, b, c, m)  _mm512_maskz_fmadd_ps(m, a, b, c)
#define gmx_simd_rcp_mask_f(a, m)          _mm512_maskz_rcp14_ps(m, a)
#define gmx_simd_rsqrt_mask_f(a, m)        _mm512_maskz_rsqrt14_ps(m, a)
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
#define gmx_simd_cmpnz_f(a)        _mm512_test_epi32_mask(_mm512_castps_si512(a), _mm512_castps_si512(a))
#define gmx_simd_cmplt_f(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_LT_OS)
#define gmx_simd_cmple_f(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_LE_OS)
#define gmx_simd_and_fb            _mm512_kand
#define gmx_simd_andnot_fb(a, b)   _mm512_kandn(a, b)
#define gmx_simd_or_fb             _mm512_kor
#define gmx_simd_anytrue_fb        _mm512_mask2int
#define gmx_simd_blendzero_f(a, sel)    _mm512_maskz_mov_ps(sel, a)
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
#define gmx_simd_blendzero_fi(a, sel)    _mm512_maskz_mov_epi32(sel, a)
#define gmx_simd_blendnotzero_fi(a, sel) _mm512_mask_mov_epi32(a, sel, _mm512_setzero_epi32())
#define gmx_simd_blendv_fi(a, b, sel)    _mm512_mask_blend_epi32(sel, a, b)
/* Conversions between different booleans */
#define gmx_simd_cvt_fb2fib(x)     (x)
#define gmx_simd_cvt_fib2fb(x)     (x)
/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_f             gmx_simd_gather_load_transpose_f_x86_avx_512f
#define gmx_simd_best_pair_alignment_f               gmx_simd_best_pair_alignment_f_x86_avx_512f
#define gmx_simd_gather_loadu_transpose_f            gmx_simd_gather_loadu_transpose_f_x86_avx_512f
#define gmx_simd_transpose_scatter_storeu_f          gmx_simd_transpose_scatter_storeu_f_x86_avx_512f
#define gmx_simd_transpose_scatter_incru_f           gmx_simd_transpose_scatter_incru_f_x86_avx_512f
#define gmx_simd_transpose_scatter_decru_f           gmx_simd_transpose_scatter_decru_f_x86_avx_512f
#define gmx_simd_expand_scalars_to_triplets_f        gmx_simd_expand_scalars_to_triplets_f_x86_avx_512f
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_x86_avx_512f
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_x86_avx_512f
/* Half-simd-width utilities */
#define gmx_simd_load_dual_hsimd_f                   gmx_simd_load_dual_hsimd_f_x86_avx_512f
#define gmx_simd_loaddup_hsimd_f                     gmx_simd_loaddup_hsimd_f_x86_avx_512f
#define gmx_simd_load1_dual_hsimd_f                  gmx_simd_load1_dual_hsimd_f_x86_avx_512f
#define gmx_simd_store_dual_hsimd_f                  gmx_simd_store_dual_hsimd_f_x86_avx_512f
#define gmx_simd_decr_hsimd_f                        gmx_simd_decr_hsimd_f_x86_avx_512f
#define gmx_simd_gather_load_transpose_hsimd_f       gmx_simd_gather_load_transpose_hsimd_f_x86_avx_512f
#define gmx_simd_reduce_incr_4_return_sum_hsimd_f    gmx_simd_reduce_incr_4_return_sum_hsimd_f_x86_avx_512f





/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_double_t          __m512d
#define gmx_simd_load_d            _mm512_load_pd
/* Avoid using _mm512_extload_pd() since it is not available on gcc-4.9 */
#define gmx_simd_load1_d(m)        _mm512_set1_pd(*m)
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
#define gmx_simd_mul_mask_d(a, b, m)       _mm512_maskz_mul_pd(m, a, b)
#define gmx_simd_fmadd_mask_d(a, b, c, m)  _mm512_maskz_fmadd_pd(m, a, b, c)
#define gmx_simd_rcp_mask_d(a, m)          _mm512_maskz_rcp14_pd(m, a)
#define gmx_simd_rsqrt_mask_d(a, m)        _mm512_maskz_rsqrt14_pd(m, a)
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
#define gmx_simd_dint32_t          __m512i

#define gmx_simd_load_di(m)        _mm512_maskz_loadu_epi32(_mm512_int2mask(0x00FF), m)
#define gmx_simd_set1_di           _mm512_set1_epi32
#define gmx_simd_store_di(m, x)    _mm512_mask_storeu_epi32(m, _mm512_int2mask(0x00FF), x)
#define gmx_simd_loadu_di(m)       _mm512_maskz_loadu_epi32(_mm512_int2mask(0x00FF), m)
#define gmx_simd_storeu_di(m, x)   _mm512_mask_storeu_epi32(m, _mm512_int2mask(0x00FF), x)
#undef  gmx_simd_extract_di
#define gmx_simd_setzero_di        _mm512_setzero_si512
#define gmx_simd_cvt_d2i(x)        _mm512_castsi256_si512(_mm512_cvtpd_epi32(x))
#define gmx_simd_cvtt_d2i(x)       _mm512_castsi256_si512(_mm512_cvttpd_epi32(x))
#define gmx_simd_cvt_i2d(x)        _mm512_cvtepi32_pd(_mm512_castsi512_si256(x))
/* Integer logical ops on gmx_simd_dint32_t */
#define gmx_simd_slli_di           _mm512_slli_epi32
#define gmx_simd_srli_di           _mm512_srli_epi32
#define gmx_simd_and_di            _mm512_and_si512
#define gmx_simd_andnot_di         _mm512_andnot_si512
#define gmx_simd_or_di             _mm512_or_si512
#define gmx_simd_xor_di            _mm512_xor_si512
/* Integer arithmetic ops on gmx_simd_dint32_t */
#define gmx_simd_add_di            _mm512_add_epi32
#define gmx_simd_sub_di            _mm512_sub_epi32
#define gmx_simd_mul_di            _mm512_mullo_epi32
/* Boolean & comparison operations on gmx_simd_double_t */
#define gmx_simd_dbool_t           __mmask8
#define gmx_simd_cmpeq_d(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_EQ_OQ)
#define gmx_simd_cmpnz_d(a)        _mm512_test_epi64_mask(_mm512_castpd_si512(a), _mm512_castpd_si512(a))
#define gmx_simd_cmplt_d(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_LT_OS)
#define gmx_simd_cmple_d(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_LE_OS)
#define gmx_simd_and_db            _mm512_kand
#define gmx_simd_or_db             _mm512_kor
#define gmx_simd_anytrue_db(x)     _mm512_mask2int(x)
#define gmx_simd_blendzero_d(a, sel)    _mm512_maskz_mov_pd(sel, a)
#define gmx_simd_blendnotzero_d(a, sel) _mm512_mask_mov_pd(a, sel, _mm512_setzero_pd())
#define gmx_simd_blendv_d(a, b, sel)    _mm512_mask_blend_pd(sel, a, b)
#define gmx_simd_reduce_d(a)       gmx_simd_reduce_d_x86_avx_512f(a)
/* Boolean & comparison operations on gmx_simd_dint32_t */
#define gmx_simd_dibool_t          __mmask16
#define gmx_simd_cmpeq_di(a, b)    _mm512_mask_cmp_epi32_mask(_mm512_int2mask(0x00FF), a, b, _MM_CMPINT_EQ)
#define gmx_simd_cmplt_di(a, b)    _mm512_mask_cmp_epi32_mask(_mm512_int2mask(0x00FF), a, b, _MM_CMPINT_LT)
#define gmx_simd_and_dib           _mm512_kand
#define gmx_simd_or_dib            _mm512_kor
#define gmx_simd_anytrue_dib(x)    (_mm512_mask2int(x)&0x00FF)
#define gmx_simd_blendzero_di(a, sel)    _mm512_maskz_mov_epi32(sel, a)
#define gmx_simd_blendnotzero_di(a, sel) _mm512_mask_mov_epi32(a, sel, _mm512_setzero_si512())
#define gmx_simd_blendv_di(a, b, sel)    _mm512_mask_blend_epi32(sel, a, b)
/* Conversions between booleans. Double & dint stuff is stored in low bits */
#define gmx_simd_cvt_db2dib(x)     (x)
#define gmx_simd_cvt_dib2db(x)     (x)

/* Float/double conversion */
#define gmx_simd_cvt_f2dd          gmx_simd_cvt_f2dd_x86_avx_512f
#define gmx_simd_cvt_dd2f          gmx_simd_cvt_dd2f_x86_avx_512f
/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_d             gmx_simd_gather_load_transpose_d_x86_avx_512f
#define gmx_simd_best_pair_alignment_d               gmx_simd_best_pair_alignment_d_x86_avx_512f
#define gmx_simd_gather_loadu_transpose_d            gmx_simd_gather_loadu_transpose_d_x86_avx_512f
#define gmx_simd_transpose_scatter_storeu_d          gmx_simd_transpose_scatter_storeu_d_x86_avx_512f
#define gmx_simd_transpose_scatter_incru_d           gmx_simd_transpose_scatter_incru_d_x86_avx_512f
#define gmx_simd_transpose_scatter_decru_d           gmx_simd_transpose_scatter_decru_d_x86_avx_512f
#define gmx_simd_expand_scalars_to_triplets_d        gmx_simd_expand_scalars_to_triplets_d_x86_avx_512f
#define gmx_simd_gather_load_bysimdint_transpose_d   gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f
#define gmx_simd_gather_loadu_bysimdint_transpose_d  gmx_simd_gather_loadu_bysimdint_transpose_d_x86_avx_512f
#define gmx_simd_reduce_incr_4_return_sum_d          gmx_simd_reduce_incr_4_return_sum_d_x86_avx_512f
/* Half-simd-width utilities */
#define gmx_simd_load_dual_hsimd_d                   gmx_simd_load_dual_hsimd_d_x86_avx_512f
#define gmx_simd_loaddup_hsimd_d                     gmx_simd_loaddup_hsimd_d_x86_avx_512f
#define gmx_simd_load1_dual_hsimd_d                  gmx_simd_load1_dual_hsimd_d_x86_avx_512f
#define gmx_simd_store_dual_hsimd_d                  gmx_simd_store_dual_hsimd_d_x86_avx_512f
#define gmx_simd_decr_hsimd_d                        gmx_simd_decr_hsimd_d_x86_avx_512f
#define gmx_simd_gather_load_transpose_hsimd_d       gmx_simd_gather_load_transpose_hsimd_d_x86_avx_512f
#define gmx_simd_reduce_incr_4_return_sum_hsimd_d    gmx_simd_reduce_incr_4_return_sum_hsimd_d_x86_avx_512f

/****************************************************
 *      SINGLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
#define gmx_simd4_float_t           __m512
#define gmx_simd4_mask              _mm512_int2mask(0xF)
#define gmx_simd4_load_f(m)         _mm512_maskz_loadu_ps(gmx_simd4_mask, m)
#define gmx_simd4_load1_f(m)        _mm512_permute_ps(_mm512_maskz_load_ps(_mm512_int2mask(0x0001), m), 0b00000000)
#define gmx_simd4_set1_f(x)         _mm512_set1_ps(x)
#define gmx_simd4_store_f(m, x)     _mm512_mask_storeu_ps(m, gmx_simd4_mask, x)
#define gmx_simd4_loadu_f(m)        _mm512_maskz_loadu_ps(gmx_simd4_mask, m)
#define gmx_simd4_storeu_f(m, x)    _mm512_mask_storeu_ps(m, gmx_simd4_mask, x)
#define gmx_simd4_setzero_f         _mm512_setzero_ps
#define gmx_simd4_add_f(a, b)        _mm512_maskz_add_ps(gmx_simd4_mask, a, b)
#define gmx_simd4_sub_f(a, b)        _mm512_maskz_sub_ps(gmx_simd4_mask, a, b)
#define gmx_simd4_mul_f(a, b)        _mm512_maskz_mul_ps(gmx_simd4_mask, a, b)
#define gmx_simd4_fmadd_f(a, b, c)    _mm512_maskz_fmadd_ps(gmx_simd4_mask, a, b, c)
#define gmx_simd4_fmsub_f(a, b, c)    _mm512_maskz_fmsub_ps(gmx_simd4_mask, a, b, c)
#define gmx_simd4_fnmadd_f(a, b, c)   _mm512_maskz_fnmadd_ps(gmx_simd4_mask, a, b, c)
#define gmx_simd4_fnmsub_f(a, b, c)   _mm512_maskz_fnmsub_ps(gmx_simd4_mask, a, b, c)
#define gmx_simd4_and_f(a, b)       _mm512_castsi512_ps(_mm512_maskz_and_epi32(gmx_simd4_mask, _mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd4_andnot_f(a, b)    _mm512_castsi512_ps(_mm512_maskz_andnot_epi32(gmx_simd4_mask, _mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd4_or_f(a, b)        _mm512_castsi512_ps(_mm512_maskz_or_epi32(gmx_simd4_mask, _mm512_castps_si512(a), _mm512_castps_si512(b)))
#define gmx_simd4_xor_f(a, b)       _mm512_castsi512_ps(_mm512_maskz_xor_epi32(gmx_simd4_mask, _mm512_castps_si512(a), _mm512_castps_si512(b)))
/* We need to use the new table lookup instructions since we have specified
 * 14 bits of accuracy for AVX-512F.
 */
#define gmx_simd4_rsqrt_f(x)        _mm512_maskz_rsqrt14_ps(gmx_simd4_mask, x)
#define gmx_simd4_fabs_f(x)         _mm512_abs_ps(x)
#define gmx_simd4_fneg_f(x)         gmx_simd4_xor_f(x, gmx_simd4_set1_f(-0.0f))
#define gmx_simd4_max_f(a, b)        _mm512_maskz_max_ps(gmx_simd4_mask, a, b)
#define gmx_simd4_min_f(a, b)        _mm512_maskz_min_ps(gmx_simd4_mask, a, b)
#define gmx_simd4_round_f(x)        _mm512_maskz_roundscale_ps(gmx_simd4_mask, x, 0)
#define gmx_simd4_trunc_f(x)        _mm512_mask_trunc_ps(gmx_simd4_setzero_f(), gmx_simd4_mask, x)
#define gmx_simd4_dotproduct3_f(a, b) gmx_simd4_dotproduct3_f_x86_avx_512f(a, b)
#define gmx_simd4_transpose_f       gmx_simd4_transpose_f_x86_avx_512f
#define gmx_simd4_fbool_t           __mmask16
#define gmx_simd4_cmpeq_f(a, b)     _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), a, b, _CMP_EQ_OQ)
#define gmx_simd4_cmplt_f(a, b)     _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), a, b, _CMP_LT_OS)
#define gmx_simd4_cmple_f(a, b)     _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), a, b, _CMP_LE_OS)
#define gmx_simd4_and_fb            _mm512_kand
#define gmx_simd4_or_fb             _mm512_kor
#define gmx_simd4_anytrue_fb(x)     (_mm512_mask2int(x)&0xF)
#define gmx_simd4_blendzero_f(a, sel)    _mm512_maskz_mov_ps(sel, a)
#define gmx_simd4_blendnotzero_f(a, sel) _mm512_maskz_mov_ps(_mm512_knot(sel), a)
#define gmx_simd4_blendv_f(a, b, sel)    _mm512_mask_blend_ps(sel, a, b)
#define gmx_simd4_reduce_f(x)       gmx_simd4_reduce_f_x86_avx_512f(x)

/****************************************************
 *      DOUBLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
#define gmx_simd4_double_t          __m512d
#define gmx_simd4_load_d(m)         _mm512_maskz_loadu_pd(gmx_simd4_mask, m)
#define gmx_simd4_load1_d(m)        _mm512_permutex_pd(_mm512_maskz_load_pd(_mm512_int2mask(0x01), m), 0b00000000)
#define gmx_simd4_set1_d            _mm512_set1_pd
#define gmx_simd4_store_d(m, x)      _mm512_mask_storeu_pd(m, gmx_simd4_mask, x)
#define gmx_simd4_loadu_d(m)        _mm512_maskz_loadu_pd(gmx_simd4_mask, m)
#define gmx_simd4_storeu_d(m, x)     _mm512_mask_storeu_pd(m, gmx_simd4_mask, x)
#define gmx_simd4_setzero_d         _mm512_setzero_pd
#define gmx_simd4_add_d(a, b)        _mm512_maskz_add_pd(gmx_simd4_mask, a, b)
#define gmx_simd4_sub_d(a, b)        _mm512_maskz_sub_pd(gmx_simd4_mask, a, b)
#define gmx_simd4_mul_d(a, b)        _mm512_maskz_mul_pd(gmx_simd4_mask, a, b)
#define gmx_simd4_fmadd_d(a, b, c)    _mm512_maskz_fmadd_pd(gmx_simd4_mask, a, b, c)
#define gmx_simd4_fmsub_d(a, b, c)    _mm512_maskz_fmsub_pd(gmx_simd4_mask, a, b, c)
#define gmx_simd4_fnmadd_d(a, b, c)   _mm512_maskz_fnmadd_pd(gmx_simd4_mask, a, b, c)
#define gmx_simd4_fnmsub_d(a, b, c)   _mm512_maskz_fnmsub_pd(gmx_simd4_mask, a, b, c)
#define gmx_simd4_and_d(a, b)       _mm512_castsi512_pd(_mm512_maskz_and_epi64(gmx_simd4_mask, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd4_andnot_d(a, b)    _mm512_castsi512_pd(_mm512_maskz_andnot_epi64(gmx_simd4_mask, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd4_or_d(a, b)        _mm512_castsi512_pd(_mm512_maskz_or_epi64(gmx_simd4_mask, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd4_xor_d(a, b)       _mm512_castsi512_pd(_mm512_maskz_xor_epi64(gmx_simd4_mask, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
/* We need to use the new table lookup instructions since we have specified
 * 14 bits of accuracy for AVX-512F.
 */
#define gmx_simd4_rsqrt_d(x)        _mm512_rsqrt14_pd(x)
/* abs/neg cannot cause FP exceptions, so we can operate on entire register */
#define gmx_simd4_fabs_d(x)         _mm512_abs_pd(x)
#define gmx_simd4_fneg_d(x)         gmx_simd4_xor_d(x, _mm512_set1_pd(-0.0))
#define gmx_simd4_max_d(a, b)       _mm512_maskz_max_pd(gmx_simd4_mask, a, b)
#define gmx_simd4_min_d(a, b)       _mm512_maskz_min_pd(gmx_simd4_mask, a, b)
#define gmx_simd4_round_d(a)        _mm512_maskz_roundscale_pd(gmx_simd4_mask, a, 0)
#define gmx_simd4_trunc_d(a)        _mm512_mask_trunc_pd(gmx_simd4_setzero_d(), gmx_simd4_mask, a)
#define gmx_simd4_dotproduct3_d(a, b) gmx_simd4_dotproduct3_d_x86_avx_512f(a, b)
#define gmx_simd4_transpose_d       gmx_simd4_transpose_d_x86_avx_512f
#define gmx_simd4_dbool_t           __mmask16
#define gmx_simd4_cmpeq_d(a, b)     _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), a, b, _CMP_EQ_OQ)
#define gmx_simd4_cmplt_d(a, b)     _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), a, b, _CMP_LT_OS)
#define gmx_simd4_cmple_d(a, b)     _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), a, b, _CMP_LE_OS)
#define gmx_simd4_and_db            _mm512_kand
#define gmx_simd4_or_db             _mm512_kor
#define gmx_simd4_anytrue_db(x)     (_mm512_mask2int(x)&0xF)
#define gmx_simd4_blendzero_d(a, sel)    _mm512_maskz_mov_pd(sel, a)
#define gmx_simd4_blendnotzero_d(a, sel) _mm512_maskz_mov_pd(_mm512_knot(sel), a)
#define gmx_simd4_blendv_d(a, b, sel)    _mm512_mask_blend_pd(sel, a, b)
#define gmx_simd4_reduce_d(x)       gmx_simd4_reduce_d_x86_avx_512f(x)


/*********************************************************
 * SIMD SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 *********************************************************/
static gmx_inline __m512
gmx_simd_set_exponent_f_x86_avx_512f(__m512 a)
{
    const __m512i expbias      = _mm512_set1_epi32(127);
    __m512i       iexp         = gmx_simd_cvt_f2i(a);

    /* This is likely faster than the built in scale operation (lat 8, t-put 3)
     * since we only work on the integer part and use shifts, meaning it will run
     * on the integer ports that are typically less utilized in our kernels.
     */
    iexp = _mm512_slli_epi32(_mm512_add_epi32(iexp, expbias), 23);
    return _mm512_castsi512_ps(iexp);
}

static gmx_inline float
gmx_simd_reduce_f_x86_avx_512f(__m512 a)
{
    float f;
    a = _mm512_add_ps(a, _mm512_shuffle_f32x4(a, a, 0b11101110));
    a = _mm512_add_ps(a, _mm512_shuffle_f32x4(a, a, 0b00010001));
    a = _mm512_add_ps(a, _mm512_permute_ps(a, 0b11101110));
    a = _mm512_add_ps(a, _mm512_permute_ps(a, 0b00010001));
    _mm512_mask_storeu_ps(&f, _mm512_int2mask(0x0001), a);
    return f;
}


/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

/* Declare functions implemented further down that are used here */
template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f(const float *       base,
                                                        gmx_simd_fint32_t   simdoffset,
                                                        gmx_simd_float_t   &v0,
                                                        gmx_simd_float_t   &v1,
                                                        gmx_simd_float_t   &v2,
                                                        gmx_simd_float_t   &v3);

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f(const float *       base,
                                                        gmx_simd_fint32_t   simdoffset,
                                                        gmx_simd_float_t   &v0,
                                                        gmx_simd_float_t   &v1);


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_x86_avx_512f(const float *        base,
                                              const gmx_int32_t    offset[],
                                              gmx_simd_float_t    &v0,
                                              gmx_simd_float_t    &v1,
                                              gmx_simd_float_t    &v2,
                                              gmx_simd_float_t    &v3)
{
    assert((size_t)offset % 64 == 0);
    gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f<align>(base,
                                                                   gmx_simd_load_fi(offset),
                                                                   v0, v1, v2, v3);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_x86_avx_512f(const float *        base,
                                              const gmx_int32_t    offset[],
                                              gmx_simd_float_t    &v0,
                                              gmx_simd_float_t    &v1)
{
    assert((size_t)offset % 64 == 0);
    gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f<align>(base,
                                                                   gmx_simd_load_fi(offset),
                                                                   v0, v1);
}


static const int gmx_simd_best_pair_alignment_f_x86_avx_512f = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_f_x86_avx_512f(const float *        base,
                                               const gmx_int32_t    offset[],
                                               gmx_simd_float_t    &v0,
                                               gmx_simd_float_t    &v1,
                                               gmx_simd_float_t    &v2)
{
    __m512i simdoffset;

    assert((size_t)offset % 64 == 0);

    simdoffset = gmx_simd_load_fi(offset);

    if (align == 4)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm512_mullo_epi32(simdoffset, _mm512_set1_epi32(align));
    }

    v0  = _mm512_i32gather_ps(simdoffset, base,   sizeof(float));
    v1  = _mm512_i32gather_ps(simdoffset, base+1, sizeof(float));
    v2  = _mm512_i32gather_ps(simdoffset, base+2, sizeof(float));
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_f_x86_avx_512f(float *              base,
                                                 const gmx_int32_t    offset[],
                                                 gmx_simd_float_t     v0,
                                                 gmx_simd_float_t     v1,
                                                 gmx_simd_float_t     v2)
{
    __m512i simdoffset;

    assert((size_t)offset % 64 == 0);

    simdoffset = gmx_simd_load_fi(offset);

    if (align == 4)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm512_mullo_epi32(simdoffset, _mm512_set1_epi32(align));
    }

    _mm512_i32scatter_ps(base,   simdoffset, v0, sizeof(float));
    _mm512_i32scatter_ps(base+1, simdoffset, v1, sizeof(float));
    _mm512_i32scatter_ps(base+2, simdoffset, v2, sizeof(float));
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_f_x86_avx_512f(float *              base,
                                                const gmx_int32_t    offset[],
                                                gmx_simd_float_t     v0,
                                                gmx_simd_float_t     v1,
                                                gmx_simd_float_t     v2)
{
    __m512 t0, t1, t2;
    gmx_simd_gather_loadu_transpose_f_x86_avx_512f<align>(base, offset, t0, t1, t2);
    t0 = _mm512_add_ps(t0, v0);
    t1 = _mm512_add_ps(t1, v1);
    t2 = _mm512_add_ps(t2, v2);
    gmx_simd_transpose_scatter_storeu_f_x86_avx_512f<align>(base, offset, t0, t1, t2);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_f_x86_avx_512f(float *              base,
                                                const gmx_int32_t    offset[],
                                                gmx_simd_float_t     v0,
                                                gmx_simd_float_t     v1,
                                                gmx_simd_float_t     v2)
{
    __m512 t0, t1, t2;
    gmx_simd_gather_loadu_transpose_f_x86_avx_512f<align>(base, offset, t0, t1, t2);
    t0 = _mm512_sub_ps(t0, v0);
    t1 = _mm512_sub_ps(t1, v1);
    t2 = _mm512_sub_ps(t2, v2);
    gmx_simd_transpose_scatter_storeu_f_x86_avx_512f<align>(base, offset, t0, t1, t2);
}

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_f_x86_avx_512f(gmx_simd_float_t   scalar,
                                                   gmx_simd_float_t  &triplets0,
                                                   gmx_simd_float_t  &triplets1,
                                                   gmx_simd_float_t  &triplets2)
{
    triplets0 = _mm512_castsi512_ps(_mm512_permutevar_epi32(_mm512_set_epi32(5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, 0),
                                                            _mm512_castps_si512(scalar)));
    triplets1 = _mm512_castsi512_ps(_mm512_permutevar_epi32(_mm512_set_epi32(10, 10, 9, 9, 9, 8, 8, 8, 7, 7, 7, 6, 6, 6, 5, 5),
                                                            _mm512_castps_si512(scalar)));
    triplets2 = _mm512_castsi512_ps(_mm512_permutevar_epi32(_mm512_set_epi32(15, 15, 15, 14, 14, 14, 13, 13, 13, 12, 12, 12, 11, 11, 11, 10),
                                                            _mm512_castps_si512(scalar)));
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f(const float *       base,
                                                        gmx_simd_fint32_t   simdoffset,
                                                        gmx_simd_float_t   &v0,
                                                        gmx_simd_float_t   &v1,
                                                        gmx_simd_float_t   &v2,
                                                        gmx_simd_float_t   &v3)
{
    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);

    if (align == 4)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm512_mullo_epi32(simdoffset, _mm512_set1_epi32(align));
    }

    v0  = _mm512_i32gather_ps(simdoffset, base,   sizeof(float));
    v1  = _mm512_i32gather_ps(simdoffset, base+1, sizeof(float));
    v2  = _mm512_i32gather_ps(simdoffset, base+2, sizeof(float));
    v3  = _mm512_i32gather_ps(simdoffset, base+3, sizeof(float));
}

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_f_x86_avx_512f(const float *       base,
                                                         gmx_simd_fint32_t   simdoffset,
                                                         gmx_simd_float_t   &v0,
                                                         gmx_simd_float_t   &v1)
{
    if (align == 2)
    {
        v0  = _mm512_i32gather_ps(simdoffset, base,   align * sizeof(float));
        v1  = _mm512_i32gather_ps(simdoffset, base+1, align * sizeof(float));
    }
    else
    {
        if (align == 4)
        {
            simdoffset = _mm512_slli_epi32(simdoffset, 2);
        }
        else if (align == 8)
        {
            simdoffset = _mm512_slli_epi32(simdoffset, 3);
        }
        else
        {
            simdoffset = _mm512_mullo_epi32(simdoffset, _mm512_set1_epi32(align));
        }
        v0  = _mm512_i32gather_ps(simdoffset, base,   sizeof(float));
        v1  = _mm512_i32gather_ps(simdoffset, base+1, sizeof(float));
    }
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f(const float *       base,
                                                        gmx_simd_fint32_t   simdoffset,
                                                        gmx_simd_float_t   &v0,
                                                        gmx_simd_float_t   &v1)
{
    assert((size_t)base % 8 == 0);
    assert(align % 2 == 0);
    gmx_simd_gather_loadu_bysimdint_transpose_f_x86_avx_512f<align>(base, simdoffset, v0, v1);
}

static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_f_x86_avx_512f(float *           m,
                                                 gmx_simd_float_t  v0,
                                                 gmx_simd_float_t  v1,
                                                 gmx_simd_float_t  v2,
                                                 gmx_simd_float_t  v3)
{
    float  f;
    __m512 t0, t1, t2, t3;

    assert((size_t)m % 16 == 0);

    t0 = _mm512_add_ps(v0, _mm512_permute_ps(v0, 0b01001110));
    t0 = _mm512_mask_add_ps(t0, _mm512_int2mask(0xCCCC), v2, _mm512_permute_ps(v2, 0b01001110));
    t1 = _mm512_add_ps(v1, _mm512_permute_ps(v1, 0b01001110));
    t1 = _mm512_mask_add_ps(t1, _mm512_int2mask(0xCCCC), v3, _mm512_permute_ps(v3, 0b01001110));
    t2 = _mm512_add_ps(t0, _mm512_permute_ps(t0, 0b10110001));
    t2 = _mm512_mask_add_ps(t2, _mm512_int2mask(0xAAAA), t1, _mm512_permute_ps(t1, 0b10110001));

    t2 = _mm512_add_ps(t2, _mm512_shuffle_f32x4(t2, t2, 0b01001110));
    t2 = _mm512_add_ps(t2, _mm512_shuffle_f32x4(t2, t2, 0b10110001));

    t0 = _mm512_maskz_loadu_ps(_mm512_int2mask(0xF), m);
    t0 = _mm512_add_ps(t0, t2);
    _mm512_mask_storeu_ps(m, _mm512_int2mask(0xF), t0);

    t2 = _mm512_add_ps(t2, _mm512_permute_ps(t2, 0b01001110));
    t2 = _mm512_add_ps(t2, _mm512_permute_ps(t2, 0b10110001));

    _mm512_mask_storeu_ps(&f, _mm512_int2mask(0x0001), t2);
    return f;
}

/******************************************************
 * Single precision half-simd-width utility functions *
 ******************************************************/
static gmx_inline gmx_simd_float_t
gmx_simd_load_dual_hsimd_f_x86_avx_512f(const float * m0,
                                        const float * m1)
{
    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    return _mm512_mask_loadu_ps(_mm512_maskz_loadu_ps(_mm512_int2mask(0x00FF), m0), _mm512_int2mask(0xFF00), m1-8);
}

static gmx_inline gmx_simd_float_t
gmx_simd_loaddup_hsimd_f_x86_avx_512f(const float * m)
{
    assert((size_t)m % 32 == 0);

    __m512 tmp;

    tmp = _mm512_maskz_loadu_ps(_mm512_int2mask(0x00FF), m);
    return _mm512_shuffle_f32x4(tmp, tmp, 0b01000100);
}

static gmx_inline gmx_simd_float_t
gmx_simd_load1_dual_hsimd_f_x86_avx_512f(const float * m)
{
    __m512 tmp;
    tmp = _mm512_maskz_expandloadu_ps(_mm512_int2mask(0x0101), m);
    tmp = _mm512_permute_ps(tmp, 0b00000000);
    return _mm512_shuffle_f32x4(tmp, tmp, 0b10100000);
}


static gmx_inline void
gmx_simd_store_dual_hsimd_f_x86_avx_512f(float *           m0,
                                         float *           m1,
                                         gmx_simd_float_t  a)
{
    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    _mm512_mask_storeu_ps(m0,   _mm512_int2mask(0x00FF), a);
    _mm512_mask_storeu_ps(m1-8, _mm512_int2mask(0xFF00), a);
}

static gmx_inline void
gmx_simd_decr_hsimd_f_x86_avx_512f(float *          m,
                                   gmx_simd_float_t a)
{
    __m512 t;

    assert((size_t)m % 32 == 0);

    a = _mm512_add_ps(a, _mm512_shuffle_f32x4(a, a, 0b11101110));
    t = _mm512_maskz_loadu_ps(_mm512_int2mask(0x00FF), m);
    t = _mm512_sub_ps(t, a);
    _mm512_mask_storeu_ps(m, _mm512_int2mask(0x00FF), t);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_hsimd_f_x86_avx_512f(const float *       base0,
                                                    const float *       base1,
                                                    gmx_int32_t         offset[],
                                                    gmx_simd_float_t   &v0,
                                                    gmx_simd_float_t   &v1)
{
    __m512i idx0, idx1, idx;
    __m512  tmp1, tmp2;

    assert((size_t)offset % 32 == 0);
    assert((size_t)base0 % 8 == 0);
    assert((size_t)base1 % 8 == 0);
    assert((size_t)align % 2 == 0);

    idx0 = _mm512_maskz_loadu_epi32(_mm512_int2mask(0x00FF), offset);

    idx0 = _mm512_mullo_epi32(idx0, _mm512_set1_epi32(align));
    idx1 = _mm512_add_epi32(idx0, _mm512_set1_epi32(1));

    idx = _mm512_mask_shuffle_i32x4(idx0, _mm512_int2mask(0xFF00), idx1, idx1, 0b01000100);

    tmp1 = _mm512_i32gather_ps(idx, base0, sizeof(float));
    tmp2 = _mm512_i32gather_ps(idx, base1, sizeof(float));

    v0 = _mm512_shuffle_f32x4(tmp1, tmp2, 0b01000100 );
    v1 = _mm512_shuffle_f32x4(tmp1, tmp2, 0b11101110 );
}

static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_hsimd_f_x86_avx_512f(float *            m,
                                                       gmx_simd_float_t   v0,
                                                       gmx_simd_float_t   v1)
{
    float  f;
    __m512 t0, t1;

    assert((size_t)m % 32 == 0);

    /* This is not optimal, but no point optimizing until we know AVX-512 latencies */
    t0 = _mm512_add_ps(v0, _mm512_permute_ps(v0, 0b01001110));
    t1 = _mm512_add_ps(v1, _mm512_permute_ps(v1, 0b01001110));
    t0 = _mm512_add_ps(t0, _mm512_permute_ps(t0, 0b10110001));
    t0 = _mm512_mask_add_ps(t0, _mm512_int2mask(0xCCCC), t1, _mm512_permute_ps(t1, 0b10110001));
    t0 = _mm512_add_ps(t0, _mm512_shuffle_f32x4(t0, t0, 0b10110001));
    t0 = _mm512_mask_shuffle_f32x4(t0, _mm512_int2mask(0b1010101010101010), t0, t0, 0b11101110);

    t1 = _mm512_maskz_loadu_ps(_mm512_int2mask(0xF), m);
    t1 = _mm512_add_ps(t1, t0);
    _mm512_mask_storeu_ps(m, _mm512_int2mask(0xF), t1);

    t0 = _mm512_add_ps(t0, _mm512_permute_ps(t0, 0b01001110));
    t0 = _mm512_add_ps(t0, _mm512_permute_ps(t0, 0b10110001));

    _mm512_mask_storeu_ps(&f, _mm512_int2mask(0x0001), t0);
    return f;
}
#endif


/*********************************************************
 * SIMD DOUBLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 *********************************************************/


static gmx_inline __m512d
gmx_simd_set_exponent_d_x86_avx_512f(__m512d a)
{
    const __m512i expbias      = _mm512_set1_epi32(1023);
    __m512i       iexp         = gmx_simd_cvt_d2i(a);

    iexp = _mm512_permutevar_epi32(_mm512_set_epi32(7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 0, 0), iexp);
    iexp = _mm512_mask_slli_epi32(_mm512_setzero_epi32(), _mm512_int2mask(0xAAAA), _mm512_add_epi32(iexp, expbias), 20);
    return _mm512_castsi512_pd(iexp);
}

static gmx_inline double
gmx_simd_reduce_d_x86_avx_512f(__m512d a)
{
    double d;
    a = _mm512_add_pd(a, _mm512_shuffle_f64x2(a, a, 0b11101110));
    a = _mm512_add_pd(a, _mm512_shuffle_f64x2(a, a, 0b00010001));
    a = _mm512_add_pd(a, _mm512_permute_pd(a, 0b00000001));
    _mm512_mask_storeu_pd(&d, _mm512_int2mask(0x0001), a);
    return d;
}

static gmx_inline void
gmx_simd_cvt_f2dd_x86_avx_512f(__m512 f, __m512d * d0, __m512d * d1)
{
    *d0 = _mm512_cvtpslo_pd(f);
    *d1 = _mm512_cvtpslo_pd(_mm512_shuffle_f32x4(f, f, 0b11101110));
}

static gmx_inline __m512
gmx_simd_cvt_dd2f_x86_avx_512f(__m512d d0, __m512d d1)
{
    __m512 f0 = _mm512_cvtpd_pslo(d0);
    __m512 f1 = _mm512_cvtpd_pslo(d1);
    return _mm512_shuffle_f32x4(f0, f1, 0b01000100);
}


/****************************************************
 * Double precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

/* Declare functions implemented further down that are used here */
template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f(const double *       base,
                                                        gmx_simd_dint32_t    simdoffset,
                                                        gmx_simd_double_t   &v0,
                                                        gmx_simd_double_t   &v1,
                                                        gmx_simd_double_t   &v2,
                                                        gmx_simd_double_t   &v3);

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f(const double *       base,
                                                        gmx_simd_dint32_t    simdoffset,
                                                        gmx_simd_double_t   &v0,
                                                        gmx_simd_double_t   &v1);

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d_x86_avx_512f(const double *        base,
                                              const gmx_int32_t     offset[],
                                              gmx_simd_double_t    &v0,
                                              gmx_simd_double_t    &v1,
                                              gmx_simd_double_t    &v2,
                                              gmx_simd_double_t    &v3)
{
    assert((size_t)offset % 32 == 0);
    gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f<align>(base,
                                                                   gmx_simd_load_di(offset),
                                                                   v0, v1, v2, v3);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d_x86_avx_512f(const double *        base,
                                              const gmx_int32_t     offset[],
                                              gmx_simd_double_t    &v0,
                                              gmx_simd_double_t    &v1)
{
    assert((size_t)offset % 32 == 0);
    gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f<align>(base,
                                                                   gmx_simd_load_di(offset),
                                                                   v0, v1);
}

static const int gmx_simd_best_pair_alignment_d_x86_avx_512f = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_d_x86_avx_512f(const double *        base,
                                               const gmx_int32_t     offset[],
                                               gmx_simd_double_t    &v0,
                                               gmx_simd_double_t    &v1,
                                               gmx_simd_double_t    &v2)
{
    __m512i simdoffset;

    assert((size_t)offset % 32 == 0);

    simdoffset = gmx_simd_load_di(offset);

    if (align == 4)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm512_mullo_epi32(simdoffset, _mm512_set1_epi32(align));
    }

    v0  = _mm512_i32logather_pd(simdoffset, base,   sizeof(double));
    v1  = _mm512_i32logather_pd(simdoffset, base+1, sizeof(double));
    v2  = _mm512_i32logather_pd(simdoffset, base+2, sizeof(double));
}



template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_d_x86_avx_512f(double *              base,
                                                 const gmx_int32_t     offset[],
                                                 gmx_simd_double_t     v0,
                                                 gmx_simd_double_t     v1,
                                                 gmx_simd_double_t     v2)
{
    __m512i simdoffset;

    assert((size_t)offset % 32 == 0);

    simdoffset = gmx_simd_load_di(offset);

    if (align == 4)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm512_mullo_epi32(simdoffset, _mm512_set1_epi32(align));
    }

    _mm512_i32loscatter_pd(base,   simdoffset, v0, sizeof(double));
    _mm512_i32loscatter_pd(base+1, simdoffset, v1, sizeof(double));
    _mm512_i32loscatter_pd(base+2, simdoffset, v2, sizeof(double));
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_d_x86_avx_512f(double *              base,
                                                const gmx_int32_t     offset[],
                                                gmx_simd_double_t     v0,
                                                gmx_simd_double_t     v1,
                                                gmx_simd_double_t     v2)
{
    __m512d t0, t1, t2;
    gmx_simd_gather_loadu_transpose_d_x86_avx_512f<align>(base, offset, t0, t1, t2);
    t0 = _mm512_add_pd(t0, v0);
    t1 = _mm512_add_pd(t1, v1);
    t2 = _mm512_add_pd(t2, v2);
    gmx_simd_transpose_scatter_storeu_d_x86_avx_512f<align>(base, offset, t0, t1, t2);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_d_x86_avx_512f(double *              base,
                                                const gmx_int32_t     offset[],
                                                gmx_simd_double_t     v0,
                                                gmx_simd_double_t     v1,
                                                gmx_simd_double_t     v2)
{
    __m512d t0, t1, t2;
    gmx_simd_gather_loadu_transpose_d_x86_avx_512f<align>(base, offset, t0, t1, t2);
    t0 = _mm512_sub_pd(t0, v0);
    t1 = _mm512_sub_pd(t1, v1);
    t2 = _mm512_sub_pd(t2, v2);
    gmx_simd_transpose_scatter_storeu_d_x86_avx_512f<align>(base, offset, t0, t1, t2);
}


static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_d_x86_avx_512f(gmx_simd_double_t   scalar,
                                                   gmx_simd_double_t  &triplets0,
                                                   gmx_simd_double_t  &triplets1,
                                                   gmx_simd_double_t  &triplets2)
{
    triplets0 = _mm512_castsi512_pd(_mm512_permutevar_epi32(_mm512_set_epi32(5, 4, 5, 4, 3, 2, 3, 2, 3, 2, 1, 0, 1, 0, 1, 0),
                                                            _mm512_castpd_si512(scalar)));
    triplets1 = _mm512_castsi512_pd(_mm512_permutevar_epi32(_mm512_set_epi32(11, 10, 9, 8, 9, 8, 9, 8, 7, 6, 7, 6, 7, 6, 5, 4),
                                                            _mm512_castpd_si512(scalar)));
    triplets2 = _mm512_castsi512_pd(_mm512_permutevar_epi32(_mm512_set_epi32(15, 14, 15, 14, 15, 14, 13, 12, 13, 12, 13, 12, 11, 10, 11, 10),
                                                            _mm512_castpd_si512(scalar)));
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f(const double *       base,
                                                        gmx_simd_dint32_t    simdoffset,
                                                        gmx_simd_double_t   &v0,
                                                        gmx_simd_double_t   &v1,
                                                        gmx_simd_double_t   &v2,
                                                        gmx_simd_double_t   &v3)
{
    assert((size_t)base % 32 == 0);
    assert(align % 4 == 0);

    if (align == 4)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm512_mullo_epi32(simdoffset, _mm512_set1_epi32(align));
    }

    v0  = _mm512_i32logather_pd(simdoffset, base,   sizeof(double));
    v1  = _mm512_i32logather_pd(simdoffset, base+1, sizeof(double));
    v2  = _mm512_i32logather_pd(simdoffset, base+2, sizeof(double));
    v3  = _mm512_i32logather_pd(simdoffset, base+3, sizeof(double));
}


template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_d_x86_avx_512f(const double *       base,
                                                         gmx_simd_dint32_t    simdoffset,
                                                         gmx_simd_double_t   &v0,
                                                         gmx_simd_double_t   &v1)
{
    if (align == 2)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 1);
    }
    else if (align == 4)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm512_mullo_epi32(simdoffset, _mm512_set1_epi32(align));
    }

    v0  = _mm512_i32logather_pd(simdoffset, base,   sizeof(double));
    v1  = _mm512_i32logather_pd(simdoffset, base+1, sizeof(double));
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f(const double *       base,
                                                        gmx_simd_dint32_t    simdoffset,
                                                        gmx_simd_double_t   &v0,
                                                        gmx_simd_double_t   &v1)
{
    assert((size_t)base % 32 == 0);
    assert(align % 2 == 0);
    gmx_simd_gather_loadu_bysimdint_transpose_d_x86_avx_512f<align>(base, simdoffset, v0, v1);
}


static gmx_inline double gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_d_x86_avx_512f(double * m,
                                                 gmx_simd_double_t v0, gmx_simd_double_t v1,
                                                 gmx_simd_double_t v2, gmx_simd_double_t v3)
{
    double  d;
    __m512d t0, t2;

    assert((size_t)m % 32 == 0);

    /* Not optimal. Will optimize when we know more about AVX-512 latencies */
    t0 = _mm512_add_pd(v0, _mm512_permute_pd(v0, 0b01010101));
    t2 = _mm512_add_pd(v2, _mm512_permute_pd(v2, 0b01010101));
    t0 = _mm512_mask_add_pd(t0, _mm512_int2mask(0b10101010), v1, _mm512_permute_pd(v1, 0b01010101));
    t2 = _mm512_mask_add_pd(t2, _mm512_int2mask(0b10101010), v3, _mm512_permute_pd(v3, 0b01010101));
    t0 = _mm512_add_pd(t0, _mm512_shuffle_f64x2(t0, t0, 0b01001110));
    t0 = _mm512_mask_add_pd(t0, _mm512_int2mask(0xF0), t2, _mm512_shuffle_f64x2(t2, t2, 0b01001110));
    t0 = _mm512_add_pd(t0, _mm512_shuffle_f64x2(t0, t0, 0b10110001));
    t0 = _mm512_mask_shuffle_f64x2(t0, _mm512_int2mask(0x0C), t0, t0, 0b11101110);

    t2 = _mm512_maskz_loadu_pd(_mm512_int2mask(0xF), m);
    t2 = _mm512_add_pd(t2, t0);
    _mm512_mask_storeu_pd(m, _mm512_int2mask(0xF), t2);

    t0 = _mm512_add_pd(t0, _mm512_permutex_pd(t0, 0b01001110));
    t0 = _mm512_add_pd(t0, _mm512_permutex_pd(t0, 0b10110001));

    _mm512_mask_storeu_pd(&d, _mm512_int2mask(0x0001), t0);
    return d;
}

/******************************************************
 * Double precision half-simd-width utility functions *
 ******************************************************/
static gmx_inline gmx_simd_double_t
gmx_simd_load_dual_hsimd_d_x86_avx_512f(const double * m0,
                                        const double * m1)
{
    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    return _mm512_mask_loadu_pd(_mm512_maskz_loadu_pd(_mm512_int2mask(0x0F), m0), _mm512_int2mask(0xF0), m1-4);
}

static gmx_inline gmx_simd_double_t
gmx_simd_loaddup_hsimd_d_x86_avx_512f(const double * m)
{
    __m512d tmp;

    assert((size_t)m % 32 == 0);

    tmp = _mm512_maskz_loadu_pd(_mm512_int2mask(0x0F), m);
    return _mm512_shuffle_f64x2(tmp, tmp, 0b01000100);
}

static gmx_inline gmx_simd_double_t
gmx_simd_load1_dual_hsimd_d_x86_avx_512f(const double * m)
{
    __m512d tmp;
    tmp = _mm512_maskz_expandloadu_pd(_mm512_int2mask(0x11), m);
    return _mm512_permutex_pd(tmp, 0);
}


static gmx_inline void
gmx_simd_store_dual_hsimd_d_x86_avx_512f(double *           m0,
                                         double *           m1,
                                         gmx_simd_double_t  a)
{
    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    _mm512_mask_storeu_pd(m0,   _mm512_int2mask(0x0F), a);
    _mm512_mask_storeu_pd(m1-4, _mm512_int2mask(0xF0), a);
}

static gmx_inline void
gmx_simd_decr_hsimd_d_x86_avx_512f(double *          m,
                                   gmx_simd_double_t a)
{
    __m512d t;

    assert((size_t)m % 32 == 0);

    a = _mm512_add_pd(a, _mm512_shuffle_f64x2(a, a, 0b11101110));
    t = _mm512_maskz_loadu_pd(_mm512_int2mask(0x0F), m);
    t = _mm512_sub_pd(t, a);
    _mm512_mask_storeu_pd(m, _mm512_int2mask(0x0F), t);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_hsimd_d_x86_avx_512f(const double *       base0,
                                                    const double *       base1,
                                                    gmx_int32_t          offset[],
                                                    gmx_simd_double_t   &v0,
                                                    gmx_simd_double_t   &v1)
{
    __m512i  idx0, idx1, idx;
    __m512d  tmp1, tmp2;

    assert((size_t)offset % 16 == 0);
    assert((size_t)base0 % 16 == 0);
    assert((size_t)base1 % 16 == 0);
    assert((size_t)align % 2  == 0);

    idx0 = _mm512_maskz_loadu_epi32(_mm512_int2mask(0x000F), offset);

    idx0 = _mm512_mullo_epi32(idx0, _mm512_set1_epi32(align));
    idx1 = _mm512_add_epi32(idx0, _mm512_set1_epi32(1));

    idx = _mm512_mask_shuffle_i32x4(idx0, _mm512_int2mask(0x00F0), idx1, idx1, 0b00000000);

    tmp1 = _mm512_i32logather_pd(idx, base0, sizeof(double));
    tmp2 = _mm512_i32logather_pd(idx, base1, sizeof(double));

    v0 = _mm512_shuffle_f64x2(tmp1, tmp2, 0b01000100 );
    v1 = _mm512_shuffle_f64x2(tmp1, tmp2, 0b11101110 );
}

static gmx_inline double
gmx_simd_reduce_incr_4_return_sum_hsimd_d_x86_avx_512f(double *           m,
                                                       gmx_simd_double_t  v0,
                                                       gmx_simd_double_t  v1)
{
    double   d;
    __m512d  t0, t1;

    assert((size_t)m % 32 == 0);

    t0 = _mm512_add_pd(v0, _mm512_permutex_pd(v0, 0b01001110));
    t0 = _mm512_mask_add_pd(t0, _mm512_int2mask(0xCC), v1, _mm512_permutex_pd(v1, 0b01001110));
    t0 = _mm512_add_pd(t0, _mm512_permutex_pd(t0, 0b10110001));
    t0 = _mm512_mask_shuffle_f64x2(t0, _mm512_int2mask(0b10101010), t0, t0, 0b11101110);

    t1 = _mm512_maskz_loadu_pd(_mm512_int2mask(0xF), m);
    t1 = _mm512_add_pd(t1, t0);
    _mm512_mask_storeu_pd(m, _mm512_int2mask(0xF), t1);

    t0 = _mm512_add_pd(t0, _mm512_permutex_pd(t0, 0b01001110));
    t0 = _mm512_add_pd(t0, _mm512_permutex_pd(t0, 0b10110001));

    _mm512_mask_storeu_pd(&d, _mm512_int2mask(0x01), t0);
    return d;
}
#endif


/**********************************************************
 * SIMD4 SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 **********************************************************/

static gmx_inline float
gmx_simd4_reduce_f_x86_avx_512f(__m512 a)
{
    float f;
    a = _mm512_add_ps(a, _mm512_permute_ps(a, 0b00111001));
    a = _mm512_add_ps(a, _mm512_permute_ps(a, 0b01001110));
    _mm512_mask_storeu_ps(&f, _mm512_int2mask(0x0001), a);
    return f;
}

static gmx_inline float
gmx_simd4_dotproduct3_f_x86_avx_512f(__m512 a, __m512 b)
{
    float  f;
    __m512 c;
    a = _mm512_mul_ps(a, b);
    c = _mm512_add_ps(a, _mm512_permute_ps(a, 0b00111001));
    c = _mm512_add_ps(c, _mm512_permute_ps(a, 0b01001110));
    _mm512_mask_storeu_ps(&f, _mm512_int2mask(0x0001), c);
    return f;
}

#ifdef __cplusplus
static gmx_inline void gmx_simdcall
gmx_simd4_transpose_f_x86_avx_512f(gmx_simd4_float_t &v0, gmx_simd4_float_t &v1,
                                   gmx_simd4_float_t &v2, gmx_simd4_float_t &v3)
{
    __m512 t0, t1, t2, t3;

    t0 = _mm512_unpacklo_ps(v0, v2);
    t1 = _mm512_unpackhi_ps(v0, v2);
    t2 = _mm512_unpacklo_ps(v1, v3);
    t3 = _mm512_unpackhi_ps(v1, v3);
    v0 = _mm512_unpacklo_ps(t0, t2);
    v1 = _mm512_unpackhi_ps(t0, t2);
    v2 = _mm512_unpacklo_ps(t1, t3);
    v3 = _mm512_unpackhi_ps(t1, t3);
}
#endif

/**********************************************************
 * SIMD4 DOUBLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 **********************************************************/

static gmx_inline double
gmx_simd4_reduce_d_x86_avx_512f(__m512d a)
{
    double d;
    a = _mm512_add_pd(a, _mm512_permute_pd(a, 0b01010101));
    a = _mm512_add_pd(a, _mm512_shuffle_f64x2(a, a, 0b00010001));
    _mm512_mask_storeu_pd(&d, _mm512_int2mask(0x01), a);
    return d;
}

static gmx_inline double
gmx_simd4_dotproduct3_d_x86_avx_512f(__m512d a, __m512d b)
{
    double  d;
    __m512d c;
    a = _mm512_mul_pd(a, b);
    c = _mm512_add_pd(a, _mm512_permutex_pd(a, 0b00111001));
    c = _mm512_add_pd(c, _mm512_permutex_pd(a, 0b01001110));
    _mm512_mask_storeu_pd(&d, _mm512_int2mask(0x01), c);
    return d;
}

#ifdef __cplusplus
static gmx_inline void gmx_simdcall
gmx_simd4_transpose_d_x86_avx_512f(gmx_simd4_double_t &v0, gmx_simd4_double_t &v1,
                                   gmx_simd4_double_t &v2, gmx_simd4_double_t &v3)
{
    __m512d t0, t1, t2, t3;

    t0 = _mm512_unpacklo_pd(v0, v1);
    t1 = _mm512_unpackhi_pd(v0, v1);
    t2 = _mm512_unpacklo_pd(v2, v3);
    t3 = _mm512_unpackhi_pd(v2, v3);

    v0 = _mm512_mask_permutex_pd(t0, _mm512_int2mask(0x0C), t2, 0b01000100);
    v1 = _mm512_mask_permutex_pd(t1, _mm512_int2mask(0x0C), t3, 0b01000100);
    v2 = _mm512_mask_permutex_pd(t2, _mm512_int2mask(0x03), t0, 0b11101110);
    v3 = _mm512_mask_permutex_pd(t3, _mm512_int2mask(0x03), t1, 0b11101110);
}
#endif



/* Function to check whether SIMD operations have resulted in overflow */
static int
gmx_simd_check_and_reset_overflow(void)
{
    int                MXCSR;
    int                sse_overflow;
    /* The overflow flag is bit 3 in the register */
    const unsigned int flag = 0x8;

    MXCSR = _mm_getcsr();
    if (MXCSR & flag)
    {
        sse_overflow = 1;
        /* Set the overflow flag to zero */
        MXCSR = MXCSR & ~flag;
        _mm_setcsr(MXCSR);
    }
    else
    {
        sse_overflow = 0;
    }
    return sse_overflow;
}

static void
gmx_simd_prefetch(const void * m)
{
    _mm_prefetch((const char *)m, _MM_HINT_T0);
}

#endif /* GMX_SIMD_IMPL_X86_AVX_512F_H */
