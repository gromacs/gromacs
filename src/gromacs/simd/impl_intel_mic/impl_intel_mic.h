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

#ifndef GMX_SIMD_IMPL_INTEL_MIC_H
#define GMX_SIMD_IMPL_INTEL_MIC_H

#include "config.h"

#include <assert.h>
#include <math.h>

#include <immintrin.h>
/* Intel Xeon Phi, or
 * the-artist-formerly-known-as-Knight's-corner, or
 * the-artist-formerly-formerly-known-as-MIC, or
 * the artist formerly-formerly-formerly-known-as-Larrabee
 * 512-bit SIMD instruction wrappers.
 */

#define GMX_SIMD_V2

/* Capability definitions for Xeon Phi SIMD */
#define GMX_SIMD_HAVE_FLOAT
#define GMX_SIMD_HAVE_DOUBLE
#define GMX_SIMD_HAVE_SIMD_HARDWARE
#define GMX_SIMD_HAVE_LOADU
#define GMX_SIMD_HAVE_STOREU
#define GMX_SIMD_HAVE_LOGICAL
#define GMX_SIMD_HAVE_FMA
#undef  GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD_HAVE_FINT32
#define GMX_SIMD_HAVE_FINT32_EXTRACT
#define GMX_SIMD_HAVE_FINT32_LOGICAL
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS
#define GMX_SIMD_HAVE_DINT32
#define GMX_SIMD_HAVE_DINT32_EXTRACT
#define GMX_SIMD_HAVE_DINT32_LOGICAL
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS
#define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT
#define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_DOUBLE
#define GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT
#define GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE
#define GMX_SIMD4_HAVE_FLOAT
#define GMX_SIMD4_HAVE_DOUBLE

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH        16
#define GMX_SIMD_DOUBLE_WIDTH        8
#define GMX_SIMD_FINT32_WIDTH       16
#define GMX_SIMD_DINT32_WIDTH        8
#define GMX_SIMD_RSQRT_BITS         23
#define GMX_SIMD_RCP_BITS           23

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
#define gmx_simd_mul_mask_f(a, b, m)       _mm512_mask_mul_ps(_mm512_setzero_ps(), m, a, b)
#define gmx_simd_fmadd_mask_f(a, b, c, m)  _mm512_mask_mov_ps(_mm512_setzero_ps(), m, _mm512_fmadd_ps(a, b, c))
#define gmx_simd_rsqrt_mask_f(a, m)        _mm512_mask_rsqrt23_ps(_mm512_setzero_ps(), m, a)
#define gmx_simd_rcp_mask_f(a, m)          _mm512_mask_rcp23_ps(_mm512_setzero_ps(), m, a)
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
#define gmx_simd_cmpnz_f(a)        _mm512_test_epi32_mask( _mm512_castps_si512(a), _mm512_castps_si512(a) )
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
/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_f             gmx_simd_gather_load_transpose_f_mic
#define gmx_simd_best_pair_alignment_f               gmx_simd_best_pair_alignment_f_mic
#define gmx_simd_gather_loadu_transpose_f            gmx_simd_gather_loadu_transpose_f_mic
#define gmx_simd_transpose_scatter_storeu_f          gmx_simd_transpose_scatter_storeu_f_mic
#define gmx_simd_transpose_scatter_incru_f           gmx_simd_transpose_scatter_incru_f_mic
#define gmx_simd_transpose_scatter_decru_f           gmx_simd_transpose_scatter_decru_f_mic
#define gmx_simd_expand_scalars_to_triplets_f        gmx_simd_expand_scalars_to_triplets_f_mic
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_mic
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_mic
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_mic
/* Half-simd-width utilities */
#define gmx_simd_load_dual_hsimd_f                   gmx_simd_load_dual_hsimd_f_mic
#define gmx_simd_loaddup_hsimd_f                     gmx_simd_loaddup_hsimd_f_mic
#define gmx_simd_load1_dual_hsimd_f                  gmx_simd_load1_dual_hsimd_f_mic
#define gmx_simd_store_dual_hsimd_f                  gmx_simd_store_dual_hsimd_f_mic
#define gmx_simd_decr_hsimd_f                        gmx_simd_decr_hsimd_f_mic
#define gmx_simd_gather_load_transpose_hsimd_f       gmx_simd_gather_load_transpose_hsimd_f_mic
#define gmx_simd_reduce_incr_4_return_sum_hsimd_f    gmx_simd_reduce_incr_4_return_sum_hsimd_f_mic

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_double_t          __m512d
#define gmx_simd_load_d            _mm512_load_pd
#define gmx_simd_load1_d(m)        _mm512_extload_pd(m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE)
#define gmx_simd_set1_d            _mm512_set1_pd
#define gmx_simd_store_d           _mm512_store_pd
#define gmx_simd_loadu_d           gmx_simd_loadu_d_mic
#define gmx_simd_storeu_d          gmx_simd_storeu_d_mic
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
#define gmx_simd_rsqrt_d(x)        _mm512_cvtpslo_pd(_mm512_rsqrt23_ps(_mm512_cvtpd_pslo(x)))
#define gmx_simd_rcp_d(x)          _mm512_cvtpslo_pd(_mm512_rcp23_ps(_mm512_cvtpd_pslo(x)))
#define gmx_simd_mul_mask_d(a, b, m)       _mm512_mask_mul_pd(_mm512_setzero_pd(), m, a, b)
#define gmx_simd_fmadd_mask_d(a, b, c, m)  _mm512_mask_mov_pd(_mm512_setzero_pd(), m, _mm512_fmadd_pd(a, b, c))
#define gmx_simd_rsqrt_mask_d(a, m)        _mm512_cvtpslo_pd(_mm512_mask_rsqrt23_ps(_mm512_setzero_ps(), m, _mm512_cvtpd_pslo(a)))
#define gmx_simd_rcp_mask_d(a, m)          _mm512_cvtpslo_pd(_mm512_mask_rcp23_ps(_mm512_setzero_ps(), m, _mm512_cvtpd_pslo(a)))
#define gmx_simd_fabs_d(x)         gmx_simd_andnot_d(_mm512_set1_pd(GMX_DOUBLE_NEGZERO), x)
#define gmx_simd_fneg_d(x)         _mm512_addn_pd(x, _mm512_setzero_pd())
#define gmx_simd_max_d             _mm512_gmax_pd
#define gmx_simd_min_d             _mm512_gmin_pd
#define gmx_simd_round_d(a)        _mm512_roundfxpnt_adjust_pd(a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define gmx_simd_trunc_d(a)        _mm512_roundfxpnt_adjust_pd(a, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define gmx_simd_fraction_d(x)     _mm512_sub_pd(x, gmx_simd_trunc_d(x))
#define gmx_simd_get_exponent_d(x) _mm512_getexp_pd(x)
#define gmx_simd_get_mantissa_d(x) _mm512_getmant_pd(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero)
#define gmx_simd_set_exponent_d(x) gmx_simd_set_exponent_d_mic(x)
/* integer datatype corresponding to float: gmx_simd_fint32_t
   Doesn't use mask other than where required. No side effect expected for operating on the (unused) upper 8.
 */
#define gmx_simd_dint32_t          __m512i
#define gmx_simd_load_di(m)        _mm512_extload_epi64(m, _MM_UPCONV_EPI64_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE)
#define gmx_simd_set1_di           _mm512_set1_epi32
#define gmx_simd_store_di          gmx_simd_store_di_mic
#define gmx_simd_loadu_di          gmx_simd_loadu_di_mic
#define gmx_simd_storeu_di         gmx_simd_storeu_di_mic
#define gmx_simd_extract_di        gmx_simd_extract_di_mic
#define gmx_simd_setzero_di        _mm512_setzero_epi32
#define gmx_simd_cvt_d2i(a)        _mm512_cvtfxpnt_roundpd_epi32lo(a, _MM_FROUND_TO_NEAREST_INT)
#define gmx_simd_cvtt_d2i(a)       _mm512_cvtfxpnt_roundpd_epi32lo(a, _MM_FROUND_TO_ZERO)
#define gmx_simd_cvt_i2d           _mm512_cvtepi32lo_pd
/* Integer logical ops on gmx_simd_fint32_t */
#define gmx_simd_slli_di           _mm512_slli_epi32
#define gmx_simd_srli_di           _mm512_srli_epi32
#define gmx_simd_and_di            _mm512_and_epi32
#define gmx_simd_andnot_di         _mm512_andnot_epi32
#define gmx_simd_or_di             _mm512_or_epi32
#define gmx_simd_xor_di            _mm512_xor_epi32
/* Integer arithmetic ops on gmx_simd_fint32_t */
#define gmx_simd_add_di            _mm512_add_epi32
#define gmx_simd_sub_di            _mm512_sub_epi32
#define gmx_simd_mul_di            _mm512_mullo_epi32
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_dbool_t           __mmask8
#define gmx_simd_cmpeq_d(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_EQ_OQ)
#define gmx_simd_cmpnz_d(a)        _mm512_cmp_pd_mask(a, _mm512_setzero_pd(), _CMP_NEQ_UQ)
#define gmx_simd_cmplt_d(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_LT_OS)
#define gmx_simd_cmple_d(a, b)     _mm512_cmp_pd_mask(a, b, _CMP_LE_OS)
#define gmx_simd_and_db            _mm512_kand
#define gmx_simd_or_db             _mm512_kor
#define gmx_simd_anytrue_db(x)     _mm512_mask2int(x)
#define gmx_simd_blendzero_d(a, sel)    _mm512_mask_mov_pd(_mm512_setzero_pd(), sel, a)
#define gmx_simd_blendnotzero_d(a, sel) _mm512_mask_mov_pd(_mm512_setzero_pd(), _mm512_knot(sel), a)
#define gmx_simd_blendv_d(a, b, sel)    _mm512_mask_blend_pd(sel, a, b)
#define gmx_simd_reduce_d(a)       _mm512_reduce_add_pd(a)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_dibool_t          __mmask16
#define gmx_simd_cmpeq_di(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_EQ)
#define gmx_simd_cmplt_di(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_LT)
#define gmx_simd_and_dib           _mm512_kand
#define gmx_simd_or_dib            _mm512_kor
#define gmx_simd_anytrue_dib(x)    (_mm512_mask2int(x)&0xFF)
#define gmx_simd_blendzero_di(a, sel)    _mm512_mask_mov_epi32(_mm512_setzero_epi32(), sel, a)
#define gmx_simd_blendnotzero_di(a, sel) _mm512_mask_mov_epi32(_mm512_setzero_epi32(), _mm512_knot(sel), a)
#define gmx_simd_blendv_di(a, b, sel)    _mm512_mask_blend_epi32(sel, a, b)
/* Conversions between booleans. Double & dint stuff is stored in low bits */
#define gmx_simd_cvt_db2dib(x)     (x)
#define gmx_simd_cvt_dib2db(x)     (x)
/* Float/double conversion */
#define gmx_simd_cvt_f2dd          gmx_simd_cvt_f2dd_mic
#define gmx_simd_cvt_dd2f          gmx_simd_cvt_dd2f_mic
/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_d             gmx_simd_gather_load_transpose_d_mic
#define gmx_simd_best_pair_alignment_d               gmx_simd_best_pair_alignment_d_mic
#define gmx_simd_gather_loadu_transpose_d            gmx_simd_gather_loadu_transpose_d_mic
#define gmx_simd_transpose_scatter_storeu_d          gmx_simd_transpose_scatter_storeu_d_mic
#define gmx_simd_transpose_scatter_incru_d           gmx_simd_transpose_scatter_incru_d_mic
#define gmx_simd_transpose_scatter_decru_d           gmx_simd_transpose_scatter_decru_d_mic
#define gmx_simd_expand_scalars_to_triplets_d        gmx_simd_expand_scalars_to_triplets_d_mic
#define gmx_simd_gather_load_bysimdint_transpose_d   gmx_simd_gather_load_bysimdint_transpose_d_mic
#define gmx_simd_gather_loadu_bysimdint_transpose_d  gmx_simd_gather_loadu_bysimdint_transpose_d_mic
#define gmx_simd_reduce_incr_4_return_sum_d          gmx_simd_reduce_incr_4_return_sum_d_mic
/* Half-simd-width utilities */
#define gmx_simd_load_dual_hsimd_d                   gmx_simd_load_dual_hsimd_d_mic
#define gmx_simd_loaddup_hsimd_d                     gmx_simd_loaddup_hsimd_d_mic
#define gmx_simd_load1_dual_hsimd_d                  gmx_simd_load1_dual_hsimd_d_mic
#define gmx_simd_store_dual_hsimd_d                  gmx_simd_store_dual_hsimd_d_mic
#define gmx_simd_decr_hsimd_d                        gmx_simd_decr_hsimd_d_mic
#define gmx_simd_gather_load_transpose_hsimd_d       gmx_simd_gather_load_transpose_hsimd_d_mic
#define gmx_simd_reduce_incr_4_return_sum_hsimd_d    gmx_simd_reduce_incr_4_return_sum_hsimd_d_mic


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
#define gmx_simd4_dotproduct3_f(a, b) gmx_simd4_dotproduct3_f_mic(a, b)
#define gmx_simd4_transpose_f       gmx_simd4_transpose_f_mic
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
#define gmx_simd4_reduce_f(x)        gmx_simd4_reduce_f_mic(x)

/****************************************************
 *      DOUBLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
#define gmx_simd4_double_t          __m512d
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
#define gmx_simd4_and_d(a, b)       _mm512_castsi512_pd(_mm512_mask_and_epi32(_mm512_undefined_epi32(), mask_lo, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd4_andnot_d(a, b)    _mm512_castsi512_pd(_mm512_mask_andnot_epi32(_mm512_undefined_epi32(), mask_lo, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd4_or_d(a, b)        _mm512_castsi512_pd(_mm512_mask_or_epi32(_mm512_undefined_epi32(), mask_lo, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd4_xor_d(a, b)       _mm512_castsi512_pd(_mm512_mask_xor_epi32(_mm512_undefined_epi32(), mask_lo, _mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define gmx_simd4_rsqrt_d(a)        _mm512_mask_cvtpslo_pd(_mm512_undefined_pd(), gmx_simd4_mask, _mm512_mask_rsqrt23_ps(_mm512_undefined_ps(), gmx_simd4_mask, _mm512_mask_cvtpd_pslo(_mm512_undefined_ps(), gmx_simd4_mask, x)))
#define gmx_simd4_fabs_d(x)         gmx_simd4_andnot_d(_mm512_set1_pd(GMX_DOUBLE_NEGZERO), x)
#define gmx_simd4_fneg_d(x)         _mm512_mask_addn_pd(_mm512_undefined_pd(), gmx_simd4_mask, x, _mm512_setzero_pd())
#define gmx_simd4_max_d(a, b)       _mm512_mask_gmax_pd(_mm512_undefined_pd(), gmx_simd4_mask, a, b)
#define gmx_simd4_min_d(a, b)       _mm512_mask_gmin_pd(_mm512_undefined_pd(), gmx_simd4_mask, a, b)
#define gmx_simd4_round_d(a)        _mm512_mask_roundfxpnt_adjust_pd(_mm512_undefined_pd(), gmx_simd4_mask, a, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
#define gmx_simd4_trunc_d(a)        _mm512_mask_roundfxpnt_adjust_pd(_mm512_undefined_pd(), gmx_simd4_mask, a, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
#define gmx_simd4_dotproduct3_d(a, b) gmx_simd4_dotproduct3_d_mic(a, b)
#define gmx_simd4_transpose_d       gmx_simd4_transpose_d_mic
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
#define gmx_simd4_reduce_d(x)       gmx_simd4_reduce_d_mic(x)

#define PERM_LOW2HIGH _MM_PERM_BABA
#define PERM_HIGH2LOW _MM_PERM_DCDC

#define mask_lo _mm512_int2mask(0x00FF) /* would be better a constant - but can't initialize with a function call. */
#define mask_hi _mm512_int2mask(0xFF00)


/* Declare the simd4 store routines implmented further down so they can be used */
static gmx_inline void gmx_simdcall
gmx_simd4_store_f_mic(float * m, __m512 s);

static gmx_inline void gmx_simdcall
gmx_simd4_store_d_mic(double * m, __m512d s);





/*********************************************************
 * SIMD SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 *********************************************************/

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


/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

/* Declare functions implemented further down that are used here */
template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_mic(const float *       base,
                                               gmx_simd_fint32_t   simdoffset,
                                               gmx_simd_float_t   &v0,
                                               gmx_simd_float_t   &v1,
                                               gmx_simd_float_t   &v2,
                                               gmx_simd_float_t   &v3);

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_mic(const float *       base,
                                               gmx_simd_fint32_t   simdoffset,
                                               gmx_simd_float_t   &v0,
                                               gmx_simd_float_t   &v1);


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_mic(const float *        base,
                                     const gmx_int32_t    offset[],
                                     gmx_simd_float_t    &v0,
                                     gmx_simd_float_t    &v1,
                                     gmx_simd_float_t    &v2,
                                     gmx_simd_float_t    &v3)
{
    assert((size_t)offset % 64 == 0);
    gmx_simd_gather_load_bysimdint_transpose_f_mic<align>(base,
                                                          gmx_simd_load_fi(offset),
                                                          v0, v1, v2, v3);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_mic(const float *        base,
                                     const gmx_int32_t    offset[],
                                     gmx_simd_float_t    &v0,
                                     gmx_simd_float_t    &v1)
{
    assert((size_t)offset % 64 == 0);
    gmx_simd_gather_load_bysimdint_transpose_f_mic<align>(base,
                                                          gmx_simd_load_fi(offset),
                                                          v0, v1);
}


static const int gmx_simd_best_pair_alignment_f_mic = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_f_mic(const float *        base,
                                      const gmx_int32_t    offset[],
                                      gmx_simd_float_t    &v0,
                                      gmx_simd_float_t    &v1,
                                      gmx_simd_float_t    &v2)
{
    __m512i simdoffset;

    assert((size_t)offset % 64 == 0);

    simdoffset = gmx_simd_load_fi(offset);

    /* All instructions might be latency ~4 on MIC, so we use shifts where we
     * only need a single instruction (since the shift parameter is an immediate),
     * but multiplication otherwise.
     */
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
gmx_simd_transpose_scatter_storeu_f_mic(float *              base,
                                        const gmx_int32_t    offset[],
                                        gmx_simd_float_t     v0,
                                        gmx_simd_float_t     v1,
                                        gmx_simd_float_t     v2)
{
    __m512i simdoffset;

    assert((size_t)offset % 64 == 0);

    simdoffset = gmx_simd_load_fi(offset);

    /* All instructions might be latency ~4 on MIC, so we use shifts where we
     * only need a single instruction (since the shift parameter is an immediate),
     * but multiplication otherwise.
     */
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
gmx_simd_transpose_scatter_incru_f_mic(float *              base,
                                       const gmx_int32_t    offset[],
                                       gmx_simd_float_t     v0,
                                       gmx_simd_float_t     v1,
                                       gmx_simd_float_t     v2)
{
    __m512 t0, t1, t2;
    gmx_simd_gather_loadu_transpose_f_mic<align>(base, offset, t0, t1, t2);
    t0 = _mm512_add_ps(t0, v0);
    t1 = _mm512_add_ps(t1, v1);
    t2 = _mm512_add_ps(t2, v2);
    gmx_simd_transpose_scatter_storeu_f_mic<align>(base, offset, t0, t1, t2);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_f_mic(float *              base,
                                       const gmx_int32_t    offset[],
                                       gmx_simd_float_t     v0,
                                       gmx_simd_float_t     v1,
                                       gmx_simd_float_t     v2)
{
    __m512 t0, t1, t2;
    gmx_simd_gather_loadu_transpose_f_mic<align>(base, offset, t0, t1, t2);
    t0 = _mm512_sub_ps(t0, v0);
    t1 = _mm512_sub_ps(t1, v1);
    t2 = _mm512_sub_ps(t2, v2);
    gmx_simd_transpose_scatter_storeu_f_mic<align>(base, offset, t0, t1, t2);
}

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_f_mic(gmx_simd_float_t   scalar,
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
gmx_simd_gather_load_bysimdint_transpose_f_mic(const float *       base,
                                               gmx_simd_fint32_t   simdoffset,
                                               gmx_simd_float_t   &v0,
                                               gmx_simd_float_t   &v1,
                                               gmx_simd_float_t   &v2,
                                               gmx_simd_float_t   &v3)
{
    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);

    /* All instructions might be latency ~4 on MIC, so we use shifts where we
     * only need a single instruction (since the shift parameter is an immediate),
     * but multiplication otherwise.
     */
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
gmx_simd_gather_loadu_bysimdint_transpose_f_mic(const float *       base,
                                                gmx_simd_fint32_t   simdoffset,
                                                gmx_simd_float_t   &v0,
                                                gmx_simd_float_t   &v1)
{
    /* All instructions might be latency ~4 on MIC, so we use shifts where we
     * only need a single instruction (since the shift parameter is an immediate),
     * but multiplication otherwise.
     * For align == 2 we can merge the constant into the scale parameter,
     * which can take constants up to 8 in total.
     */
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
gmx_simd_gather_load_bysimdint_transpose_f_mic(const float *       base,
                                               gmx_simd_fint32_t   simdoffset,
                                               gmx_simd_float_t   &v0,
                                               gmx_simd_float_t   &v1)
{
    assert((size_t)base % 8 == 0);
    assert(align % 2 == 0);
    gmx_simd_gather_loadu_bysimdint_transpose_f_mic<align>(base, simdoffset, v0, v1);
}

static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_f_mic(float *           m,
                                        gmx_simd_float_t  v0,
                                        gmx_simd_float_t  v1,
                                        gmx_simd_float_t  v2,
                                        gmx_simd_float_t  v3)
{
    float  f;
    __m512 t0, t1, t2, t3;

    assert((size_t)m % 16 == 0);

    t0 = _mm512_add_ps(v0, _mm512_swizzle_ps(v0, _MM_SWIZ_REG_BADC));
    t0 = _mm512_mask_add_ps(t0, _mm512_int2mask(0xCCCC), v2, _mm512_swizzle_ps(v2, _MM_SWIZ_REG_BADC));
    t1 = _mm512_add_ps(v1, _mm512_swizzle_ps(v1, _MM_SWIZ_REG_BADC));
    t1 = _mm512_mask_add_ps(t1, _mm512_int2mask(0xCCCC), v3, _mm512_swizzle_ps(v3, _MM_SWIZ_REG_BADC));
    t2 = _mm512_add_ps(t0, _mm512_swizzle_ps(t0, _MM_SWIZ_REG_CDAB));
    t2 = _mm512_mask_add_ps(t2, _mm512_int2mask(0xAAAA), t1, _mm512_swizzle_ps(t1, _MM_SWIZ_REG_CDAB));

    t2 = _mm512_add_ps(t2, _mm512_permute4f128_ps(t2, _MM_PERM_BADC));
    t2 = _mm512_add_ps(t2, _mm512_permute4f128_ps(t2, _MM_PERM_CDAB));

    t0 = gmx_simd4_load_f(m);
    t0 = gmx_simd4_add_f(t0, t2);
    gmx_simd4_store_f(m, t0);

    t2 = _mm512_add_ps(t2, _mm512_swizzle_ps(t2, _MM_SWIZ_REG_BADC));
    t2 = _mm512_add_ps(t2, _mm512_swizzle_ps(t2, _MM_SWIZ_REG_CDAB));

    _mm512_mask_packstorelo_ps(&f, _mm512_mask2int(0x1), t2);
    return f;
}

/******************************************************
 * Single precision half-simd-width utility functions *
 ******************************************************/
static gmx_inline gmx_simd_float_t
gmx_simd_load_dual_hsimd_f_mic(const float * m0,
                               const float * m1)
{
    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    return _mm512_castpd_ps(_mm512_mask_extload_pd(_mm512_extload_pd((const double*)m0, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE),
                                                   mask_hi, (const double*)m1, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE));
}

static gmx_inline gmx_simd_float_t
gmx_simd_loaddup_hsimd_f_mic(const float * m)
{
    assert((size_t)m % 32 == 0);

    return _mm512_castpd_ps(_mm512_extload_pd((const double*)m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE));
}

static gmx_inline gmx_simd_float_t
gmx_simd_load1_dual_hsimd_f_mic(const float * m)
{
    return _mm512_mask_extload_ps(_mm512_extload_ps(m, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE), _mm512_int2mask(0xFF00),
                                  m+1, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
}


static gmx_inline void
gmx_simd_store_dual_hsimd_f_mic(float *           m0,
                                float *           m1,
                                gmx_simd_float_t  a)
{
    __m512 t0;

    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    _mm512_mask_packstorelo_ps(m0, _mm512_int2mask(0x00FF), a);
    _mm512_mask_packstorelo_ps(m1, _mm512_int2mask(0xFF00), a);
}

static gmx_inline void
gmx_simd_decr_hsimd_f_mic(float *          m,
                          gmx_simd_float_t a)
{
    __m512 t;

    assert((size_t)m % 32 == 0);

    t = _mm512_castpd_ps(_mm512_extload_pd((const double*)m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE));
    a = _mm512_add_ps(a, _mm512_permute4f128_ps(a, _MM_PERM_BADC));
    t = _mm512_sub_ps(t, a);
    _mm512_mask_packstorelo_ps(m, _mm512_int2mask(0x00FF), t);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_hsimd_f_mic(const float *       base0,
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

    idx0 = _mm512_loadunpacklo_epi32(_mm512_undefined_epi32(), offset);

    idx0 = _mm512_mullo_epi32(idx0, _mm512_set1_epi32(align));
    idx1 = _mm512_add_epi32(idx0, _mm512_set1_epi32(1));

    idx = _mm512_mask_permute4f128_epi32(idx0, _mm512_int2mask(0xFF00), idx1, PERM_LOW2HIGH);

    tmp1 = _mm512_i32gather_ps(idx, base0, sizeof(float));
    tmp2 = _mm512_i32gather_ps(idx, base1, sizeof(float));

    v0 = _mm512_mask_permute4f128_ps(tmp1, _mm512_int2mask(0xFF00), tmp2, PERM_LOW2HIGH);
    v1 = _mm512_mask_permute4f128_ps(tmp2, _mm512_int2mask(0x00FF), tmp1, PERM_HIGH2LOW);
}


static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_hsimd_f_mic(float *            m,
                                              gmx_simd_float_t   v0,
                                              gmx_simd_float_t   v1)
{
    float  f;
    __m512 t0, t1;

    assert((size_t)m % 32 == 0);

    t0 = _mm512_add_ps(v0, _mm512_swizzle_ps(v0, _MM_SWIZ_REG_BADC));
    t0 = _mm512_mask_add_ps(t0, _mm512_int2mask(0xCCCC), v1, _mm512_swizzle_ps(v1, _MM_SWIZ_REG_BADC));
    t0 =  _mm512_add_ps(t0, _mm512_swizzle_ps(t0, _MM_SWIZ_REG_CDAB));
    t0 = _mm512_add_ps(t0, _mm512_permute4f128_ps(t0, _MM_PERM_CDAB));
    t0 = _mm512_mask_mov_ps(t0, _mm512_int2mask(0b0101010101010101), t0);

    t1 = gmx_simd4_load_f(m);
    t1 = _mm512_add_ps(t1, t0);
    gmx_simd4_store_f(m, t1);

    t0 = _mm512_add_ps(t0, _mm512_swizzle_ps(t0, _MM_SWIZ_REG_BADC));
    t0 = _mm512_add_ps(t0, _mm512_swizzle_ps(t0, _MM_SWIZ_REG_CDAB));

    _mm512_mask_packstorelo_ps(&f, _mm512_mask2int(0x1), t0);
    return f;
}
#endif


/*********************************************************
 * SIMD DOUBLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 *********************************************************/

/* load store double */
static gmx_inline __m512d gmx_simdcall
gmx_simd_loadu_d_mic(const double * m)
{
    return _mm512_loadunpackhi_pd(_mm512_loadunpacklo_pd(_mm512_undefined_pd(), m), m+8);
}

static gmx_inline void gmx_simdcall
gmx_simd_storeu_d_mic(double * m, __m512d s)
{
    _mm512_packstorelo_pd(m, s);
    _mm512_packstorehi_pd(m+8, s);
}

static gmx_inline __m512d gmx_simdcall
gmx_simd_set_exponent_d_mic(__m512d a)
{
    const __m512i expbias      = _mm512_set1_epi32(1023);
    __m512i       iexp         = _mm512_cvtfxpnt_roundpd_epi32lo(a, _MM_FROUND_TO_NEAREST_INT);
    iexp = _mm512_permutevar_epi32(_mm512_set_epi32(7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 0, 0), iexp);
    iexp = _mm512_mask_slli_epi32(_mm512_setzero_epi32(), _mm512_int2mask(0xAAAA), _mm512_add_epi32(iexp, expbias), 20);
    return _mm512_castsi512_pd(iexp);
}

/* load store dint32 */
static gmx_inline void gmx_simdcall
gmx_simd_store_di_mic(gmx_int32_t * m, __m512i s)
{
    assert((size_t)m%32 == 0);
    _mm512_mask_packstorelo_epi32(m, mask_lo, s);
}

static gmx_inline __m512i gmx_simdcall
gmx_simd_loadu_di_mic(const gmx_int32_t * m)
{
    return _mm512_mask_loadunpackhi_epi32(_mm512_mask_loadunpacklo_epi32(_mm512_undefined_epi32(), mask_lo, m), mask_lo, m+16);
}

static gmx_inline void gmx_simdcall
gmx_simd_storeu_di_mic(gmx_int32_t * m, __m512i s)
{
    _mm512_mask_packstorelo_epi32(m, mask_lo, s);
    _mm512_mask_packstorehi_epi32(m+16, mask_lo, s);
}

static gmx_inline gmx_int32_t gmx_simdcall
gmx_simd_extract_di_mic(gmx_simd_dint32_t a, int index)
{
    int r;
    _mm512_mask_packstorelo_epi32(&r, _mm512_mask2int(1<<index), a);
    return r;
}

static gmx_inline void gmx_simdcall
gmx_simd_cvt_f2dd_mic(__m512 f, __m512d * d0, __m512d * d1)
{
    __m512i i1 = _mm512_permute4f128_epi32(_mm512_castps_si512(f), PERM_HIGH2LOW);

    *d0 = _mm512_cvtpslo_pd(f);
    *d1 = _mm512_cvtpslo_pd(_mm512_castsi512_ps(i1));
}

static gmx_inline __m512 gmx_simdcall
gmx_simd_cvt_dd2f_mic(__m512d d0, __m512d d1)
{
    __m512 f0 = _mm512_cvtpd_pslo(d0);
    __m512 f1 = _mm512_cvtpd_pslo(d1);
    return _mm512_mask_permute4f128_ps(f0, mask_hi, f1, PERM_LOW2HIGH);
}

/****************************************************
 * Double precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

/* Declare functions implemented further down that are used here */
template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_mic(const double *       base,
                                               gmx_simd_dint32_t    simdoffset,
                                               gmx_simd_double_t   &v0,
                                               gmx_simd_double_t   &v1,
                                               gmx_simd_double_t   &v2,
                                               gmx_simd_double_t   &v3);

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_mic(const double *       base,
                                               gmx_simd_dint32_t    simdoffset,
                                               gmx_simd_double_t   &v0,
                                               gmx_simd_double_t   &v1);

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d_mic(const double *        base,
                                     const gmx_int32_t     offset[],
                                     gmx_simd_double_t    &v0,
                                     gmx_simd_double_t    &v1,
                                     gmx_simd_double_t    &v2,
                                     gmx_simd_double_t    &v3)
{
    assert((size_t)offset % 32 == 0);
    gmx_simd_gather_load_bysimdint_transpose_d_mic<align>(base,
                                                          gmx_simd_load_di(offset),
                                                          v0, v1, v2, v3);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d_mic(const double *        base,
                                     const gmx_int32_t     offset[],
                                     gmx_simd_double_t    &v0,
                                     gmx_simd_double_t    &v1)
{
    assert((size_t)offset % 32 == 0);
    gmx_simd_gather_load_bysimdint_transpose_d_mic<align>(base,
                                                          gmx_simd_load_di(offset),
                                                          v0, v1);
}

static const int gmx_simd_best_pair_alignment_d_mic = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_d_mic(const double *        base,
                                      const gmx_int32_t     offset[],
                                      gmx_simd_double_t    &v0,
                                      gmx_simd_double_t    &v1,
                                      gmx_simd_double_t    &v2)
{
    __m512i simdoffset;

    assert((size_t)offset % 32 == 0);

    simdoffset = gmx_simd_load_di(offset);

    /* All instructions might be latency ~4 on MIC, so we use shifts where we
     * only need a single instruction (since the shift parameter is an immediate),
     * but multiplication otherwise.
     */
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
gmx_simd_transpose_scatter_storeu_d_mic(double *              base,
                                        const gmx_int32_t     offset[],
                                        gmx_simd_double_t     v0,
                                        gmx_simd_double_t     v1,
                                        gmx_simd_double_t     v2)
{
    __m512i simdoffset;

    assert((size_t)offset % 32 == 0);

    simdoffset = gmx_simd_load_di(offset);

    /* All instructions might be latency ~4 on MIC, so we use shifts where we
     * only need a single instruction (since the shift parameter is an immediate),
     * but multiplication otherwise.
     */
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
gmx_simd_transpose_scatter_incru_d_mic(double *              base,
                                       const gmx_int32_t     offset[],
                                       gmx_simd_double_t     v0,
                                       gmx_simd_double_t     v1,
                                       gmx_simd_double_t     v2)
{
    __m512d t0, t1, t2;
    gmx_simd_gather_loadu_transpose_d_mic<align>(base, offset, t0, t1, t2);
    t0 = _mm512_add_pd(t0, v0);
    t1 = _mm512_add_pd(t1, v1);
    t2 = _mm512_add_pd(t2, v2);
    gmx_simd_transpose_scatter_storeu_d_mic<align>(base, offset, t0, t1, t2);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_d_mic(double *              base,
                                       const gmx_int32_t     offset[],
                                       gmx_simd_double_t     v0,
                                       gmx_simd_double_t     v1,
                                       gmx_simd_double_t     v2)
{
    __m512d t0, t1, t2;
    gmx_simd_gather_loadu_transpose_d_mic<align>(base, offset, t0, t1, t2);
    t0 = _mm512_sub_pd(t0, v0);
    t1 = _mm512_sub_pd(t1, v1);
    t2 = _mm512_sub_pd(t2, v2);
    gmx_simd_transpose_scatter_storeu_d_mic<align>(base, offset, t0, t1, t2);
}


static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_d_mic(gmx_simd_double_t   scalar,
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
gmx_simd_gather_load_bysimdint_transpose_d_mic(const double *       base,
                                               gmx_simd_dint32_t    simdoffset,
                                               gmx_simd_double_t   &v0,
                                               gmx_simd_double_t   &v1,
                                               gmx_simd_double_t   &v2,
                                               gmx_simd_double_t   &v3)
{
    assert((size_t)base % 32 == 0);
    assert(align % 4 == 0);

    /* All instructions might be latency ~4 on MIC, so we use shifts where we
     * only need a single instruction (since the shift parameter is an immediate),
     * but multiplication otherwise.
     */
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
gmx_simd_gather_loadu_bysimdint_transpose_d_mic(const double *       base,
                                                gmx_simd_dint32_t    simdoffset,
                                                gmx_simd_double_t   &v0,
                                                gmx_simd_double_t   &v1)
{
    /* All instructions might be latency ~4 on MIC, so we use shifts where we
     * only need a single instruction (since the shift parameter is an immediate),
     * but multiplication otherwise.
     */
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
gmx_simd_gather_load_bysimdint_transpose_d_mic(const double *       base,
                                               gmx_simd_dint32_t    simdoffset,
                                               gmx_simd_double_t   &v0,
                                               gmx_simd_double_t   &v1)
{
    assert((size_t)base % 32 == 0);
    assert(align % 2 == 0);
    gmx_simd_gather_loadu_bysimdint_transpose_d_mic<align>(base, simdoffset, v0, v1);
}


static gmx_inline double gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_d_mic(double * m,
                                        gmx_simd_double_t v0, gmx_simd_double_t v1,
                                        gmx_simd_double_t v2, gmx_simd_double_t v3)
{
    double  d;
    __m512d t0, t1, t2, t3;

    assert((size_t)m % 32 == 0);

    t0 = _mm512_swizzle_pd(_mm512_mask_blend_pd(_mm512_int2mask(0x33), v0, v2), _MM_SWIZ_REG_BADC);
    t2 = _mm512_mask_blend_pd(_mm512_int2mask(0x33), v2, v0);
    t1 = _mm512_swizzle_pd(_mm512_mask_blend_pd(_mm512_int2mask(0x33), v1, v3), _MM_SWIZ_REG_BADC);
    t3 = _mm512_mask_blend_pd(_mm512_int2mask(0x33), v3, v1);
    t0 = _mm512_add_pd(t0, t2);
    t1 = _mm512_add_pd(t1, t3);

    t2 = _mm512_swizzle_pd(_mm512_mask_blend_pd(_mm512_int2mask(0b01010101), t0, t1), _MM_SWIZ_REG_CDAB);
    t3 = _mm512_mask_blend_pd(_mm512_int2mask(0b01010101), t1, t0);
    t2 = _mm512_add_pd(t2, t3);

    t2 = _mm512_add_pd(t2, _mm512_castps_pd(_mm512_permute4f128_ps(_mm512_castpd_ps(t2), _MM_PERM_BADC)));

    t0 = gmx_simd4_load_d(m);
    t0 = gmx_simd4_add_d(t0, t2);
    gmx_simd4_store_d(m, t0);

    t2 = _mm512_add_pd(t2, _mm512_swizzle_pd(t2, _MM_SWIZ_REG_BADC));
    t2 = _mm512_add_pd(t2, _mm512_swizzle_pd(t2, _MM_SWIZ_REG_CDAB));

    _mm512_mask_packstorelo_pd(&d, _mm512_mask2int(0x01), t2);
    return d;
}

/******************************************************
 * Double precision half-simd-width utility functions *
 ******************************************************/
static gmx_inline gmx_simd_double_t
gmx_simd_load_dual_hsimd_d_mic(const double * m0,
                               const double * m1)
{
    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    return _mm512_mask_extload_pd(_mm512_extload_pd(m0, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE), _mm512_int2mask(0xF0),
                                  m1, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE);
}

static gmx_inline gmx_simd_double_t
gmx_simd_loaddup_hsimd_d_mic(const double * m)
{
    assert((size_t)m % 32 == 0);

    return _mm512_extload_pd(m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE);
}

static gmx_inline gmx_simd_double_t
gmx_simd_load1_dual_hsimd_d_mic(const double * m)
{
    return _mm512_mask_extload_pd(_mm512_extload_pd(m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE), _mm512_int2mask(0xF0),
                                  m+1, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
}


static gmx_inline void
gmx_simd_store_dual_hsimd_d_mic(double *           m0,
                                double *           m1,
                                gmx_simd_double_t  a)
{
    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    _mm512_mask_packstorelo_pd(m0, _mm512_int2mask(0x0F), a);
    _mm512_mask_packstorelo_pd(m1, _mm512_int2mask(0xF0), a);
}

static gmx_inline void
gmx_simd_decr_hsimd_d_mic(double *          m,
                          gmx_simd_double_t a)
{
    __m512d t;

    assert((size_t)m % 32 == 0);

    t = _mm512_extload_pd(m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE);
    a = _mm512_add_pd(a, _mm512_castps_pd(_mm512_permute4f128_ps(_mm512_castpd_ps(a), _MM_PERM_BADC)));
    t = _mm512_sub_pd(t, a);
    _mm512_mask_packstorelo_pd(m, _mm512_int2mask(0x0F), t);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_hsimd_d_mic(const double *       base0,
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

    idx0 = _mm512_extload_epi32(offset, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_4X16, _MM_HINT_NONE);

    idx0 = _mm512_mullo_epi32(idx0, _mm512_set1_epi32(align));
    idx1 = _mm512_add_epi32(idx0, _mm512_set1_epi32(1));

    idx = _mm512_mask_permute4f128_epi32(idx0, _mm512_int2mask(0x00F0), idx1, _MM_PERM_AAAA);

    tmp1 = _mm512_i32logather_pd(idx, base0, sizeof(double));
    tmp2 = _mm512_i32logather_pd(idx, base1, sizeof(double));

    v0 = _mm512_castps_pd(_mm512_mask_permute4f128_ps(_mm512_castpd_ps(tmp1), _mm512_int2mask(0xFF00), _mm512_castpd_ps(tmp2), PERM_LOW2HIGH));
    v1 = _mm512_castps_pd(_mm512_mask_permute4f128_ps(_mm512_castpd_ps(tmp2), _mm512_int2mask(0x00FF), _mm512_castpd_ps(tmp1), PERM_HIGH2LOW));
}

static gmx_inline double
gmx_simd_reduce_incr_4_return_sum_hsimd_d_mic(double *           m,
                                              gmx_simd_double_t  v0,
                                              gmx_simd_double_t  v1)
{
    double   d;
    __m512d  t0, t1;

    assert((size_t)m % 32 == 0);

    t0 = _mm512_add_pd(v0, _mm512_swizzle_pd(v0, _MM_SWIZ_REG_BADC));
    t0 = _mm512_mask_add_pd(t0, _mm512_int2mask(0xCC), v1, _mm512_swizzle_pd(v1, _MM_SWIZ_REG_BADC));
    t0 = _mm512_add_pd(t0, _mm512_swizzle_pd(t0, _MM_SWIZ_REG_CDAB));
    t0 = _mm512_castps_pd(_mm512_mask_permute4f128_ps(_mm512_castpd_ps(t0), _mm512_int2mask(0xCCCC),
                                                      _mm512_castpd_ps(t0), _MM_PERM_DCDC));

    t1 = gmx_simd4_load_d(m);
    t1 = gmx_simd4_add_d(t1, t0);
    gmx_simd4_store_d(m, t1);

    t0 = _mm512_add_pd(t0, _mm512_swizzle_pd(t0, _MM_SWIZ_REG_BADC));
    t0 = _mm512_add_pd(t0, _mm512_swizzle_pd(t0, _MM_SWIZ_REG_CDAB));

    _mm512_mask_packstorelo_pd(&d, _mm512_mask2int(0x03), t0);
    return d;
}
#endif



/*****************************************
 * SIMD4 IMPLEMENTATION HELPER FUNCTIONS *
 *****************************************/

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

static gmx_inline float gmx_simdcall
gmx_simd4_reduce_f_mic(__m512 a)
{
    float f;
    a = _mm512_add_ps(a, _mm512_swizzle_ps(a, _MM_SWIZ_REG_BADC));
    a = _mm512_add_ps(a, _mm512_swizzle_ps(a, _MM_SWIZ_REG_CDAB));
    _mm512_mask_packstorelo_ps(&f, _mm512_mask2int(0x1), a);
    return f;
}

static gmx_inline float gmx_simdcall
gmx_simd4_dotproduct3_f_mic(__m512 a, __m512 b)
{
    a = _mm512_mask_mul_ps(_mm512_setzero_ps(), _mm512_int2mask(0x7), a, b);
    return gmx_simd4_reduce_f(a);
}

#ifdef __cplusplus
static gmx_inline void gmx_simdcall
gmx_simd4_transpose_f_mic(gmx_simd4_float_t &v0, gmx_simd4_float_t &v1,
                          gmx_simd4_float_t &v2, gmx_simd4_float_t &v3)
{
    v0 = _mm512_mask_permute4f128_ps(v0, _mm512_int2mask(0x00F0), v1, _MM_PERM_AAAA);
    v2 = _mm512_mask_permute4f128_ps(v2, _mm512_int2mask(0x00F0), v3, _MM_PERM_AAAA);
    v0 = _mm512_mask_permute4f128_ps(v0, _mm512_int2mask(0xFF00), v2, _MM_PERM_BABA);
    v0 = _mm512_castsi512_ps(_mm512_permutevar_epi32(_mm512_set_epi32(15, 11, 7, 3, 14, 10, 6, 2, 13, 9, 5, 1, 12, 8, 4, 0),
                                                     _mm512_castps_si512(v0)));
    v1 = _mm512_mask_permute4f128_ps(_mm512_setzero_ps(), _mm512_int2mask(0x000F), v0, _MM_PERM_BBBB);
    v2 = _mm512_mask_permute4f128_ps(_mm512_setzero_ps(), _mm512_int2mask(0x000F), v0, _MM_PERM_CCCC);
    v3 = _mm512_mask_permute4f128_ps(_mm512_setzero_ps(), _mm512_int2mask(0x000F), v0, _MM_PERM_DDDD);
}
#endif


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

static gmx_inline double gmx_simdcall
gmx_simd4_reduce_d_mic(__m512d a)
{
    double f;
    a = _mm512_add_pd(a, _mm512_swizzle_pd(a, _MM_SWIZ_REG_BADC));
    a = _mm512_add_pd(a, _mm512_swizzle_pd(a, _MM_SWIZ_REG_CDAB));
    _mm512_mask_packstorelo_pd(&f, _mm512_mask2int(0x1), a);
    return f;
}

static gmx_inline double gmx_simdcall
gmx_simd4_dotproduct3_d_mic(__m512d a, __m512d b)
{
    a = _mm512_mask_mul_pd(_mm512_setzero_pd(), _mm512_int2mask(0x7), a, b);
    return gmx_simd4_reduce_d(a);
}


#ifdef __cplusplus


static gmx_inline void gmx_simdcall
gmx_simd4_transpose_d_mic(gmx_simd4_double_t &v0, gmx_simd4_double_t &v1,
                          gmx_simd4_double_t &v2, gmx_simd4_double_t &v3)
{
    __m512i t0, t1;

    t0 = _mm512_mask_permute4f128_epi32(_mm512_castpd_si512(v0), 0xFF00,
                                        _mm512_castpd_si512(v1), _MM_PERM_BABA);
    t1 = _mm512_mask_permute4f128_epi32(_mm512_castpd_si512(v2), 0xFF00,
                                        _mm512_castpd_si512(v3), _MM_PERM_BABA);

    t0 = _mm512_permutevar_epi32(_mm512_set_epi32(15, 14, 7, 6, 13, 12, 5, 4, 11, 10, 3, 2, 9, 8, 1, 0), t0);
    t1 = _mm512_permutevar_epi32(_mm512_set_epi32(15, 14, 7, 6, 13, 12, 5, 4, 11, 10, 3, 2, 9, 8, 1, 0), t1);

    v0 = _mm512_mask_swizzle_pd(_mm512_castsi512_pd(t0), _mm512_int2mask(0xCC),
                                _mm512_castsi512_pd(t1), _MM_SWIZ_REG_BADC);
    v1 = _mm512_mask_swizzle_pd(_mm512_castsi512_pd(t1), _mm512_int2mask(0x33),
                                _mm512_castsi512_pd(t0), _MM_SWIZ_REG_BADC);

    v2 = _mm512_castps_pd(_mm512_permute4f128_ps(_mm512_castpd_ps(v0), _MM_PERM_DCDC));
    v3 = _mm512_castps_pd(_mm512_permute4f128_ps(_mm512_castpd_ps(v1), _MM_PERM_DCDC));
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

#endif /* GMX_SIMD_IMPL_INTEL_MIC_H */
