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

#ifndef GMX_SIMD_IMPL_X86_AVX_512F_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_512F_SIMD_DOUBLE_H

#include <math.h>

#include <immintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_avx_512f_common.h"

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
#define gmx_simd_fneg_d(x)         gmx_simd_xor_d(x, _mm512_set1_pd(GMX_DOUBLE_NEGZERO))
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


/* Implementation helpers */

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

#endif /* GMX_SIMD_IMPL_X86_AVX_512F_SIMD_DOUBLE_H */
