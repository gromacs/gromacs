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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_256_SIMD_DOUBLE_H

#include "config.h"

#include <immintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_avx_256_common.h"

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
#define gmx_simd_fabs_d(x)         _mm256_andnot_pd(_mm256_set1_pd(GMX_DOUBLE_NEGZERO), x)
#define gmx_simd_fneg_d(x)         _mm256_xor_pd(x, _mm256_set1_pd(GMX_DOUBLE_NEGZERO))
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

#endif /* GMX_SIMD_IMPL_X86_AVX_256_SIMD_DOUBLE_H */
