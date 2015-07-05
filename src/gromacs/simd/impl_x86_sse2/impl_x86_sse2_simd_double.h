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

#ifndef GMX_SIMD_IMPL_X86_SSE2_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_X86_SSE2_SIMD_DOUBLE_H

#include "config.h"

#include <emmintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_sse2_common.h"

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_double_t          __m128d
#define gmx_simd_load_d            _mm_load_pd
#define gmx_simd_load1_d           _mm_load1_pd
#define gmx_simd_set1_d            _mm_set1_pd
#define gmx_simd_store_d           _mm_store_pd
#define gmx_simd_loadu_d           _mm_loadu_pd
#define gmx_simd_storeu_d          _mm_storeu_pd
#define gmx_simd_setzero_d         _mm_setzero_pd
#define gmx_simd_add_d             _mm_add_pd
#define gmx_simd_sub_d             _mm_sub_pd
#define gmx_simd_mul_d             _mm_mul_pd
#define gmx_simd_fmadd_d(a, b, c)    _mm_add_pd(_mm_mul_pd(a, b), c)
#define gmx_simd_fmsub_d(a, b, c)    _mm_sub_pd(_mm_mul_pd(a, b), c)
#define gmx_simd_fnmadd_d(a, b, c)   _mm_sub_pd(c, _mm_mul_pd(a, b))
#define gmx_simd_fnmsub_d(a, b, c)   _mm_sub_pd(_mm_setzero_pd(), gmx_simd_fmadd_d(a, b, c))
#define gmx_simd_and_d             _mm_and_pd
#define gmx_simd_andnot_d          _mm_andnot_pd
#define gmx_simd_or_d              _mm_or_pd
#define gmx_simd_xor_d             _mm_xor_pd
#define gmx_simd_rsqrt_d(x)        _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(x)))
/* Don't use FMA for sqrt N-R iterations - this saves 1 instruction without FMA hardware */
#define gmx_simd_rcp_d(x)          _mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(x)))
#define gmx_simd_fabs_d(x)         _mm_andnot_pd(_mm_set1_pd(GMX_DOUBLE_NEGZERO), x)
#define gmx_simd_fneg_d(x)         _mm_xor_pd(x, _mm_set1_pd(GMX_DOUBLE_NEGZERO))
#define gmx_simd_max_d             _mm_max_pd
#define gmx_simd_min_d             _mm_min_pd
#define gmx_simd_round_d(x)        _mm_cvtepi32_pd(_mm_cvtpd_epi32(x))
#define gmx_simd_trunc_d(x)        _mm_cvtepi32_pd(_mm_cvttpd_epi32(x))
#define gmx_simd_fraction_d(x)     _mm_sub_pd(x, gmx_simd_trunc_d(x))
#define gmx_simd_get_exponent_d    gmx_simd_get_exponent_d_sse2
#define gmx_simd_get_mantissa_d    gmx_simd_get_mantissa_d_sse2
#define gmx_simd_set_exponent_d    gmx_simd_set_exponent_d_sse2
/* integer datatype corresponding to double: gmx_simd_dint32_t */
#define gmx_simd_dint32_t          __m128i
#define gmx_simd_load_di(m)        _mm_loadl_epi64((const __m128i *)(m))
#define gmx_simd_set1_di           _mm_set1_epi32
#define gmx_simd_store_di(m, x)     _mm_storel_epi64((__m128i *)(m), x)
#define gmx_simd_loadu_di(m)       _mm_loadl_epi64((const __m128i *)(m))
#define gmx_simd_storeu_di(m, x)    _mm_storel_epi64((__m128i *)(m), x)
#define gmx_simd_setzero_di        _mm_setzero_si128
#define gmx_simd_cvt_d2i           _mm_cvtpd_epi32
#define gmx_simd_cvtt_d2i          _mm_cvttpd_epi32
#define gmx_simd_cvt_i2d           _mm_cvtepi32_pd
#define gmx_simd_extract_di(x, i)   _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (i)))
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
#define gmx_simd_mul_di            gmx_simd_mul_di_sse2
/* Boolean & comparison operations on gmx_simd_double_t */
#define gmx_simd_dbool_t            __m128d
#define gmx_simd_cmpeq_d            _mm_cmpeq_pd
#define gmx_simd_cmplt_d            _mm_cmplt_pd
#define gmx_simd_cmple_d            _mm_cmple_pd
#define gmx_simd_and_db             _mm_and_pd
#define gmx_simd_or_db              _mm_or_pd
#define gmx_simd_anytrue_db         _mm_movemask_pd
#define gmx_simd_blendzero_d        _mm_and_pd
#define gmx_simd_blendnotzero_d(a, sel) _mm_andnot_pd(sel, a)
#define gmx_simd_blendv_d(a, b, sel)  _mm_or_pd(_mm_andnot_pd(sel, a), _mm_and_pd(sel, b))
#define gmx_simd_reduce_d(a)        gmx_simd_reduce_d_sse2(a)

/* Boolean & comparison operations on gmx_simd_dint32_t */
#define gmx_simd_dibool_t           __m128i
#define gmx_simd_cmpeq_di           _mm_cmpeq_epi32
#define gmx_simd_cmplt_di           _mm_cmplt_epi32
#define gmx_simd_and_dib            _mm_and_si128
#define gmx_simd_or_dib             _mm_or_si128
#define gmx_simd_anytrue_dib(x)     _mm_movemask_epi8(_mm_shuffle_epi32(x, _MM_SHUFFLE(1, 0, 1, 0)))
#define gmx_simd_blendzero_di       _mm_and_si128
#define gmx_simd_blendnotzero_di(a, sel)  _mm_andnot_si128(sel, a)
#define gmx_simd_blendv_di(a, b, sel) _mm_or_si128(_mm_andnot_si128(sel, a), _mm_and_si128(sel, b))
#define gmx_simd_cvt_db2dib(x)      _mm_shuffle_epi32(_mm_castpd_si128(x), _MM_SHUFFLE(2, 0, 2, 0))
#define gmx_simd_cvt_dib2db(x)      _mm_castsi128_pd(_mm_shuffle_epi32(x, _MM_SHUFFLE(1, 1, 0, 0)))
/* Float/double conversion */
#define gmx_simd_cvt_f2dd(f, d0, d1)  { *d0 = _mm_cvtps_pd(f); *d1 = _mm_cvtps_pd(_mm_movehl_ps(f, f)); }
#define gmx_simd_cvt_dd2f(d0, d1)    _mm_movelh_ps(_mm_cvtpd_ps(d0), _mm_cvtpd_ps(d1))


/****************************************************
 * DOUBLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static gmx_inline __m128d gmx_simdcall
gmx_simd_get_exponent_d_sse2(__m128d x)
{
    /* Don't use _mm_set1_epi64x() - on MSVC it is only supported for 64-bit builds */
    const __m128d expmask      = _mm_castsi128_pd( _mm_set_epi32(0x7FF00000, 0x00000000, 0x7FF00000, 0x00000000) );
    const __m128i expbias      = _mm_set1_epi32(1023);
    __m128i       iexp;

    iexp   = _mm_castpd_si128(_mm_and_pd(x, expmask));
    iexp   = _mm_sub_epi32(_mm_srli_epi64(iexp, 52), expbias);
    iexp   = _mm_shuffle_epi32(iexp, _MM_SHUFFLE(3, 1, 2, 0) );
    return _mm_cvtepi32_pd(iexp);
}

static gmx_inline __m128d gmx_simdcall
gmx_simd_get_mantissa_d_sse2(__m128d x)
{
    /* Don't use _mm_set1_epi64x() - on MSVC it is only supported for 64-bit builds */
    const __m128d mantmask = _mm_castsi128_pd( _mm_set_epi32(0x000FFFFF, 0xFFFFFFFF, 0x000FFFFF, 0xFFFFFFFF) );
    const __m128d one      = _mm_set1_pd(1.0);

    x = _mm_and_pd(x, mantmask);
    return _mm_or_pd(x, one);
}

static gmx_inline __m128d gmx_simdcall
gmx_simd_set_exponent_d_sse2(__m128d x)
{
    const __m128i  expbias      = _mm_set1_epi32(1023);
    __m128i        iexp         = _mm_cvtpd_epi32(x);

    /* After conversion integers will be in slot 0,1. Move them to 0,2 so
     * we can do a 64-bit shift and get them to the dp exponents. */
    iexp = _mm_shuffle_epi32(iexp, _MM_SHUFFLE(3, 1, 2, 0));
    iexp = _mm_slli_epi64(_mm_add_epi32(iexp, expbias), 52);
    return _mm_castsi128_pd(iexp);
}

static gmx_inline __m128i gmx_simdcall
gmx_simd_mul_di_sse2(__m128i a, __m128i b)
{
    __m128i c;

    a = _mm_unpacklo_epi32(a, _mm_setzero_si128());       /* 0 a[1] 0 a[0] */
    b = _mm_unpacklo_epi32(b, _mm_setzero_si128());       /* 0 b[1] 0 b[0] */

    c  = _mm_mul_epu32(a, b);                             /* 0 a[1]*b[1] 0 a[0]*b[0] */
    return _mm_shuffle_epi32(c, _MM_SHUFFLE(3, 1, 2, 0)); /* 0 0 a[1]*b[1] a[0]*b[0] */
}

static gmx_inline double gmx_simdcall
gmx_simd_reduce_d_sse2(__m128d a)
{
    __m128d b;
    double  f;

    b = _mm_add_sd(a, _mm_shuffle_pd(a, a, _MM_SHUFFLE2(1, 1)));
    _mm_store_sd(&f, b);
    return f;
}

#endif /* GMX_SIMD_IMPL_X86_SSE2_SIMD_DOUBLE_H */
