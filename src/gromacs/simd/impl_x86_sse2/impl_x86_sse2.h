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

#ifndef GMX_SIMD_IMPL_X86_SSE2_H
#define GMX_SIMD_IMPL_X86_SSE2_H

#include "config.h"

#include <math.h>

#include <emmintrin.h>

/* Set capabilities that can be inherited */
#define GMX_SIMD_X86_SSE2_OR_HIGHER

/* x86 SSE2 SIMD instruction wrappers
 *
 * Please see documentation in gromacs/simd/simd.h for defines.
 */

/* Capability definitions for SSE2 */
#define GMX_SIMD_HAVE_FLOAT
#define GMX_SIMD_HAVE_DOUBLE
#define GMX_SIMD_HAVE_HARDWARE
#define GMX_SIMD_HAVE_LOADU
#define GMX_SIMD_HAVE_STOREU
#define GMX_SIMD_HAVE_LOGICAL
#undef  GMX_SIMD_HAVE_FMA
#undef  GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD_HAVE_FINT32
#define GMX_SIMD_HAVE_FINT32_EXTRACT   /* No SSE2 instruction, but use shifts */
#define GMX_SIMD_HAVE_FINT32_LOGICAL
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS
#define GMX_SIMD_HAVE_DINT32
#define GMX_SIMD_HAVE_DINT32_EXTRACT   /* No SSE2 instruction, but use shifts */
#define GMX_SIMD_HAVE_DINT32_LOGICAL
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS
#define GMX_SIMD4_HAVE_FLOAT
#undef  GMX_SIMD4_HAVE_DOUBLE

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH         4
#define GMX_SIMD_DOUBLE_WIDTH        2
#define GMX_SIMD_FINT32_WIDTH        4
#define GMX_SIMD_DINT32_WIDTH        2
#define GMX_SIMD_RSQRT_BITS         11
#define GMX_SIMD_RCP_BITS           11

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_float_t          __m128
#define gmx_simd_load_f           _mm_load_ps
#define gmx_simd_load1_f          _mm_load1_ps
#define gmx_simd_set1_f           _mm_set1_ps
#define gmx_simd_store_f          _mm_store_ps
#define gmx_simd_loadu_f          _mm_loadu_ps
#define gmx_simd_storeu_f         _mm_storeu_ps
#define gmx_simd_setzero_f        _mm_setzero_ps
#define gmx_simd_add_f            _mm_add_ps
#define gmx_simd_sub_f            _mm_sub_ps
#define gmx_simd_mul_f            _mm_mul_ps
#define gmx_simd_fmadd_f(a, b, c)   _mm_add_ps(_mm_mul_ps(a, b), c)
#define gmx_simd_fmsub_f(a, b, c)   _mm_sub_ps(_mm_mul_ps(a, b), c)
#define gmx_simd_fnmadd_f(a, b, c)  _mm_sub_ps(c, _mm_mul_ps(a, b))
#define gmx_simd_fnmsub_f(a, b, c)  _mm_sub_ps(_mm_setzero_ps(), gmx_simd_fmadd_f(a, b, c))
#define gmx_simd_and_f            _mm_and_ps
#define gmx_simd_andnot_f         _mm_andnot_ps
#define gmx_simd_or_f             _mm_or_ps
#define gmx_simd_xor_f            _mm_xor_ps
#define gmx_simd_rsqrt_f          _mm_rsqrt_ps
#define gmx_simd_rcp_f            _mm_rcp_ps
#define gmx_simd_fabs_f(x)        _mm_andnot_ps(_mm_set1_ps(GMX_FLOAT_NEGZERO), x)
#define gmx_simd_fneg_f(x)        _mm_xor_ps(x, _mm_set1_ps(GMX_FLOAT_NEGZERO))
#define gmx_simd_max_f            _mm_max_ps
#define gmx_simd_min_f            _mm_min_ps
#define gmx_simd_round_f(x)       _mm_cvtepi32_ps(_mm_cvtps_epi32(x))
#define gmx_simd_trunc_f(x)       _mm_cvtepi32_ps(_mm_cvttps_epi32(x))
#define gmx_simd_fraction_f(x)    _mm_sub_ps(x, gmx_simd_trunc_f(x))
#define gmx_simd_get_exponent_f   gmx_simd_get_exponent_f_sse2
#define gmx_simd_get_mantissa_f   gmx_simd_get_mantissa_f_sse2
#define gmx_simd_set_exponent_f   gmx_simd_set_exponent_f_sse2
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_fint32_t         __m128i
#define gmx_simd_load_fi(m)       _mm_load_si128((const __m128i *)(m))
#define gmx_simd_set1_fi          _mm_set1_epi32
#define gmx_simd_store_fi(m, x)    _mm_store_si128((__m128i *)(m), x)
#define gmx_simd_loadu_fi(m)      _mm_loadu_si128((const __m128i *)(m))
#define gmx_simd_storeu_fi(m, x)   _mm_storeu_si128((__m128i *)(m), x)
#define gmx_simd_setzero_fi       _mm_setzero_si128
#define gmx_simd_cvt_f2i          _mm_cvtps_epi32
#define gmx_simd_cvtt_f2i         _mm_cvttps_epi32
#define gmx_simd_cvt_i2f          _mm_cvtepi32_ps
#define gmx_simd_extract_fi(x, i)  _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (i)))
/* Integer logical ops on gmx_simd_fint32_t */
#define gmx_simd_slli_fi          _mm_slli_epi32
#define gmx_simd_srli_fi          _mm_srli_epi32
#define gmx_simd_and_fi           _mm_and_si128
#define gmx_simd_andnot_fi        _mm_andnot_si128
#define gmx_simd_or_fi            _mm_or_si128
#define gmx_simd_xor_fi           _mm_xor_si128
/* Integer arithmetic ops on gmx_simd_fint32_t */
#define gmx_simd_add_fi           _mm_add_epi32
#define gmx_simd_sub_fi           _mm_sub_epi32
#define gmx_simd_mul_fi           gmx_simd_mul_fi_sse2
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_fbool_t          __m128
#define gmx_simd_cmpeq_f          _mm_cmpeq_ps
#define gmx_simd_cmplt_f          _mm_cmplt_ps
#define gmx_simd_cmple_f          _mm_cmple_ps
#define gmx_simd_and_fb           _mm_and_ps
#define gmx_simd_or_fb            _mm_or_ps
#define gmx_simd_anytrue_fb       _mm_movemask_ps
#define gmx_simd_blendzero_f      _mm_and_ps
#define gmx_simd_blendnotzero_f(a, sel)   _mm_andnot_ps(sel, a)
#define gmx_simd_blendv_f(a, b, s)  _mm_or_ps(_mm_andnot_ps(s, a), _mm_and_ps(s, b))
#define gmx_simd_reduce_f(a)      gmx_simd_reduce_f_sse2(a)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_fibool_t         __m128i
#define gmx_simd_cmpeq_fi         _mm_cmpeq_epi32
#define gmx_simd_cmplt_fi         _mm_cmplt_epi32
#define gmx_simd_and_fib          _mm_and_si128
#define gmx_simd_or_fib           _mm_or_si128
#define gmx_simd_anytrue_fib      _mm_movemask_epi8
#define gmx_simd_blendzero_fi     _mm_and_si128
#define gmx_simd_blendnotzero_fi(a, sel) _mm_andnot_si128(sel, a)
#define gmx_simd_blendv_fi(a, b, s) _mm_or_si128(_mm_andnot_si128(s, a), _mm_and_si128(s, b))
/* Conversions between different booleans */
#define gmx_simd_cvt_fb2fib       _mm_castps_si128
#define gmx_simd_cvt_fib2fb       _mm_castsi128_ps

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
 * SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static gmx_inline __m128 gmx_simdcall
gmx_simd_get_exponent_f_sse2(__m128 x)
{
    const __m128  expmask      = _mm_castsi128_ps(_mm_set1_epi32(0x7F800000));
    const __m128i expbias      = _mm_set1_epi32(127);
    __m128i       iexp;

    iexp = _mm_castps_si128(_mm_and_ps(x, expmask));
    iexp = _mm_sub_epi32(_mm_srli_epi32(iexp, 23), expbias);
    return _mm_cvtepi32_ps(iexp);
}

static gmx_inline __m128 gmx_simdcall
gmx_simd_get_mantissa_f_sse2(__m128 x)
{
    const __m128 mantmask = _mm_castsi128_ps(_mm_set1_epi32(0x007FFFFF));
    const __m128 one      = _mm_set1_ps(1.0f);

    x = _mm_and_ps(x, mantmask);
    return _mm_or_ps(x, one);
}

static gmx_inline __m128 gmx_simdcall
gmx_simd_set_exponent_f_sse2(__m128 x)
{
    const __m128i expbias      = _mm_set1_epi32(127);
    __m128i       iexp         = _mm_cvtps_epi32(x);

    iexp = _mm_slli_epi32(_mm_add_epi32(iexp, expbias), 23);
    return _mm_castsi128_ps(iexp);
}

static gmx_inline __m128i gmx_simdcall
gmx_simd_mul_fi_sse2(__m128i a, __m128i b)
{
    __m128i a1 = _mm_srli_si128(a, 4); /* - a[3] a[2] a[1] */
    __m128i b1 = _mm_srli_si128(b, 4); /* - b[3] b[2] b[1] */
    __m128i c  = _mm_mul_epu32(a, b);
    __m128i c1 = _mm_mul_epu32(a1, b1);

    c  = _mm_shuffle_epi32(c, _MM_SHUFFLE(3, 1, 2, 0));  /* - - a[2]*b[2] a[0]*b[0] */
    c1 = _mm_shuffle_epi32(c1, _MM_SHUFFLE(3, 1, 2, 0)); /* - - a[3]*b[3] a[1]*b[1] */

    return _mm_unpacklo_epi32(c, c1);
}

static gmx_inline float gmx_simdcall
gmx_simd_reduce_f_sse2(__m128 a)
{
    __m128 b;
    float  f;
    b = _mm_add_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(1, 0, 3, 2)));
    b = _mm_add_ss(b, _mm_shuffle_ps(b, b, _MM_SHUFFLE(0, 3, 2, 1)));
    _mm_store_ss(&f, b);
    return f;
}

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

/* SSE2 is already 4-wide in single, so we just reuse float datatype for SIMD4.
 * SSE2 cannot do double-precision SIMD4.
 */
#define gmx_simd4_float_t                gmx_simd_float_t
#define gmx_simd4_load_f                 gmx_simd_load_f
#define gmx_simd4_load1_f                gmx_simd_load1_f
#define gmx_simd4_set1_f                 gmx_simd_set1_f
#define gmx_simd4_store_f                gmx_simd_store_f
#define gmx_simd4_loadu_f                gmx_simd_loadu_f
#define gmx_simd4_storeu_f               gmx_simd_storeu_f
#define gmx_simd4_setzero_f              gmx_simd_setzero_f
#define gmx_simd4_add_f                  gmx_simd_add_f
#define gmx_simd4_sub_f                  gmx_simd_sub_f
#define gmx_simd4_mul_f                  gmx_simd_mul_f
#define gmx_simd4_fmadd_f                gmx_simd_fmadd_f
#define gmx_simd4_fmsub_f                gmx_simd_fmsub_f
#define gmx_simd4_fnmadd_f               gmx_simd_fnmadd_f
#define gmx_simd4_fnmsub_f               gmx_simd_fnmsub_f
#define gmx_simd4_and_f                  gmx_simd_and_f
#define gmx_simd4_andnot_f               gmx_simd_andnot_f
#define gmx_simd4_or_f                   gmx_simd_or_f
#define gmx_simd4_xor_f                  gmx_simd_xor_f
#define gmx_simd4_rsqrt_f                gmx_simd_rsqrt_f
#define gmx_simd4_fabs_f                 gmx_simd_fabs_f
#define gmx_simd4_fneg_f                 gmx_simd_fneg_f
#define gmx_simd4_max_f                  gmx_simd_max_f
#define gmx_simd4_min_f                  gmx_simd_min_f
#define gmx_simd4_round_f                gmx_simd_round_f
#define gmx_simd4_trunc_f                gmx_simd_trunc_f
#define gmx_simd4_dotproduct3_f          gmx_simd4_dotproduct3_f_sse2
#define gmx_simd4_fbool_t                gmx_simd_fbool_t
#define gmx_simd4_cmpeq_f                gmx_simd_cmpeq_f
#define gmx_simd4_cmplt_f                gmx_simd_cmplt_f
#define gmx_simd4_cmple_f                gmx_simd_cmple_f
#define gmx_simd4_and_fb                 gmx_simd_and_fb
#define gmx_simd4_or_fb                  gmx_simd_or_fb
#define gmx_simd4_anytrue_fb             gmx_simd_anytrue_fb
#define gmx_simd4_blendzero_f            gmx_simd_blendzero_f
#define gmx_simd4_blendnotzero_f         gmx_simd_blendnotzero_f
#define gmx_simd4_blendv_f               gmx_simd_blendv_f
#define gmx_simd4_reduce_f               gmx_simd_reduce_f

/* SIMD4 Dotproduct helper function */
static gmx_inline float gmx_simdcall
gmx_simd4_dotproduct3_f_sse2(__m128 a, __m128 b)
{
    float  f;
    __m128 c;
    a = _mm_mul_ps(a, b);
    c = _mm_add_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(2, 1, 2, 1)));
    c = _mm_add_ps(c, _mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 2, 3, 2)));
    _mm_store_ss(&f, c);
    return f;
}

#endif /* GMX_SIMD_IMPL_X86_SSE2_H */
