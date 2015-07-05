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

#ifndef GMX_SIMD_IMPL_X86_SSE2_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_SSE2_SIMD_FLOAT_H

#include "config.h"

#include <emmintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_sse2_common.h"

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

#endif /* GMX_SIMD_IMPL_X86_SSE2_SIMD_FLOAT_H */
