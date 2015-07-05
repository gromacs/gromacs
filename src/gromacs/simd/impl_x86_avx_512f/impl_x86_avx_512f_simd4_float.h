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

#ifndef GMX_SIMD_IMPL_X86_AVX_512F_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_512F_SIMD4_FLOAT_H

#include <math.h>

#include <immintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_avx_512f_common.h"
#include "impl_x86_avx_512f_simd_float.h"

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
#define gmx_simd4_fneg_f(x)         _mm_xor_ps(x, _mm_set1_ps(GMX_FLOAT_NEGZERO))
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


/* Implementation helpers */
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

#endif /* GMX_SIMD_IMPL_X86_AVX_512F_SIMD4_FLOAT_H */
