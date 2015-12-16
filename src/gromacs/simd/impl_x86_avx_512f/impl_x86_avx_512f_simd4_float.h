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
#define Simd4Float           __m128
#define simd4LoadF            _mm_load_ps
#define simd4Load1F           _mm_load1_ps
#define simd4Set1F            _mm_set1_ps
#define simd4StoreF           _mm_store_ps
#define simd4LoadUF           _mm_loadu_ps
#define simd4StoreUF          _mm_storeu_ps
#define simd4SetZeroF         _mm_setzero_ps
#define simd4AddF             _mm_add_ps
#define simd4SubF             _mm_sub_ps
#define simd4MulF             _mm_mul_ps
#define simd4FmaddF           _mm_fmadd_ps
#define simd4FmsubF           _mm_fmsub_ps
#define simd4FnmaddF          _mm_fnmadd_ps
#define simd4FnmsubF          _mm_fnmsub_ps
#define simd4AndF             _mm_and_ps
#define simd4AndNotF          _mm_andnot_ps
#define simd4OrF              _mm_or_ps
#define simd4XorF             _mm_xor_ps
/* We need to use the new table lookup instructions since we have specified
 * 14 bits of accuracy for AVX-512F.
 */
#define simd4RsqrtF(x)        _mm512_castps512_ps128(_mm512_rsqrt14_ps(_mm512_castps128_ps512(x)))
/* abs/neg cannot cause FP exceptions, so we can operate on entire register */
#define simd4AbsF(x)         _mm512_castps512_ps128(_mm512_abs_ps(_mm512_castps128_ps512(x)))
#define simd4NegF(x)         _mm_xor_ps(x, _mm_set1_ps(GMX_FLOAT_NEGZERO))
#define simd4MaxF             _mm_max_ps
#define simd4MinF             _mm_min_ps
#define simd4RoundF(x)        _mm_round_ps(x, _MM_FROUND_NINT)
#define simd4TruncF(x)        _mm_round_ps(x, _MM_FROUND_TRUNC)
#define simd4DotProductF(a, b) simd4DotProductF_x86_avx_512f(a, b)
#define Simd4FBool           __mmask16
#define simd4CmpEqF(a, b)     _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), _mm512_castps128_ps512(a), _mm512_castps128_ps512(b), _CMP_EQ_OQ)
#define simd4CmpLtF(a, b)     _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), _mm512_castps128_ps512(a), _mm512_castps128_ps512(b), _CMP_LT_OS)
#define simd4CmpLeF(a, b)     _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), _mm512_castps128_ps512(a), _mm512_castps128_ps512(b), _CMP_LE_OS)
#define simd4AndFB            _mm512_kand
#define simd4OrFB             _mm512_kor
#define simd4AnyTrueFB(x)     (_mm512_mask2int(x)&0xF)
#define simd4MaskF(a, sel)    _mm512_castps512_ps128(_mm512_mask_mov_ps(_mm512_setzero_ps(), sel, _mm512_castps128_ps512(a)))
#define simd4MaskNotF(a, sel) _mm512_castps512_ps128(_mm512_mask_mov_ps(_mm512_setzero_ps(), _mm512_knot(sel), _mm512_castps128_ps512(a)))
#define simd4BlendF(a, b, sel)    _mm512_castps512_ps128(_mm512_mask_blend_ps(sel, _mm512_castps128_ps512(a), _mm512_castps128_ps512(b)))
#define simd4ReduceF(x)       simd4ReduceF_x86_avx_512f(x)


/* Implementation helpers */
static inline float
simd4ReduceF_x86_avx_512f(__m128 a)
{
    float f;
    a = _mm_hadd_ps(a, a);
    a = _mm_hadd_ps(a, a);
    _mm_store_ss(&f, a);
    return f;
}

static inline float
simd4DotProductF_x86_avx_512f(__m128 a, __m128 b)
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
