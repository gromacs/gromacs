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
#define SimdFloat          __m128
#define simdLoadF           _mm_load_ps
#define simdLoad1F          _mm_load1_ps
#define simdSet1F           _mm_set1_ps
#define simdStoreF          _mm_store_ps
#define simdLoadUF          _mm_loadu_ps
#define simdStoreUF         _mm_storeu_ps
#define simdSetZeroF        _mm_setzero_ps
#define simdAddF            _mm_add_ps
#define simdSubF            _mm_sub_ps
#define simdMulF            _mm_mul_ps
#define simdFmaddF(a, b, c)   _mm_add_ps(_mm_mul_ps(a, b), c)
#define simdFmsubF(a, b, c)   _mm_sub_ps(_mm_mul_ps(a, b), c)
#define simdFnmaddF(a, b, c)  _mm_sub_ps(c, _mm_mul_ps(a, b))
#define simdFnmsubF(a, b, c)  _mm_sub_ps(_mm_setzero_ps(), simdFmaddF(a, b, c))
#define simdAndF            _mm_and_ps
#define simdAndNotF         _mm_andnot_ps
#define simdOrF             _mm_or_ps
#define simdXorF            _mm_xor_ps
#define simdRsqrtF          _mm_rsqrt_ps
#define simdRcpF            _mm_rcp_ps
#define simdAbsF(x)        _mm_andnot_ps(_mm_set1_ps(GMX_FLOAT_NEGZERO), x)
#define simdNegF(x)        _mm_xor_ps(x, _mm_set1_ps(GMX_FLOAT_NEGZERO))
#define simdMaxF            _mm_max_ps
#define simdMinF            _mm_min_ps
#define simdRoundF(x)       _mm_cvtepi32_ps(_mm_cvtps_epi32(x))
#define simdTruncF(x)       _mm_cvtepi32_ps(_mm_cvttps_epi32(x))
#define simdFractionF(x)    _mm_sub_ps(x, simdTruncF(x))
#define simdGetExponentF   simdGetExponentF_sse2
#define simdGetMantissaF   simdGetMantissaF_sse2
#define simdSetExponentF   simdSetExponentF_sse2
/* integer datatype corresponding to float: SimdFInt32 */
#define SimdFInt32         __m128i
#define simdLoadFI(m)       _mm_load_si128((const __m128i *)(m))
#define simdSet1FI          _mm_set1_epi32
#define simdStoreFI(m, x)    _mm_store_si128((__m128i *)(m), x)
#define simdLoadUFI(m)      _mm_loadu_si128((const __m128i *)(m))
#define simdStoreUFI(m, x)   _mm_storeu_si128((__m128i *)(m), x)
#define simdSetZeroFI       _mm_setzero_si128
#define simdCvtF2I          _mm_cvtps_epi32
#define simdCvttF2I         _mm_cvttps_epi32
#define simdCvtI2F          _mm_cvtepi32_ps
#define simdExtractFI(x, i)  _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (i)))
/* Integer logical ops on SimdFInt32 */
#define simdSlliFI          _mm_slli_epi32
#define simdSrliFI          _mm_srli_epi32
#define simdAndFI           _mm_and_si128
#define simdAndNotFI        _mm_andnot_si128
#define simdOrFI            _mm_or_si128
#define simdXorFI           _mm_xor_si128
/* Integer arithmetic ops on SimdFInt32 */
#define simdAddFI           _mm_add_epi32
#define simdSubFI           _mm_sub_epi32
#define simdMulFI           simdMulFI_sse2
/* Boolean & comparison operations on SimdFloat */
#define SimdFBool          __m128
#define simdCmpEqF          _mm_cmpeq_ps
#define simdCmpLtF          _mm_cmplt_ps
#define simdCmpLeF          _mm_cmple_ps
#define simdAndFB           _mm_and_ps
#define simdOrFB            _mm_or_ps
#define simdAnyTrueFB       _mm_movemask_ps
#define simdMaskF      _mm_and_ps
#define simdMaskNotF(a, sel)   _mm_andnot_ps(sel, a)
#define simdBlendF(a, b, s)  _mm_or_ps(_mm_andnot_ps(s, a), _mm_and_ps(s, b))
#define simdReduceF(a)      simdReduceF_sse2(a)
/* Boolean & comparison operations on SimdFInt32 */
#define SimdFIBool         __m128i
#define simdCmpEqFI         _mm_cmpeq_epi32
#define simdCmpLtFI         _mm_cmplt_epi32
#define simdAndFIB          _mm_and_si128
#define simdOrFIB           _mm_or_si128
#define simdAnyTrueFIB      _mm_movemask_epi8
#define simdMaskFI     _mm_and_si128
#define simdMaskNotFI(a, sel) _mm_andnot_si128(sel, a)
#define simdBlendFI(a, b, s) _mm_or_si128(_mm_andnot_si128(s, a), _mm_and_si128(s, b))
/* Conversions between different booleans */
#define simdCvtFB2FIB       _mm_castps_si128
#define simdCvtFIB2FB       _mm_castsi128_ps

/****************************************************
 * SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static inline __m128 gmx_simdcall
simdGetExponentF_sse2(__m128 x)
{
    const __m128  expmask      = _mm_castsi128_ps(_mm_set1_epi32(0x7F800000));
    const __m128i expbias      = _mm_set1_epi32(127);
    __m128i       iexp;

    iexp = _mm_castps_si128(_mm_and_ps(x, expmask));
    iexp = _mm_sub_epi32(_mm_srli_epi32(iexp, 23), expbias);
    return _mm_cvtepi32_ps(iexp);
}

static inline __m128 gmx_simdcall
simdGetMantissaF_sse2(__m128 x)
{
    const __m128 mantmask = _mm_castsi128_ps(_mm_set1_epi32(0x007FFFFF));
    const __m128 one      = _mm_set1_ps(1.0f);

    x = _mm_and_ps(x, mantmask);
    return _mm_or_ps(x, one);
}

static inline __m128 gmx_simdcall
simdSetExponentF_sse2(__m128 x)
{
    const __m128i expbias      = _mm_set1_epi32(127);
    __m128i       iexp         = _mm_cvtps_epi32(x);

    iexp = _mm_slli_epi32(_mm_add_epi32(iexp, expbias), 23);
    return _mm_castsi128_ps(iexp);
}

static inline __m128i gmx_simdcall
simdMulFI_sse2(__m128i a, __m128i b)
{
    __m128i a1 = _mm_srli_si128(a, 4); /* - a[3] a[2] a[1] */
    __m128i b1 = _mm_srli_si128(b, 4); /* - b[3] b[2] b[1] */
    __m128i c  = _mm_mul_epu32(a, b);
    __m128i c1 = _mm_mul_epu32(a1, b1);

    c  = _mm_shuffle_epi32(c, _MM_SHUFFLE(3, 1, 2, 0));  /* - - a[2]*b[2] a[0]*b[0] */
    c1 = _mm_shuffle_epi32(c1, _MM_SHUFFLE(3, 1, 2, 0)); /* - - a[3]*b[3] a[1]*b[1] */

    return _mm_unpacklo_epi32(c, c1);
}

static inline float gmx_simdcall
simdReduceF_sse2(__m128 a)
{
    __m128 b;
    float  f;
    b = _mm_add_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(1, 0, 3, 2)));
    b = _mm_add_ss(b, _mm_shuffle_ps(b, b, _MM_SHUFFLE(0, 3, 2, 1)));
    _mm_store_ss(&f, b);
    return f;
}

#endif /* GMX_SIMD_IMPL_X86_SSE2_SIMD_FLOAT_H */
