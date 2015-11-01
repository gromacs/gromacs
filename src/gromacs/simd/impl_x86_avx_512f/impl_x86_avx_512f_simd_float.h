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

#ifndef GMX_SIMD_IMPL_X86_AVX_512F_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_512F_SIMD_FLOAT_H

#include <math.h>

#include <immintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_avx_512f_common.h"

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdFloat           __m512
#define simdLoadF            _mm512_load_ps
/* Avoid using _mm512_extload_ps() since it is not available on gcc-4.9 */
#define simdLoad1F(m)        _mm512_broadcastss_ps(_mm_broadcast_ss(m))
#define simdSet1F            _mm512_set1_ps
#define simdStoreF           _mm512_store_ps
#define simdLoadUF           _mm512_loadu_ps
#define simdStoreUF          _mm512_storeu_ps
#define simdSetZeroF         _mm512_setzero_ps
#define simdAddF             _mm512_add_ps
#define simdSubF             _mm512_sub_ps
#define simdMulF             _mm512_mul_ps
#define simdFmaddF           _mm512_fmadd_ps
#define simdFmsubF           _mm512_fmsub_ps
#define simdFnmaddF          _mm512_fnmadd_ps
#define simdFnmsubF          _mm512_fnmsub_ps
#define simdAndF(a, b)        _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define simdAndNotF(a, b)     _mm512_castsi512_ps(_mm512_andnot_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define simdOrF(a, b)         _mm512_castsi512_ps(_mm512_or_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define simdXorF(a, b)        _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(a), _mm512_castps_si512(b)))
#define simdRsqrtF           _mm512_rsqrt14_ps
#define simdRcpF             _mm512_rcp14_ps
#define simdAbsF(x)         _mm512_abs_ps(x)
#define simdNegF(x)         simdXorF(x, _mm512_set1_ps(GMX_FLOAT_NEGZERO))
#define simdMaxF             _mm512_max_ps
#define simdMinF             _mm512_min_ps
#define simdRoundF(x)        _mm512_roundscale_ps(x, 0)
#define simdTruncF(x)        _mm512_trunc_ps(x)
#define simdFractionF(x)     _mm512_sub_ps(x, simdTruncF(x))
#define simdGetExponentF(x) _mm512_getexp_ps(x)
#define simdGetMantissaF(x) _mm512_getmant_ps(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero)
#define simdSetExponentF(x) simdSetExponentF_x86_avx_512f(x)
/* integer datatype corresponding to float: SimdFInt32 */
#define SimdFInt32          __m512i
#define simdLoadFI           _mm512_load_si512
#define simdSet1FI           _mm512_set1_epi32
#define simdStoreFI          _mm512_store_si512
#define simdLoadUFI          _mm512_loadu_si512
#define simdStoreUFI         _mm512_storeu_si512
#undef  simdExtractFI
#define simdSetZeroFI        _mm512_setzero_epi32
#define simdCvtF2I           _mm512_cvtps_epi32
#define simdCvttF2I          _mm512_cvttps_epi32
#define simdCvtI2F           _mm512_cvtepi32_ps
/* Integer logical ops on SimdFInt32 */
#define simdSlliFI           _mm512_slli_epi32
#define simdSrliFI           _mm512_srli_epi32
#define simdAndFI            _mm512_and_epi32
#define simdAndNotFI         _mm512_andnot_epi32
#define simdOrFI             _mm512_or_epi32
#define simdXorFI            _mm512_xor_epi32
/* Integer arithmetic ops on SimdFInt32 */
#define simdAddFI            _mm512_add_epi32
#define simdSubFI            _mm512_sub_epi32
#define simdMulFI            _mm512_mullo_epi32
/* Boolean & comparison operations on SimdFloat */
#define SimdFBool           __mmask16
#define simdCmpEqF(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_EQ_OQ)
#define simdCmpLtF(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_LT_OS)
#define simdCmpLeF(a, b)     _mm512_cmp_ps_mask(a, b, _CMP_LE_OS)
#define simdAndFB            _mm512_kand
#define simdAndNotFB(a, b)   _mm512_kandn(a, b)
#define simdOrFB             _mm512_kor
#define simdAnyTrueFB        _mm512_mask2int
#define simdMaskF(a, sel)    _mm512_mask_mov_ps(_mm512_setzero_ps(), sel, a)
#define simdMaskNotF(a, sel) _mm512_mask_mov_ps(a, sel, _mm512_setzero_ps())
#define simdBlendF(a, b, sel)    _mm512_mask_blend_ps(sel, a, b)
#define simdReduceF(a)       simdReduceF_x86_avx_512f(a)
/* Boolean & comparison operations on SimdFInt32 */
#define SimdFIBool          __mmask16
#define simdCmpEqFI(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_EQ)
#define simdCmpLtFI(a, b)    _mm512_cmp_epi32_mask(a, b, _MM_CMPINT_LT)
#define simdAndFIB           _mm512_kand
#define simdOrFIB            _mm512_kor
#define simdAnyTrueFIB       _mm512_mask2int
#define simdMaskFI(a, sel)    _mm512_mask_mov_epi32(_mm512_setzero_epi32(), sel, a)
#define simdMaskNotFI(a, sel) _mm512_mask_mov_epi32(a, sel, _mm512_setzero_epi32())
#define simdBlendFI(a, b, sel)    _mm512_mask_blend_epi32(sel, a, b)
/* Conversions between different booleans */
#define simdCvtFB2FIB(x)     (x)
#define simdCvtFIB2FB(x)     (x)


/* Implementation helper functions */

/* This is likely faster than the built in scale operation (lat 8, t-put 3)
 * since we only work on the integer part and use shifts, meaning it will run
 * on the integer ports that are typically less utilized in our kernels.
 */
static inline __m512
simdSetExponentF_x86_avx_512f(__m512 a)
{
    const __m512i expbias      = _mm512_set1_epi32(127);
    __m512i       iexp         = simdCvtF2I(a);

    iexp = _mm512_slli_epi32(_mm512_add_epi32(iexp, expbias), 23);
    return _mm512_castsi512_ps(iexp);
}

static inline float
simdReduceF_x86_avx_512f(__m512 a)
{
    __m128 b;
    a = _mm512_add_ps(a, _mm512_shuffle_f32x4(a, a, _MM_PERM_DCDC));
    a = _mm512_add_ps(a, _mm512_shuffle_f32x4(a, a, _MM_PERM_ABAB));
    b = _mm512_castps512_ps128(a);
    b = _mm_hadd_ps(b, b);
    b = _mm_hadd_ps(b, b);
    return _mm_cvtss_f32(b);
}

#endif /* GMX_SIMD_IMPL_X86_AVX_512F_SIMD_FLOAT_H */
