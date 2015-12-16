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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_256_SIMD_FLOAT_H

#include "config.h"

#include <immintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_avx_256_common.h"

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdFloat           __m256
#define simdLoadF            _mm256_load_ps
#define simdLoad1F           _mm256_broadcast_ss
#define simdSet1F            _mm256_set1_ps
#define simdStoreF           _mm256_store_ps
#define simdLoadUF           _mm256_loadu_ps
#define simdStoreUF          _mm256_storeu_ps
#define simdSetZeroF         _mm256_setzero_ps
#define simdAddF             _mm256_add_ps
#define simdSubF             _mm256_sub_ps
#define simdMulF             _mm256_mul_ps
#define simdFmaddF(a, b, c)    _mm256_add_ps(_mm256_mul_ps(a, b), c)
#define simdFmsubF(a, b, c)    _mm256_sub_ps(_mm256_mul_ps(a, b), c)
#define simdFnmaddF(a, b, c)   _mm256_sub_ps(c, _mm256_mul_ps(a, b))
#define simdFnmsubF(a, b, c)   _mm256_sub_ps(_mm256_setzero_ps(), simdFmaddF(a, b, c))
#define simdAndF             _mm256_and_ps
#define simdAndNotF          _mm256_andnot_ps
#define simdOrF              _mm256_or_ps
#define simdXorF             _mm256_xor_ps
#define simdRsqrtF           _mm256_rsqrt_ps
#define simdRcpF             _mm256_rcp_ps
#define simdAbsF(x)         _mm256_andnot_ps(_mm256_set1_ps(GMX_FLOAT_NEGZERO), x)
#define simdNegF(x)         _mm256_xor_ps(x, _mm256_set1_ps(GMX_FLOAT_NEGZERO))
#define simdMaxF             _mm256_max_ps
#define simdMinF             _mm256_min_ps
#define simdRoundF(x)        _mm256_round_ps(x, _MM_FROUND_NINT)
#define simdTruncF(x)        _mm256_round_ps(x, _MM_FROUND_TRUNC)
#define simdFractionF(x)     _mm256_sub_ps(x, simdTruncF(x))
#define simdGetExponentF    simdGetExponentF_avx_256
#define simdGetMantissaF    simdGetMantissaF_avx_256
#define simdSetExponentF    simdSetExponentF_avx_256
/* integer datatype corresponding to float: SimdFInt32 */
#define SimdFInt32          __m256i
#define simdLoadFI(m)        _mm256_load_si256((__m256i const*)(m))
#define simdSet1FI           _mm256_set1_epi32
#define simdStoreFI(m, x)    _mm256_store_si256((__m256i *)(m), x)
#define simdLoadUFI(m)       _mm256_loadu_si256((__m256i const*)(m))
#define simdStoreUFI(m, x)   _mm256_storeu_si256((__m256i *)(m), x)
#define simdSetZeroFI        _mm256_setzero_si256
#define simdCvtF2I           _mm256_cvtps_epi32
#define simdCvttF2I          _mm256_cvttps_epi32
#define simdCvtI2F           _mm256_cvtepi32_ps
#define simdExtractFI(x, i)   _mm_extract_epi32(_mm256_extractf128_si256(x, (i)>>2), (i)&0x3)
/* Integer logical ops on SimdFInt32 */
/* simdAddFI not supported     */
/* simdSubFI not supported     */
/* simdMulFI not supported     */
/* simdSlliFI not supported    */
/* simdSrliFI not supported    */
/* simdAndFI not supported     */
/* simdAndNotFI not supported  */
/* simdOrFI not supported      */
/* simdXorFI not supported     */
/* Integer arithmetic ops on SimdFInt32 */
/* simdAddFI not supported     */
/* simdSubFI not supported     */
/* simdMulFI not supported     */
/* Boolean & comparison operations on SimdFloat */
#define SimdFBool           __m256
#define simdCmpEqF(a, b)      _mm256_cmp_ps(a, b, _CMP_EQ_OQ)
#define simdCmpLtF(a, b)      _mm256_cmp_ps(a, b, _CMP_LT_OQ)
#define simdCmpLeF(a, b)      _mm256_cmp_ps(a, b, _CMP_LE_OQ)
#define simdAndFB            _mm256_and_ps
#define simdOrFB             _mm256_or_ps
#define simdAnyTrueFB        _mm256_movemask_ps
#define simdMaskF       _mm256_and_ps
#define simdMaskNotF(a, sel)  _mm256_andnot_ps(sel, a)
#define simdBlendF          _mm256_blendv_ps
#define simdReduceF          simdReduceF_avx_256
/* Boolean & comparison operations on SimdFInt32 */
#define SimdFIBool          __m256i
/* simdCmpEqFI not supported        */
/* simdCmpLtFI not supported        */
/* simdAndFIB not supported         */
/* simdOrFIB not supported          */
/* simdAnyTrueFIB not supported     */
/* simdMaskFI not supported    */
/* simdMaskNotFI not supported    */
/* simdBlendFI not supported       */
/* Conversions between different booleans */
#define simdCvtFB2FIB        _mm256_castps_si256
#define simdCvtFIB2FB        _mm256_castsi256_ps

/*********************************************************
 * SIMD SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 *********************************************************/
static inline __m256 gmx_simdcall
simdGetExponentF_avx_256(__m256 x)
{
    const __m256  expmask      = _mm256_castsi256_ps(_mm256_set1_epi32(0x7F800000));
    const __m128i expbias      = _mm_set1_epi32(127);
    __m256i       iexp256;
    __m128i       iexp128a, iexp128b;

    iexp256   = _mm256_castps_si256(_mm256_and_ps(x, expmask));
    iexp128b  = _mm256_extractf128_si256(iexp256, 0x1);
    iexp128a  = _mm256_castsi256_si128(iexp256);
    iexp128a  = _mm_srli_epi32(iexp128a, 23);
    iexp128b  = _mm_srli_epi32(iexp128b, 23);
    iexp128a  = _mm_sub_epi32(iexp128a, expbias);
    iexp128b  = _mm_sub_epi32(iexp128b, expbias);
    iexp256   = _mm256_castsi128_si256(iexp128a);
    iexp256   = _mm256_insertf128_si256(iexp256, iexp128b, 0x1);
    return _mm256_cvtepi32_ps(iexp256);
}

static inline __m256 gmx_simdcall
simdGetMantissaF_avx_256(__m256 x)
{
    const __m256 mantmask   = _mm256_castsi256_ps(_mm256_set1_epi32(0x007FFFFF));
    const __m256 one        = _mm256_set1_ps(1.0);

    x = _mm256_and_ps(x, mantmask);
    return _mm256_or_ps(x, one);
}

static inline __m256 gmx_simdcall
simdSetExponentF_avx_256(__m256 x)
{
    const __m128i expbias      = _mm_set1_epi32(127);
    __m256i       iexp256;
    __m128i       iexp128a, iexp128b;

    iexp256   = _mm256_cvtps_epi32(x);
    iexp128b  = _mm256_extractf128_si256(iexp256, 0x1);
    iexp128a  = _mm256_castsi256_si128(iexp256);
    iexp128a  = _mm_slli_epi32(_mm_add_epi32(iexp128a, expbias), 23);
    iexp128b  = _mm_slli_epi32(_mm_add_epi32(iexp128b, expbias), 23);
    iexp256   = _mm256_castsi128_si256(iexp128a);
    iexp256   = _mm256_insertf128_si256(iexp256, iexp128b, 0x1);
    return _mm256_castsi256_ps(iexp256);
}

static inline float gmx_simdcall
simdReduceF_avx_256(__m256 a)
{
    float  f;

    __m128 a0, a1;
    a  = _mm256_hadd_ps(a, a);
    a  = _mm256_hadd_ps(a, a);
    a0 = _mm256_castps256_ps128(a);
    a1 = _mm256_extractf128_ps(a, 0x1);
    a0 = _mm_add_ss(a0, a1);
    _mm_store_ss(&f, a0);
    return f;
}

#endif /* GMX_SIMD_IMPL_X86_AVX_256_SIMD_FLOAT_H */
