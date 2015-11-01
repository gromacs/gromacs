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

#ifndef GMX_SIMD_IMPL_X86_AVX2_256_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX2_256_SIMD_FLOAT_H

#include "config.h"

#include <immintrin.h>

#include "impl_x86_avx2_256_common.h"

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#undef  simdFmaddF
#define simdFmaddF           _mm256_fmadd_ps
#undef  simdFmsubF
#define simdFmsubF           _mm256_fmsub_ps
#undef  simdFnmaddF
#define simdFnmaddF          _mm256_fnmadd_ps
#undef  simdFnmsubF
#define simdFnmsubF          _mm256_fnmsub_ps
#undef  simdGetExponentF
#define simdGetExponentF    simdGetExponentF_avx2_256
#undef  simdSetExponentF
#define simdSetExponentF    simdSetExponentF_avx2_256
/* Previously undefined logical ops on SimdFInt32 */
#define simdSlliFI           _mm256_slli_epi32
#define simdSrliFI           _mm256_srli_epi32
#define simdAndFI            _mm256_and_si256
#define simdAndNotFI         _mm256_andnot_si256
#define simdOrFI             _mm256_or_si256
#define simdXorFI            _mm256_xor_si256
/* Previously undefined arithmetic ops on SimdFInt32 */
#define simdAddFI            _mm256_add_epi32
#define simdSubFI            _mm256_sub_epi32
#define simdMulFI            _mm256_mullo_epi32
/* Previously undefined boolean ops on SimdFInt32 */
#define simdCmpEqFI          _mm256_cmpeq_epi32
#define simdCmpLtFI(a, b)     _mm256_cmpgt_epi32(b, a)
#define simdAndFIB           _mm256_and_si256
#define simdOrFIB            _mm256_or_si256
#define simdAnyTrueFIB       _mm256_movemask_epi8
#define simdMaskFI      _mm256_and_si256
#define simdMaskNotFI(a, sel) _mm256_andnot_si256(sel, a)
#define simdBlendFI         _mm256_blendv_epi8

/*********************************************************
 * SIMD SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 *********************************************************/
static inline SimdFloat gmx_simdcall
simdGetExponentF_avx2_256(SimdFloat x)
{
    const __m256  expmask      = _mm256_castsi256_ps(_mm256_set1_epi32(0x7F800000));
    const __m256i expbias      = _mm256_set1_epi32(127);
    __m256i       iexp;

    iexp = _mm256_castps_si256(_mm256_and_ps(x, expmask));
    iexp = _mm256_sub_epi32(_mm256_srli_epi32(iexp, 23), expbias);
    return _mm256_cvtepi32_ps(iexp);
}

static inline SimdFloat gmx_simdcall
simdSetExponentF_avx2_256(SimdFloat x)
{
    const __m256i  expbias      = _mm256_set1_epi32(127);
    __m256i        iexp         = _mm256_cvtps_epi32(x);

    iexp = _mm256_slli_epi32(_mm256_add_epi32(iexp, expbias), 23);
    return _mm256_castsi256_ps(iexp);
}

#endif /* GMX_SIMD_IMPL_X86_AVX2_256_SIMD_FLOAT_H */
