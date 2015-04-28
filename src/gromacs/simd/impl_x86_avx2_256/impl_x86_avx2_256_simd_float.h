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

#include "gromacs/simd/impl_x86_avx_256/impl_x86_avx_256_simd_float.h"

namespace gmx
{

struct SimdFIBool
{
    __m256i b;
};

static inline SimdFloat gmx_simdcall
simdFmaddF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_fmadd_ps(a.r, b.r, c.r)
    };
}

static inline SimdFloat gmx_simdcall
simdFmsubF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_fmsub_ps(a.r, b.r, c.r)
    };
}

static inline SimdFloat gmx_simdcall
simdFnmaddF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_fnmadd_ps(a.r, b.r, c.r)
    };
}

static inline SimdFloat gmx_simdcall
simdFnmsubF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               _mm256_fnmsub_ps(a.r, b.r, c.r)
    };
}

static inline SimdFloat gmx_simdcall
simdGetExponentF(SimdFloat x)
{
    const __m256  exponentMask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7F800000));
    const __m256i exponentBias = _mm256_set1_epi32(127);
    __m256i       iExponent;

    iExponent = _mm256_castps_si256(_mm256_and_ps(x.r, exponentMask));
    iExponent = _mm256_sub_epi32(_mm256_srli_epi32(iExponent, 23), exponentBias);
    return {
               _mm256_cvtepi32_ps(iExponent)
    };
}

static inline SimdFloat gmx_simdcall
simdSetExponentF(SimdFloat x)
{
    const __m256i  exponentBias = _mm256_set1_epi32(127);
    __m256i        iExponent    = _mm256_cvtps_epi32(x.r);

    iExponent = _mm256_slli_epi32(_mm256_add_epi32(iExponent, exponentBias), 23);
    return {
               _mm256_castsi256_ps(iExponent)
    };
}

static inline SimdFInt32 gmx_simdcall
simdSlliFI(SimdFInt32 a, int n)
{
    return {
               _mm256_slli_epi32(a.i, n)
    };
}

static inline SimdFInt32 gmx_simdcall
simdSrliFI(SimdFInt32 a, int n)
{
    return {
               _mm256_srli_epi32(a.i, n)
    };
}

static inline SimdFInt32 gmx_simdcall
simdAndFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_and_si256(a.i, b.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdAndNotFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_andnot_si256(a.i, b.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdOrFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_or_si256(a.i, b.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdXorFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_xor_si256(a.i, b.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdAddFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_add_epi32(a.i, b.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdSubFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_sub_epi32(a.i, b.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdMulFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_mullo_epi32(a.i, b.i)
    };
}

static inline SimdFIBool gmx_simdcall
simdCmpEqFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_cmpeq_epi32(a.i, b.i)
    };
}

static inline SimdFIBool gmx_simdcall
simdCmpLtFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm256_cmpgt_epi32(b.i, a.i)
    };
}

static inline SimdFIBool gmx_simdcall
simdAndFIB(SimdFIBool a, SimdFIBool b)
{
    return {
               _mm256_and_si256(a.b, b.b)
    };
}

static inline SimdFIBool gmx_simdcall
simdOrFIB(SimdFIBool a, SimdFIBool b)
{
    return {
               _mm256_or_si256(a.b, b.b)
    };
}

static inline bool gmx_simdcall
simdAnyTrueFIB(SimdFIBool a) { return _mm256_movemask_epi8(a.b) != 0; }

static inline SimdFInt32 gmx_simdcall
simdMaskFI(SimdFInt32 a, SimdFIBool mask)
{
    return {
               _mm256_and_si256(a.i, mask.b)
    };
}

static inline SimdFInt32 gmx_simdcall
simdMaskNotFI(SimdFInt32 a, SimdFIBool mask)
{
    return {
               _mm256_andnot_si256(mask.b, a.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdBlendFI(SimdFInt32 a, SimdFInt32 b, SimdFIBool sel)
{
    return {
               _mm256_blendv_epi8(a.i, b.i, sel.b)
    };
}

static inline SimdFIBool gmx_simdcall
simdCvtFB2FIB(SimdFBool a)
{
    return {
               _mm256_castps_si256(a.b)
    };
}

static inline SimdFBool gmx_simdcall
simdCvtFIB2FB(SimdFIBool a)
{
    return {
               _mm256_castsi256_ps(a.b)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX2_256_SIMD_FLOAT_H
