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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_SIMD4_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_256_SIMD4_DOUBLE_H

#include "config.h"

#include <immintrin.h>

#include "impl_x86_avx_256_common.h"
#include "impl_x86_avx_256_simd_double.h"

/****************************************************
 *      DOUBLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
#define Simd4Double           SimdDouble
#define simd4LoadD            simdLoadD
#define simd4Load1D           simdLoad1D
#define simd4Set1D            simdSet1D
#define simd4StoreD           simdStoreD
#define simd4LoadUD           simdLoadUD
#define simd4StoreUD          simdStoreUD
#define simd4SetZeroD         simdSetZeroD
#define simd4AddD             simdAddD
#define simd4SubD             simdSubD
#define simd4MulD             simdMulD
#define simd4FmaddD           simdFmaddD
#define simd4FmsubD           simdFmsubD
#define simd4FnmaddD          simdFnmaddD
#define simd4FnmsubD          simdFnmsubD
#define simd4AndD             simdAndD
#define simd4AndNotD          simdAndNotD
#define simd4OrD              simdOrD
#define simd4XorD             simdXorD
#define simd4RsqrtD           simdRsqrtD
#define simd4AbsD            simdAbsD
#define simd4NegD            simdNegD
#define simd4MaxD             simdMaxD
#define simd4MinD             simdMinD
#define simd4RoundD           simdRoundD
#define simd4TruncD           simdTruncD
#define simd4DotProductD     simd4DotProductD_avx_256
#define Simd4DBool           SimdDBool
#define simd4CmpEqD           simdCmpEqD
#define simd4CmpLtD           simdCmpLtD
#define simd4CmpLeD           simdCmpLeD
#define simd4AndDB            simdAndDB
#define simd4OrDB             simdOrDB
#define simd4AnyTrueDB        simdAnyTrueDB
#define simd4MaskD       simdMaskD
#define simd4MaskNotD    simdMaskNotD
#define simd4BlendD          simdBlendD
#define simd4ReduceD          simdReduceD

/* Implementation helpers */
static inline double gmx_simdcall
simd4DotProductD_avx_256(__m256d a, __m256d b)
{
    double  d;
    __m128d tmp1, tmp2;
    a    = _mm256_mul_pd(a, b);
    tmp1 = _mm256_castpd256_pd128(a);
    tmp2 = _mm256_extractf128_pd(a, 0x1);

    tmp1 = _mm_add_pd(tmp1, _mm_permute_pd(tmp1, _MM_SHUFFLE2(0, 1)));
    tmp1 = _mm_add_pd(tmp1, tmp2);
    _mm_store_sd(&d, tmp1);
    return d;
}


#endif /* GMX_SIMD_IMPL_X86_AVX_256_SIMD4_DOUBLE_H */
