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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD4_DOUBLE_H
#define GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD4_DOUBLE_H

#include "config.h"

#include <math.h>
#ifdef __clang__
#include <qpxmath.h>
#endif

#include "impl_ibm_qpx_common.h"
#include "impl_ibm_qpx_simd_double.h"

/* QPX is already 4-wide both in single and double, so just reuse for SIMD4 */

/* DOUBLE */
#define Simd4Double                SimdDouble
#define simd4LoadD                 simdLoadD
#define simd4Load1D                simdLoad1D
#define simd4Set1D                 simdSet1D
#define simd4StoreD                simdStoreD
#define simd4LoadUD                simdLoadUD
#define simd4StoreUD               simdStoreUD
#define simd4SetZeroD              simdSetZeroD
#define simd4AddD                  simdAddD
#define simd4SubD                  simdSubD
#define simd4MulD                  simdMulD
#define simd4FmaddD                simdFmaddD
#define simd4FmsubD                simdFmsubD
#define simd4FnmaddD               simdFnmaddD
#define simd4FnmsubD               simdFnmsubD
#define simd4AndD                  simdAndD
#define simd4AndNotD               simdAndNotD
#define simd4OrD                   simdOrD
#define simd4XorD                  simdXorD
#define simd4RsqrtD                simdRsqrtD
#define simd4AbsD                 simdAbsD
#define simd4NegD                 simdNegD
#define simd4MaxD                  simdMaxD
#define simd4MinD                  simdMinD
#define simd4RoundD                simdRoundD
#define simd4TruncD                simdTruncD
#define simd4DotProductD          simd4DotProductD_ibm_qpx
#define Simd4DInt32                SimdDInt32
#define simd4LoadDI                simdLoadDI
#define simd4Load1DI               simdLoad1DI
#define simd4Set1DI                simdSet1DI
#define simd4StoreDI               simdStoreDI
#define simd4LoadUDI               simdLoadUDI
#define simd4StoreUDI              simdStoreUDI
#define simd4SetZeroDI             simdSetZeroDI
#define Simd4DBool                SimdDBool
#define simd4CmpEqD                simdCmpEqD
#define simd4CmpLtD                simdCmpLtD
#define simd4CmpLeD                simdCmpLeD
#define simd4AndDB                 simdAndDB
#define simd4OrDB                  simdOrDB
#define simd4AnyTrueDB             simdAnyTrueDB
#define simd4MaskD            simdMaskD
#define simd4MaskNotD         simdMaskNotD
#define simd4BlendD               simdBlendD
#define simd4ReduceD               simdReduceD

static __attribute__((always_inline)) double gmx_simdcall
simd4DotProductD_ibm_qpx(vector4double a, vector4double b)
{
    vector4double dp_sh0 = vec_mul(a, b);
    vector4double dp_sh1 = vec_sldw(dp_sh0, dp_sh0, 1);
    vector4double dp_sh2 = vec_sldw(dp_sh0, dp_sh0, 2);
    vector4double dp     = vec_add(dp_sh2, vec_add(dp_sh0, dp_sh1));

    return vec_extract(dp, 0);
}

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD4_DOUBLE_H */
