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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_FLOAT_H

#include <math.h>
#ifdef __clang__
#include <qpxmath.h>
#endif

#include "impl_ibm_qpx_common.h"

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdFloat          vector4double
#ifdef NDEBUG
#    define simdLoadF(m)    vec_ld(0, (float *)(m))
#    define simdStoreF(m, a) vec_st(a, 0, (float *)(m))
#else
#    define simdLoadF(m)    vec_lda(0, (float *)(m))
#    define simdStoreF(m, a) vec_sta(a, 0, (float *)(m))
#endif
#    define simdLoad1F(m)   vec_lds(0, (float *)(m))
#define simdSet1F(x)        vec_splats(x)
/* No support for unaligned load/store */
#define simdSetZeroF        simdSetZeroIbm_qpx
#define simdAddF(a, b)       vec_add(a, b)
#define simdSubF(a, b)       vec_sub(a, b)
#define simdMulF(a, b)       vec_mul(a, b)
#define simdFmaddF(a, b, c)   vec_madd(a, b, c)
#define simdFmsubF(a, b, c)   vec_msub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b+c=-(a*b-c) is "nmsub" */
#define simdFnmaddF(a, b, c)  vec_nmsub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b-c=-(a*b+c) is "nmadd" */
#define simdFnmsubF(a, b, c)  vec_nmadd(a, b, c)
/* simdAndF not supported - no bitwise logical ops */
/* simdAndNotF not supported - no bitwise logical ops */
/* simdOrF not supported - no bitwise logical ops */
/* simdXorF not supported - no bitwise logical ops */
#define simdRsqrtF(a)       vec_rsqrte(a)
#define simdRcpF(a)         vec_re(a)
#define simdAbsF(a)        vec_abs(a)
#define simdNegF           gmx_simd_fneg_ibm_qpx
#define simdMaxF(a, b)       vec_sel(b, a, vec_sub(a, b))
#define simdMinF(a, b)       vec_sel(b, a, vec_sub(b, a))
/* Note: It is critical to use vec_cfid(vec_ctid(a)) for the implementation
 * of simdRoundF(), since vec_round() does not adhere to the FP control
 * word rounding scheme. We rely on float-to-float and float-to-integer
 * rounding being the same for half-way values in a few algorithms.
 */
#define simdRoundF(a)       vec_cfid(vec_ctid(a))
#define simdTruncF(a)       vec_trunc(a)
#define simdFractionF(x)    vec_sub(x, vec_trunc(x))
#define simdGetExponentF(a) gmx_simd_get_exponent_ibm_qpx(a)
#define simdGetMantissaF(a) gmx_simd_get_mantissa_ibm_qpx(a)
#define simdSetExponentF(a) gmx_simd_set_exponent_ibm_qpx(a)
/* integer datatype corresponding to float: SimdFInt32 */
#define SimdFInt32         vector4double
#ifdef NDEBUG
#    define simdLoadFI(m)   vec_ldia(0, (int *)(m))
#else
#    define simdLoadFI(m)   vec_ldiaa(0, (int *)(m))
#endif
#define simdSet1FI(i)       simdSet1Int_ibm_qpx(i)
#define simdStoreFI(m, x)    vec_st(x, 0, (int *)(m))
#define simdSetZeroFI       simdSetZeroIbm_qpx
#define simdCvtF2I(a)       vec_ctiw(a)
#define simdCvttF2I(a)      vec_ctiwz(a)
#define simdCvtI2F(a)       vec_cfid(a)
/* Integer simd extract not available */
/* Integer logical ops on SimdFInt32 not supported */
/* Integer arithmetic ops on SimdFInt32 not supported */
/* Boolean & comparison operations on SimdFloat */
#define SimdFBool          vector4double
#define simdCmpEqF(a, b)     vec_cmpeq(a, b)
#define simdCmpLtF(a, b)     vec_cmplt((a), (b))
#define simdCmpLeF(a, b)     simdOrFB(vec_cmpeq(a, b), vec_cmplt(a, b))
#define simdAndFB(a, b)      vec_and(a, b)
#define simdOrFB(a, b)       vec_or(a, b)
#define simdAnyTrueFB(a)    simdAnyTrueBool_ibm_qpx(a)
#define simdMaskF(a, sel) vec_sel(vec_splats(0.0), a, sel)
#define simdMaskNotF(a, sel) vec_sel(a, vec_splats(0.0), sel)
#define simdBlendF(a, b, sel)  vec_sel(a, b, sel)
#define simdReduceF(a)       gmx_simd_reduce_ibm_qpx(a)


/* Boolean & comparison operations on SimdFInt32 not supported */
/* Conversions between different booleans not supported */

/* Note: Since QPX registers are always double internally, the single
 * precision defines use several double precision helper functions.
 */

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_FLOAT_H */
