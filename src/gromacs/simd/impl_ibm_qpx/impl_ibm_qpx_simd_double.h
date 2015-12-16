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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_DOUBLE_H
#define GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_DOUBLE_H

#include "config.h"

#include <cmath>
#include <cstdint>

#ifdef __clang__
#include <qpxmath.h>
#endif

#include "impl_ibm_qpx_common.h"

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdDouble         vector4double
#ifdef NDEBUG
#    define simdLoadD(m)    vec_ld(0, (double *)(m))
#    define simdStoreD(m, a) vec_st(a, 0, (double *)(m))
#else
#    define simdLoadD(m)    vec_lda(0, (double *)(m))
#    define simdStoreD(m, a) vec_sta(a, 0, (double *)(m))
#endif
#    define simdLoad1D(m)   vec_lds(0, (double *)(m))
#define simdSet1D(x)        vec_splats(x)
/* No support for unaligned load/store */
#define simdSetZeroD        simdSetZeroIbm_qpx
#define simdAddD(a, b)       vec_add(a, b)
#define simdSubD(a, b)       vec_sub(a, b)
#define simdMulD(a, b)       vec_mul(a, b)
#define simdFmaddD(a, b, c)   vec_madd(a, b, c)
#define simdFmsubD(a, b, c)   vec_msub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b+c=-(a*b-c) is "nmsub" */
#define simdFnmaddD(a, b, c)  vec_nmsub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b-c=-(a*b+c) is "nmadd" */
#define simdFnmsubD(a, b, c)  vec_nmadd(a, b, c)
/* simdAndD not supported - no bitwise logical ops */
/* simdAndNotD not supported - no bitwise logical ops */
/* simdOrD not supported - no bitwise logical ops */
/* simdXorD not supported - no bitwise logical ops */
#define simdRsqrtD(a)       vec_rsqrte(a)
#define simdRcpD(a)         vec_re(a)
#define simdAbsD(a)        vec_abs(a)
#define simdNegD           gmx_simd_fneg_ibm_qpx
#define simdMaxD(a, b)       vec_sel(b, a, vec_sub(a, b))
#define simdMinD(a, b)       vec_sel(b, a, vec_sub(b, a))
/* Note: It is critical to use vec_cfid(vec_ctid(a)) for the implementation
 * of simdRoundF(), since vec_round() does not adhere to the FP control
 * word rounding scheme. We rely on float-to-float and float-to-integer
 * rounding being the same for half-way values in a few algorithms.
 */
#define simdRoundD(a)       vec_cfid(vec_ctid(a))
#define simdTruncD(a)       vec_trunc(a)
#define simdFractionD(x)    vec_sub(x, vec_trunc(x))
#define simdGetExponentD(a) gmx_simd_get_exponent_ibm_qpx(a)
#define simdGetMantissaD(a) gmx_simd_get_mantissa_ibm_qpx(a)
#define simdSetExponentD(a) gmx_simd_set_exponent_ibm_qpx(a)
/* integer datatype corresponding to double: SimdDInt32 */
#define SimdDInt32         vector4double
#ifdef NDEBUG
#    define simdLoadDI(m)   vec_ldia(0, (int *)(m))
#else
#    define simdLoadDI(m)   vec_ldiaa(0, (int *)(m))
#endif
#define simdSet1DI(i)       simdSet1Int_ibm_qpx(i)
#define simdStoreDI(m, x)   vec_st(x, 0, (int *)(m))
#define simdSetZeroDI       simdSetZeroIbm_qpx
#define simdCvtD2I(a)       vec_ctiw(a)
#define simdCvttD2I(a)      vec_ctiwz(a)
#define simdCvtI2D(a)       vec_cfid(a)
/* Integer simd extract not available */
/* Integer logical ops on SimdDInt32 not supported */
/* Integer arithmetic ops on SimdDInt32 not supported */
/* Boolean & comparison operations on SimdDouble */
#define SimdDBool          vector4double
#define simdCmpEqD(a, b)     vec_cmpeq(a, b)
#define simdCmpLtD(a, b)     vec_cmplt((a), (b))
#define simdCmpLeD(a, b)     simdOrFB(vec_cmpeq(a, b), vec_cmplt(a, b))
#define simdAndDB(a, b)      vec_and(a, b)
#define simdOrDB(a, b)       vec_or(a, b)
#define simdAnyTrueDB(a)     simdAnyTrueB_ibm_qpx(a)
#define simdMaskD(a, sel)   vec_sel(vec_splats(0.0), a, sel)
#define simdMaskNotD(a, sel) vec_sel(a, vec_splats(0.0), sel)
#define simdBlendD(a, b, sel)  vec_sel(a, b, sel)
#define simdReduceD(a)      gmx_simd_reduce_ibm_qpx(a)

/* Boolean & comparison operations on SimdDInt32 not supported */
/* Conversions between different booleans not supported */


/****************************************************
 * IMPLEMENTATION HELPER FUNCTIONS                  *
 ****************************************************/
static __attribute__((always_inline)) vector4double
gmx_simd_fneg_ibm_qpx(vector4double a)
{
    return vec_neg(a);
}

static __attribute__((always_inline)) vector4double gmx_simdcall
simdSetZeroIbm_qpx(void)
{
    return vec_splats(0.0);
}

static __attribute__((always_inline)) vector4double gmx_simdcall
gmx_simd_get_exponent_ibm_qpx(vector4double x)
{
    const std::int64_t    expmask   = 0x7ff0000000000000LL;
    const std::int64_t    expbase   = 1023;
    std::int64_t          idata[4] __attribute__((aligned(32)));

    /* Store to memory */
    vec_st(x, 0, idata);
    /* Perform integer arithmetics in general registers. */
    idata[0] = ((idata[0] & expmask) >> 52) - expbase;
    idata[1] = ((idata[1] & expmask) >> 52) - expbase;
    idata[2] = ((idata[2] & expmask) >> 52) - expbase;
    idata[3] = ((idata[3] & expmask) >> 52) - expbase;
    /* Reload from memory */
    return vec_cfid(vec_ld(0, idata));
}

static __attribute__((always_inline)) vector4double gmx_simdcall
gmx_simd_get_mantissa_ibm_qpx(vector4double x)
{
    const std::int64_t    exp_and_sign_mask = 0xfff0000000000000LL;
    const std::int64_t    ione              = 0x3ff0000000000000LL;
    std::int64_t          idata[4] __attribute__((aligned(32)));

    /* Store to memory */
    vec_st(x, 0, idata);
    /* Perform integer arithmetics in general registers. */
    idata[0] = (idata[0] & (~exp_and_sign_mask)) | ione;
    idata[1] = (idata[1] & (~exp_and_sign_mask)) | ione;
    idata[2] = (idata[2] & (~exp_and_sign_mask)) | ione;
    idata[3] = (idata[3] & (~exp_and_sign_mask)) | ione;
    /* Reload from memory */
    return vec_ld(0, idata);
}

static __attribute__((always_inline)) vector4double gmx_simdcall
gmx_simd_set_exponent_ibm_qpx(vector4double x)
{
    const std::int64_t    expbase = 1023;
    std::int64_t          idata[4] __attribute__((aligned(32)));

    /* Store to memory for shifts. It is REALLY critical that we use the same
     * rounding mode as for simdRound() here. In particular, for QPX
     * this means we implement simdRound(a) as vec_cfid(vec_ctid(a)),
     * since vec_round() uses a different rounding scheme.
     */
    vec_st(vec_ctid(x), 0, idata);
    /* Perform integer arithmetics in general registers. */
    idata[0] = (idata[0] + expbase) << 52;
    idata[1] = (idata[1] + expbase) << 52;
    idata[2] = (idata[2] + expbase) << 52;
    idata[3] = (idata[3] + expbase) << 52;
    /* Reload from memory */
    return vec_ld(0, idata);
}

static __attribute__((always_inline)) double gmx_simdcall
gmx_simd_reduce_ibm_qpx(vector4double x)
{
    vector4double y = vec_sldw(x, x, 2);
    vector4double z;

    y = vec_add(y, x);
    z = vec_sldw(y, y, 1);
    y = vec_add(y, z);
    return vec_extract(y, 0);
}

static __attribute__((always_inline)) vector4double gmx_simdcall
simdSet1Int_ibm_qpx(int i)
{
    int idata[4] __attribute__((aligned(32)));

    idata[0] = i;

    /* Reload from memory */
    return vec_splat(vec_ldia(0, idata), 0);
}

/* This works in both single and double */
static __attribute__((always_inline)) int gmx_simdcall
simdAnyTrueB_ibm_qpx(vector4double a)
{
    vector4double b = vec_sldw(a, a, 2);

    a = vec_or(a, b);
    b = vec_sldw(a, a, 1);
    a = vec_or(a, b);
    return (vec_extract(a, 0) > 0);
}

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_DOUBLE_H */
