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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VMX_SIMD_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VMX_SIMD_FLOAT_H

#include <math.h>

#include <altivec.h>

#include "impl_ibm_vmx_common.h"

/* Make sure we do not screw up c++ - undefine vector/bool, and rely on __vector and __bool */
#undef vector
#undef bool

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdFloat            __vector float
#define simdLoadF(m)         vec_ld(0, (const __vector float *)(m))
#define simdLoad1F(m)        simdLoad1F_ibm_vmx(m)
#define simdSet1F(x)         simdSet1F_ibm_vmx(x)
#define simdStoreF(m, x)     vec_st(x, 0, (__vector float *)(m))
#undef  simdLoadUF
#undef  simdStoreUF
#define simdSetZeroF()       ((__vector float)vec_splat_u32(0))
#define simdAddF(a, b)       vec_add(a, b)
#define simdSubF(a, b)       vec_sub(a, b)
#define simdMulF(a, b)       vec_mul(a, b)
#define simdFmaddF(a, b, c)  vec_madd(a, b, c)
#define simdFmsubF(a, b, c)  vec_sub(vec_mul(a, b), c)
/* IBM uses an alternative FMA definition, so -a*b+c=-(a*b-c) is "nmsub" */
#define simdFnmaddF(a, b, c) vec_nmsub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b-c=-(a*b+c) is "nmadd" */
#define simdFnmsubF(a, b, c) vec_sub(simdSetZeroF(), vec_madd(a, b, c))
#define simdAndF(a, b)       vec_and(a, b)
#define simdAndNotF(a, b)    vec_andc(b, a)
#define simdOrF(a, b)        vec_or(a, b)
#define simdXorF(a, b)       vec_xor(a, b)
#define simdRsqrtF(a)        vec_rsqrte(a)
#define simdRcpF(a)          vec_re(a)
#define simdAbsF(a)         vec_abs(a)
#define simdNegF(a)         vec_xor(a, (__vector float)vec_sl(vec_splat_u32(-1), vec_splat_u32(-1)))
#define simdMaxF(a, b)       vec_max(a, b)
#define simdMinF(a, b)       vec_min(a, b)
#define simdRoundF(a)        vec_round(a)
#define simdTruncF(a)        vec_trunc(a)
#define simdFractionF(x)     vec_sub(x, vec_trunc(x))
#define simdGetExponentF(a) simdGetExponentF_ibm_vmx(a)
#define simdGetMantissaF(a) simdGetMantissaF_ibm_vmx(a)
#define simdSetExponentF(a) simdSetExponentF_ibm_vmx(a)
/* integer datatype corresponding to float: SimdFInt32 */
#define SimdFInt32          __vector int
#define simdLoadFI(m)        vec_ld(0, (const __vector int *)(m))
#define simdSet1FI(i)        simdSet1FI_ibm_vmx((int)(i))
#define simdStoreFI(m, x)    vec_st(x, 0, (__vector int *)(m))
#undef  simdLoadUFI
#undef  simdStoreUFI
#define simdSetZeroFI()      vec_splat_s32(0)
#define simdCvtF2I(a)        vec_cts(vec_round(a), 0)
#define simdCvttF2I(a)       vec_cts(a, 0)
#define simdCvtI2F(a)        vec_ctf(a, 0)
#undef  simdExtractFI
/* Integer logical ops on SimdFInt32 */
/* The shift constant magic requires an explanation:
 * VMX only allows literals up to 15 to be created directly with vec_splat_u32,
 * and we need to be able to shift up to 31 bits. The code on the right hand
 * side splits the constant in three parts with values in the range 0..15.
 * Since the argument has to be a constant (but our and VMX requirement), these
 * constants will be evaluated at compile-time, and if one or two parts evaluate
 * to zero they will be removed with -O2 or higher optimization (checked).
 */
#define simdSlliFI(a, i)      vec_sl(a, vec_add(vec_add(vec_splat_u32( (((i&0xF)+(i/16))&0xF)+i/31 ), vec_splat_u32( (i/16)*15 )), vec_splat_u32( (i/31)*15 )))
#define simdSrliFI(a, i)      vec_sr(a, vec_add(vec_add(vec_splat_u32( (((i&0xF)+(i/16))&0xF)+i/31 ), vec_splat_u32( (i/16)*15 )), vec_splat_u32( (i/31)*15 )))
#define simdAndFI(a, b)       vec_and(a, b)
#define simdAndNotFI(a, b)   vec_andc(b, a)
#define simdOrFI(a, b)        vec_or(a, b)
#define simdXorFI(a, b)       vec_xor(a, b)
/* Integer arithmetic ops on SimdFInt32 */
#define simdAddFI(a, b)       vec_add(a, b)
#define simdSubFI(a, b)       vec_sub(a, b)
#define simdMulFI(a, b)       vec_mule((__vector short)a, (__vector short)b)
/* Boolean & comparison operations on SimdFloat */
#define SimdFBool           __vector __bool int
#define simdCmpEqF(a, b)     vec_cmpeq(a, b)
#define simdCmpLtF(a, b)     vec_cmplt(a, b)
#define simdCmpLeF(a, b)     vec_cmple(a, b)
#define simdAndFB(a, b)      vec_and(a, b)
#define simdOrFB(a, b)       vec_or(a, b)
#define simdAnyTrueFB(a)     vec_any_ne(a, (__vector __bool int)vec_splat_u32(0))
#define simdMaskF(a, sel)    vec_and(a, (__vector float)sel)
#define simdMaskNotF(a, sel) vec_andc(a, (__vector float)sel)
#define simdBlendF(a, b, sel)    vec_sel(a, b, sel)
#define simdReduceF(a)       simdReduceF_ibm_vmx(a)
/* Boolean & comparison operations on SimdFInt32 */
#define SimdFIBool          __vector __bool int
#define simdCmpEqFI(a, b)     vec_cmpeq(a, b)
#define simdCmpLtFI(a, b)     vec_cmplt(a, b)
#define simdAndFIB(a, b)      vec_and(a, b)
#define simdOrFIB(a, b)       vec_or(a, b)
#define simdAnyTrueFIB(a)          vec_any_ne(a, (__vector __bool int)vec_splat_u32(0))
#define simdMaskFI(a, sel)    vec_and(a, (__vector int)sel)
#define simdMaskNotFI(a, sel) vec_andc(a, (__vector int)sel)
#define simdBlendFI(a, b, sel)    vec_sel(a, b, sel)
/* Conversions between different booleans */
#define simdCvtFB2FIB(x)     (x)
#define simdCvtFIB2FB(x)     (x)

/* Double is not available with VMX SIMD */

/****************************************************
 * IMPLEMENTATION HELPER FUNCTIONS                  *
 ****************************************************/
static inline SimdFloat
simdSet1F_ibm_vmx(const float x)
{
    /* In the old days when PPC was all big endian we could
     * use the vec_lvsl() instruction to permute bytes based on
     * a load adress. However, at least with gcc-4.8.2 the bytes
     * end up reversed on Power8 running little endian (Linux).
     * Since this is not a critical instruction we work around
     * it by first putting the data in an aligned position before
     * loading, so we can avoid vec_lvsl() entirely. We can
     * do this slightly faster on GCC with alignment attributes.
     */
    __vector float vx;
#ifdef __GNUC__
    float alignedx __attribute ((aligned (16)));
    alignedx = x;
    vx       = vec_lde(0, &alignedx);
#else
    struct {
        vector float vx; float x[4];
    } conv;
    conv.x[0] = x;
    vx        = vec_lde(0, conv.x);
#endif
    return vec_splat(vx, 0);
}

static inline SimdFloat
simdLoad1F_ibm_vmx(const float * m)
{
    return simdSet1F_ibm_vmx(*m);
}

static inline SimdFInt32
simdSet1FI_ibm_vmx(const int x)
{
    /* See comment in simdSet1F_ibm_vmx why we
     * cannot use vec_lvsl().
     */
    __vector int vx;
#ifdef __GNUC__
    int alignedx __attribute ((aligned (16)));
    alignedx = x;
    vx       = vec_lde(0, &alignedx);
#else
    struct {
        vector int vx; int x[4];
    } conv;
    conv.x[0] = x;
    vx        = vec_lde(0, conv.x);
#endif
    return vec_splat(vx, 0);
}


static inline SimdFloat
simdGetExponentF_ibm_vmx(SimdFloat x)
{
    /* Generate 0x7F800000 without memory operations */
    SimdFloat     expmask = (__vector float)simdSlliFI(vec_add(vec_splat_s32(15), vec_sl(vec_splat_s32(15), vec_splat_u32(4))), 23);
    SimdFInt32    i127    = vec_sub(vec_sl(vec_splat_s32(1), vec_splat_u32(7)), vec_splat_s32(1));
    SimdFInt32    iexp;

    iexp = (__vector int)simdAndF(x, expmask);
    iexp = vec_sub(simdSrliFI(iexp, 23), i127);
    return vec_ctf(iexp, 0);
}

static inline SimdFloat
simdGetMantissaF_ibm_vmx(SimdFloat x)
{
    SimdFloat     expmask = (__vector float)simdSlliFI(vec_add(vec_splat_s32(15), vec_sl(vec_splat_s32(15), vec_splat_u32(4))), 23);

    /* Get mantissa. By taking the absolute value (to get rid of the sign bit) we can
     * use the same mask as for simdGetExponentF() (but complement it). Since
     * these two routines are typically called together, this will save a few operations.
     */
    x = simdAndNotF(expmask, vec_abs(x));
    /* Reset zero (but correctly biased) exponent */
    return simdOrF(x, vec_ctf(vec_splat_s32(1), 0));
}

static inline SimdFloat
simdSetExponentF_ibm_vmx(SimdFloat x)
{
    SimdFInt32  iexp = simdCvtF2I(x);
    SimdFInt32  i127 = vec_sub(vec_sl(vec_splat_s32(1), vec_splat_u32(7)), vec_splat_s32(1));

    iexp = simdSlliFI(vec_add(iexp, i127), 23);
    return (__vector float)iexp;
}

static inline float
simdReduceF_ibm_vmx(SimdFloat x)
{
    float res;
    x = vec_add(x, vec_sld(x, x, 8));
    x = vec_add(x, vec_sld(x, x, 4));
    vec_ste(x, 0, &res);
    return res;
}

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_VMX_SIMD_FLOAT_H */
