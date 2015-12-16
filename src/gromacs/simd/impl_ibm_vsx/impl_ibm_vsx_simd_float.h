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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD_FLOAT_H

#include <math.h>

#include <altivec.h>

#include "impl_ibm_vsx_common.h"

/* IBM VSX SIMD instruction wrappers. Power7 and later.
 *
 * While this instruction set is similar to VMX, there are quite a few differences
 * that make it easier to understand if we start from scratch rather than by
 * including the VMX implementation and changing lots of things.
 */


/* Make sure we do not screw up c++ - undefine vector/bool, and rely on __vector,
 * which is present both on gcc and xlc.
 */
#undef vector

/* g++ is also unhappy with the clash of vector bool and the C++ reserved 'bool',
 * which is solved by undefining bool and reyling on __bool. However, that does
 * not work with xlc, which requires us to use bool. Solve the conflict by
 * defining a new vsx_bool.
 */
#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__)
#    define vsx_bool __bool
#    undef  bool
#else
#    define vsx_bool bool
#endif

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdFloat           __vector float
#define simdLoadF(m)         (*(const SimdFloat *)(m))
#define simdStoreF(m, x)      { *(SimdFloat *)(m) = (x); }
#define simdLoad1F(m)        vec_splats((float)(*m))
#define simdSet1F(x)         vec_splats((float)(x))
#if defined(__ibmxl__) || defined(__xlC__)
#    define simdLoadUF(m)    vec_xlw4(0, (float *)(m))
#    define simdStoreUF(m, x) vec_xstw4(x, 0, (m))
#else
/* GCC can handle unaligned load/store as pointer dereference */
#    define simdLoadUF       simdLoadF
#    define simdStoreUF      simdStoreF
#endif
#define simdSetZeroF()       vec_splats(0.0f)
#define simdAddF(a, b)       vec_add(a, b)
#define simdSubF(a, b)       vec_sub(a, b)
#define simdMulF(a, b)       vec_mul(a, b)
#define simdFmaddF(a, b, c)  vec_madd(a, b, c)
#define simdFmsubF(a, b, c)  vec_msub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b+c=-(a*b-c) is "nmsub" */
#define simdFnmaddF(a, b, c) vec_nmsub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b-c=-(a*b+c) is "nmadd" */
#define simdFnmsubF(a, b, c) vec_nmadd(a, b, c)
#define simdAndF(a, b)       vec_and(a, b)
#define simdAndNotF(a, b)    vec_andc(b, a)
#define simdOrF(a, b)        vec_or(a, b)
#define simdXorF(a, b)       vec_xor(a, b)
#define simdRsqrtF(a)        vec_rsqrte(a)
#define simdRcpF(a)          vec_re(a)
#define simdAbsF(a)         vec_abs(a)
#define simdNegF(a)         (-(a))
#define simdMaxF(a, b)       vec_max(a, b)
#define simdMinF(a, b)       vec_min(a, b)
#define simdRoundF(a)        vec_round(a)
#define simdTruncF(a)        vec_trunc(a)
#define simdFractionF(x)     vec_sub(x, vec_trunc(x))
#define simdGetExponentF(a) simdGetExponentF_ibm_vsx(a)
#define simdGetMantissaF(a) simdGetMantissaF_ibm_vsx(a)
#define simdSetExponentF(a) simdSetExponentF_ibm_vsx(a)
/* integer datatype corresponding to float: SimdFInt32 */
#define SimdFInt32          __vector signed int
#define simdLoadFI(m)        (*(const SimdFInt32 *)(m))
#define simdStoreFI(m, x)     { *(SimdFInt32 *)(m) = (x); }
#define simdSet1FI(i)        vec_splats((int)(i))
#if defined(__ibmxl__) || defined(__xlC__)
#    define simdLoadUFI(m)    vec_xlw4(0, (int *)(m))
#    define simdStoreUFI(m, x) vec_xstw4(x, 0, (m))
#else
/* GCC can handle unaligned load/store as pointer dereference */
#    define simdLoadUFI      simdLoadFI
#    define simdStoreUFI     simdStoreFI
#endif
#define simdSetZeroFI()      vec_splats((int)0)
#define simdCvtF2I(a)        vec_cts(vec_round(a), 0)
#define simdCvttF2I(a)       vec_cts(a, 0)
#define simdCvtI2F(a)        vec_ctf(a, 0)
#define simdExtractFI(a, i)   simdExtractFI_ibm_vsx(a, i)
/* Integer logical ops on SimdFInt32 */
#define simdSlliFI(a, i)      vec_sl(a, vec_splats((unsigned int)(i)))
#define simdSrliFI(a, i)      vec_sr(a, vec_splats((unsigned int)(i)))
#define simdAndFI(a, b)       vec_and(a, b)
#define simdAndNotFI(a, b)    vec_andc(b, a)
#define simdOrFI(a, b)        vec_or(a, b)
#define simdXorFI(a, b)       vec_xor(a, b)
/* Integer arithmetic ops on SimdFInt32 */
#define simdAddFI(a, b)       vec_add(a, b)
#define simdSubFI(a, b)       vec_sub(a, b)
#define simdMulFI(a, b)      ((a)*(b))
/* Boolean & comparison operations on SimdFloat */
#define SimdFBool           __vector vsx_bool int
#define simdCmpEqF(a, b)     vec_cmpeq(a, b)
#define simdCmpLtF(a, b)     vec_cmplt(a, b)
#define simdCmpLeF(a, b)     vec_cmple(a, b)
#define simdAndFB(a, b)      vec_and(a, b)
#define simdOrFB(a, b)       vec_or(a, b)
#define simdAnyTrueFB(a)     vec_any_ne(a, (__vector vsx_bool int)vec_splats(0))
#define simdMaskF(a, sel)    vec_and(a, (__vector float)sel)
#define simdMaskNotF(a, sel) vec_andc(a, (__vector float)sel)
#define simdBlendF(a, b, sel)    vec_sel(a, b, sel)
#define simdReduceF(a)       simdReduceF_ibm_vsx(a)
/* Boolean & comparison operations on SimdFInt32 */
#define SimdFIBool          __vector vsx_bool int
#define simdCmpEqFI(a, b)    vec_cmpeq(a, b)
#define simdCmpLtFI(a, b)    vec_cmplt(a, b)
#define simdAndFIB(a, b)     vec_and(a, b)
#define simdOrFIB(a, b)      vec_or(a, b)
#define simdAnyTrueFIB(a)          vec_any_ne(a, (__vector vsx_bool int)vec_splats(0))
#define simdMaskFI(a, sel)    vec_and(a, (__vector signed int)sel)
#define simdMaskNotFI(a, sel) vec_andc(a, (__vector signed int)sel)
#define simdBlendFI(a, b, sel)    vec_sel(a, b, sel)
/* Conversions between different booleans */
#define simdCvtFB2FIB(x)     (x)
#define simdCvtFIB2FB(x)     (x)


/****************************************************
 * SINGLE PREC. IMPLEMENTATION HELPER FUNCTIONS     *
 ****************************************************/
static inline SimdFloat
simdGetExponentF_ibm_vsx(SimdFloat x)
{
    /* Generate 0x7F800000 without memory operations */
    SimdFloat     expmask = (__vector float)vec_splats(0x7F800000);
    SimdFInt32    i127    = vec_splats(127);
    SimdFInt32    iexp;

    iexp = (__vector signed int)simdAndF(x, expmask);
    iexp = vec_sub(simdSrliFI(iexp, 23), i127);
    return vec_ctf(iexp, 0);
}

static inline SimdFloat
simdGetMantissaF_ibm_vsx(SimdFloat x)
{
    SimdFloat     expmask = (__vector float)vec_splats(0x7F800000);

    x = simdAndNotF(expmask, vec_abs(x));
    /* Reset zero (but correctly biased) exponent */
    return simdOrF(x, vec_splats(1.0f));
}

static inline SimdFloat
simdSetExponentF_ibm_vsx(SimdFloat x)
{
    SimdFInt32  iexp = simdCvtF2I(x);
    SimdFInt32  i127 = vec_splats(127);

    iexp = simdSlliFI(vec_add(iexp, i127), 23);
    return (__vector float)iexp;
}

static inline float
simdReduceF_ibm_vsx(SimdFloat x)
{
    const __vector unsigned char perm1 = { 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7 };
    const __vector unsigned char perm2 = { 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3 };

    x = vec_add(x, vec_perm(x, x, perm1));
    x = vec_add(x, vec_perm(x, x, perm2));
    return vec_extract(x, 0);
}

/* gcc-4.9 does not detect that vec_extract() uses its argument */
static inline int
simdExtractFI_ibm_vsx(SimdFInt32 gmx_unused x, unsigned int i)
{
    return vec_extract(x, i);
}

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD_FLOAT_H */
