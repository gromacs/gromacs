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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD4_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD4_FLOAT_H

#include <math.h>

#include <altivec.h>

#include "impl_ibm_vsx_common.h"
#include "impl_ibm_vsx_simd_float.h"

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


/* Single precision VSX is 4 elements wide, use for SIMD4 */
#define Simd4Float                 SimdFloat
#define simd4LoadF                 simdLoadF
#define simd4Load1F                simdLoad1F
#define simd4Set1F                 simdSet1F
#define simd4StoreF                simdStoreF
#define simd4LoadUF                simdLoadUF
#define simd4StoreUF               simdStoreUF
#define simd4SetZeroF              simdSetZeroF
#define simd4AddF                  simdAddF
#define simd4SubF                  simdSubF
#define simd4MulF                  simdMulF
#define simd4FmaddF                simdFmaddF
#define simd4FmsubF                simdFmsubF
#define simd4FnmaddF               simdFnmaddF
#define simd4FnmsubF               simdFnmsubF
#define simd4AndF                  simdAndF
#define simd4AndNotF               simdAndNotF
#define simd4OrF                   simdOrF
#define simd4XorF                  simdXorF
#define simd4RsqrtF                simdRsqrtF
#define simd4AbsF                 simdAbsF
#define simd4NegF                 simdNegF
#define simd4MaxF                  simdMaxF
#define simd4MinF                  simdMinF
#define simd4RoundF                simdRoundF
#define simd4TruncF                simdTruncF
#define simd4DotProductF           simd4DotProductF_ibm_vsx
#define Simd4FInt32                SimdFInt32
#define simd4LoadFI                simdLoadFI
#define simd4Load1FI               simdLoad1FI
#define simd4Set1FI                simdSet1FI
#define simd4StoreFI               simdStoreFI
#define simd4LoadUFI               simdLoadUFI
#define simd4StoreUFI              simdStoreUFI
#define simd4SetZeroFI             simdSetZeroFI
#define Simd4FBool                SimdFBool
#define simd4CmpEqF                simdCmpEqF
#define simd4CmpLtF                simdCmpLtF
#define simd4CmpLeF                simdCmpLeF
#define simd4AndFB                 simdAndFB
#define simd4OrFB                  simdOrFB
#define simd4AnyTrueFB             simdAnyTrueFB
#define simd4MaskF            simdMaskF
#define simd4MaskNotF         simdMaskNotF
#define simd4BlendF               simdBlendF
#define simd4ReduceF               simdReduceF

static inline float
simd4DotProductF_ibm_vsx(Simd4Float a, Simd4Float b)
{
    Simd4Float                   c     = simdMulF(a, b);
    const __vector unsigned char perm1 = { 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7 };
    const __vector unsigned char perm2 = { 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3 };
    __vector float               sum;
    sum = vec_add(c, vec_perm(c, c, (__vector unsigned char)perm1));
    sum = vec_add(sum, vec_perm(c, c, (__vector unsigned char)perm2));
    return vec_extract(sum, 0);
}

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD4_FLOAT_H */
