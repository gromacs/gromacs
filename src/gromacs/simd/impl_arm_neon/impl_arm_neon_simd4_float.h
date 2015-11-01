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

#ifndef GMX_SIMD_IMPL_ARM_NEON_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_ARM_NEON_SIMD4_FLOAT_H

#include <math.h>

#include <arm_neon.h>

#include "impl_arm_neon_common.h"
#include "impl_arm_neon_simd_float.h"

/* ARM 32-bit Neon is already 4-wide in single, so just reuse float type for SIMD4 */
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
#define simd4AbsF                  simdAbsF
#define simd4NegF                  simdNegF
#define simd4MaxF                  simdMaxF
#define simd4MinF                  simdMinF
#define simd4RoundF                simdRoundF
#define simd4TruncF                simdTruncF
#define simd4DotProductF           simd4DotProductF_arm_neon
#define Simd4FBool                 SimdFBool
#define simd4CmpEqF                simdCmpEqF
#define simd4CmpLtF                simdCmpLtF
#define simd4CmpLeF                simdCmpLeF
#define simd4AndFB                 simdAndFB
#define simd4OrFB                  simdOrFB
#define simd4AnyTrueFB             simdAnyTrueFB
#define simd4MaskF                 simdMaskF
#define simd4MaskNotF              simdMaskNotF
#define simd4BlendF                simdBlendF
#define simd4ReduceF               simdReduceF

/* SIMD4 Dot product helper function */
static inline float
simd4DotProductF_arm_neon(SimdFloat a, SimdFloat b)
{
    SimdFloat  c;
    c = simdMulF(a, b);
    /* set 4th element to 0, then add all of them */
    c = vsetq_lane_f32(0.0f, c, 3);
    return simdReduceF_arm_neon(c);
}

#endif /* GMX_SIMD_IMPL_ARM_NEON_SIMD4_FLOAT_H */
