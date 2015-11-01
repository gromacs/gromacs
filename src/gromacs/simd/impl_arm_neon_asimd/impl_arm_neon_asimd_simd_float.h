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

#ifndef GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD_FLOAT_H

#include <math.h>

#include <arm_neon.h>

#include "impl_arm_neon_asimd_common.h"

/* NEON ASIMD always has FMA support, so make sure we use that for single too. */
#undef  simdFmaddF
#define simdFmaddF(a, b, c)  vfmaq_f32(c, b, a)
#undef  simdFmsubF
#define simdFmsubF(a, b, c)  vnegq_f32(vfmsq_f32(c, b, a))
#undef  simdFnmaddF
#define simdFnmaddF(a, b, c) vfmsq_f32(c, b, a)
#undef  simdFnmsubF
#define simdFnmsubF(a, b, c) vnegq_f32(vfmaq_f32(c, b, a))

/* The rounding instructions were actually added already in ARMv8, but most
 * compilers did not add intrinsics for them. Make sure we use them for single
 * precision too when enabling NEON Advanced SIMD.
 */
#undef  simdRoundF
#define simdRoundF(x)        vrndnq_f32(x)
#undef  simdTruncF
#define simdTruncF(x)        vrndq_f32(x)

/* NEON Advanced SIMD has a real rounding conversion instruction */
#undef  simdCvtF2I
#define simdCvtF2I(x)        vcvtnq_s32_f32(x)

/* Since we redefine rounding/conversion-with-rounding, make
 * sure we use the new operations by redefining the routine
 * to set the exponent too.
 */
#undef  simdSetExponentF
#define simdSetExponentF    simdSetExponentF_arm_neon_asimd

/* We can do more efficient reduce with vector pairwise arithmetic */
#undef  simdReduceF
#define simdReduceF(a)       simdReduceF_arm_neon_asimd(a)

/* Pick the largest unsigned integer as a shortcut for any-true */
#undef  simdAnyTrueFB
#define simdAnyTrueFB(x)     (vmaxvq_u32(x) != 0)
#undef  simdAnyTrueFIB
#define simdAnyTrueFIB(x)    (vmaxvq_u32(x) != 0)

/****************************************************
 * SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static inline SimdFloat
simdSetExponentF_arm_neon_asimd(SimdFloat x)
{
    int32x4_t  iexp = vcvtnq_s32_f32(x);

    iexp = vshlq_n_s32(vaddq_s32(iexp, vdupq_n_s32(127)), 23);
    return vreinterpretq_f32_s32(iexp);
}

static inline float
simdReduceF_arm_neon_asimd(SimdFloat a)
{
    a = vpaddq_f32(a, a);
    a = vpaddq_f32(a, a);
    return vgetq_lane_f32(a, 0);
}

#endif /* GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD_FLOAT_H */
