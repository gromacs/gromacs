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

#ifndef GMX_SIMD_IMPL_ARM_NEON_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_ARM_NEON_SIMD_FLOAT_H

#include <math.h>

#include <arm_neon.h>

#include "impl_arm_neon_common.h"

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdFloat            float32x4_t
#define simdLoadF            vld1q_f32
#define simdLoad1F           vld1q_dup_f32
#define simdSet1F            vdupq_n_f32
#define simdStoreF           vst1q_f32
#define simdLoadUF           vld1q_f32
#define simdStoreUF          vst1q_f32
#define simdSetZeroF()       vdupq_n_f32(0.0f)
#define simdAddF             vaddq_f32
#define simdSubF             vsubq_f32
#define simdMulF             vmulq_f32
#ifdef __ARM_FEATURE_FMA
#    define simdFmaddF(a, b, c)  vfmaq_f32(c, b, a)
#    define simdFmsubF(a, b, c)  vnegq_f32(vfmsq_f32(c, b, a))
#    define simdFnmaddF(a, b, c) vfmaq_f32(c, b, a)
#    define simdFnmsubF(a, b, c) vnegq_f32(vfmaq_f32(c, b, a))
#else
#    define simdFmaddF(a, b, c)  vmlaq_f32(c, b, a)
#    define simdFmsubF(a, b, c)  vnegq_f32(vmlsq_f32(c, b, a))
#    define simdFnmaddF(a, b, c) vmlsq_f32(c, b, a)
#    define simdFnmsubF(a, b, c) vnegq_f32(vmlaq_f32(c, b, a))
#endif
#define simdAndF(a, b)        vreinterpretq_f32_s32(vandq_s32(vreinterpretq_s32_f32(a), vreinterpretq_s32_f32(b)))
#define simdAndNotF(a, b)     vreinterpretq_f32_s32(vbicq_s32(vreinterpretq_s32_f32(b), vreinterpretq_s32_f32(a)))
#define simdOrF(a, b)         vreinterpretq_f32_s32(vorrq_s32(vreinterpretq_s32_f32(a), vreinterpretq_s32_f32(b)))
#define simdXorF(a, b)        vreinterpretq_f32_s32(veorq_s32(vreinterpretq_s32_f32(a), vreinterpretq_s32_f32(b)))
#define simdRsqrtF            vrsqrteq_f32
#define simdRsqrtIterF(lu, x) vmulq_f32(lu, vrsqrtsq_f32(vmulq_f32(lu, lu), x))
#define simdRcpF              vrecpeq_f32
#define simdRcpIterF(lu, x)   vmulq_f32(lu, vrecpsq_f32(lu, x))
#define simdAbsF(x)         vabsq_f32(x)
#define simdNegF(x)         vnegq_f32(x)
#define simdMaxF             vmaxq_f32
#define simdMinF             vminq_f32
#define simdRoundF(x)        simdCvtI2F(simdCvtF2I(x))
#define simdTruncF(x)        simdCvtI2F(simdCvttF2I(x))
#define simdFractionF(x)     vsubq_f32(x, simdTruncF(x))
#define simdGetExponentF    simdGetExponentF_arm_neon
#define simdGetMantissaF    simdGetMantissaF_arm_neon
#define simdSetExponentF    simdSetExponentF_arm_neon
/* integer datatype corresponding to float: SimdFInt32 */
#define SimdFInt32         int32x4_t
#define simdLoadFI(m)        vld1q_s32(m)
#define simdSet1FI           vdupq_n_s32
#define simdStoreFI(m, x)    vst1q_s32(m, x)
#define simdLoadUFI(m)       vld1q_s32(m)
#define simdStoreUFI(m, x)   vst1q_s32(m, x)
#define simdSetZeroFI()      vdupq_n_s32(0)
#define simdCvttF2I          vcvtq_s32_f32
#define simdCvtF2I(x)        vcvtq_s32_f32(simdAddF(simdOrF(simdAndF(vdupq_n_f32(-0.0f), x), vdupq_n_f32(0.5f)), x))
#define simdCvtI2F           vcvtq_f32_s32
#define simdExtractFI(x, i)  vgetq_lane_s32(x, i)
/* Integer logical ops on SimdFInt32 */
#define simdSlliFI           vshlq_n_s32
#define simdSrliFI           vshrq_n_s32
#define simdAndFI            vandq_s32
#define simdAndNotFI(a, b)   vbicq_s32(b, a)
#define simdOrFI             vorrq_s32
#define simdXorFI            veorq_s32
/* Integer arithmetic ops on SimdFInt32 */
#define simdAddFI            vaddq_s32
#define simdSubFI            vsubq_s32
#define simdMulFI            vmulq_s32
/* Boolean & comparison operations on SimdFloat */
#define SimdFBool           uint32x4_t
#define simdCmpEqF           vceqq_f32
#define simdCmpLtF           vcltq_f32
#define simdCmpLeF           vcleq_f32
#define simdAndFB            vandq_u32
#define simdOrFB             vorrq_u32
#define simdAnyTrueFB        simdAnyTrueFB_arm_neon
#define simdMaskF(a, sel)     vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(a), sel))
#define simdMaskNotF(a, sel)  vreinterpretq_f32_u32(vbicq_u32(vreinterpretq_u32_f32(a), sel))
#define simdBlendF(a, b, sel)     vbslq_f32(sel, b, a)
#define simdReduceF(a)       simdReduceF_arm_neon(a)
/* Boolean & comparison operations on SimdFInt32 */
#define SimdFIBool          uint32x4_t
#define simdCmpEqFI          vceqq_s32
#define simdCmpLtFI          vcltq_s32
#define simdAndFIB           vandq_u32
#define simdOrFIB            vorrq_u32
#define simdAnyTrueFIB       simdAnyTrueFB
#define simdMaskFI(a, sel)     vandq_s32(a, vreinterpretq_s32_u32(sel))
#define simdMaskNotFI(a, sel)  vbicq_s32(a, vreinterpretq_s32_u32(sel))
#define simdBlendFI(a, b, sel)     vbslq_s32(sel, b, a)
/* Conversions between different booleans */
#define simdCvtFB2FIB(x)     (x)
#define simdCvtFIB2FB(x)     (x)

/****************************************************
 * SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static inline SimdFloat
simdGetExponentF_arm_neon(SimdFloat x)
{
    const float32x4_t expmask    = vreinterpretq_f32_s32( vdupq_n_s32(0x7F800000) );
    int32x4_t         iexp;

    iexp = vreinterpretq_s32_f32(simdAndF(x, expmask));
    iexp = vsubq_s32(vshrq_n_s32(iexp, 23), vdupq_n_s32(127));
    return vcvtq_f32_s32(iexp);
}


static inline SimdFloat
simdGetMantissaF_arm_neon(SimdFloat x)
{
    const float32x4_t mantmask   = vreinterpretq_f32_s32( vdupq_n_s32(0x007FFFFF) );
    const float32x4_t one        = vdupq_n_f32(1.0f);

    /* Get mantissa */
    x = simdAndF(mantmask, x);
    /* Reset zero (but correctly biased) exponent */
    return simdOrF(x, one);
}


static inline SimdFloat
simdSetExponentF_arm_neon(SimdFloat x)
{
    int32x4_t  iexp = simdCvtF2I(x);

    iexp = vshlq_n_s32(vaddq_s32(iexp, vdupq_n_s32(127)), 23);
    return vreinterpretq_f32_s32(iexp);
}

static inline float
simdReduceF_arm_neon(SimdFloat a)
{
    float32x4_t b = vextq_f32(a, a, 2);

    a = vaddq_f32(a, b);
    b = vextq_f32(a, a, 1);
    a = vaddq_f32(a, b);
    return vgetq_lane_f32(a, 0);
}

static inline int
simdAnyTrueFB_arm_neon(SimdFBool a)
{
    uint32x4_t b = vextq_u32(a, a, 2);

    a = simdOrFB(a, b);
    b = vextq_u32(a, a, 1);
    a = simdOrFB(a, b);
    return (vgetq_lane_u32(a, 0) != 0);
}

#endif /* GMX_SIMD_IMPL_ARM_NEON_SIMD_FLOAT_H */
