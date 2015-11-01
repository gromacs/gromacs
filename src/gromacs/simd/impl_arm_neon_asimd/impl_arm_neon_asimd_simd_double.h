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

#ifndef GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD_DOUBLE_H

#include <math.h>

#include <arm_neon.h>

#include "impl_arm_neon_asimd_common.h"

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdDouble          float64x2_t
#define simdLoadD            vld1q_f64
#define simdLoad1D           vld1q_dup_f64
#define simdSet1D            vdupq_n_f64
#define simdStoreD           vst1q_f64
#define simdLoadUD           vld1q_f64
#define simdStoreUD          vst1q_f64
#define simdSetZeroD()       vdupq_n_f64(0.0)
#define simdAddD             vaddq_f64
#define simdSubD             vsubq_f64
#define simdMulD             vmulq_f64
#define simdFmaddD(a, b, c)  vfmaq_f64(c, b, a)
#define simdFmsubD(a, b, c)  vnegq_f64(vfmsq_f64(c, b, a))
#define simdFnmaddD(a, b, c) vfmsq_f64(c, b, a)
#define simdFnmsubD(a, b, c) vnegq_f64(vfmaq_f64(c, b, a))
#define simdAndD(a, b)        (float64x2_t)(vandq_s64((int64x2_t)(a), (int64x2_t)(b)))
#define simdAndNotD(a, b)     (float64x2_t)(vbicq_s64((int64x2_t)(b), (int64x2_t)(a)))
#define simdOrD(a, b)         (float64x2_t)(vorrq_s64((int64x2_t)(a), (int64x2_t)(b)))
#define simdXorD(a, b)        (float64x2_t)(veorq_s64((int64x2_t)(a), (int64x2_t)(b)))
#define simdRsqrtD            vrsqrteq_f64
#define simdRsqrtIterD(lu, x) vmulq_f64(lu, vrsqrtsq_f64(vmulq_f64(lu, lu), x))
#define simdRcpD              vrecpeq_f64
#define simdRcpIterD(lu, x)   vmulq_f64(lu, vrecpsq_f64(lu, x))
#define simdAbsD(x)         vabsq_f64(x)
#define simdNegD(x)         vnegq_f64(x)
#define simdMaxD             vmaxq_f64
#define simdMinD             vminq_f64
#define simdRoundD(x)        vrndnq_f64(x)
#define simdTruncD(x)        vrndq_f64(x)
#define simdFractionD(x)     vsubq_f64(x, simdTruncD(x))
#define simdGetExponentD    simdGetExponentD_arm_neon_asimd
#define simdGetMantissaD    simdGetMantissaD_arm_neon_asimd
#define simdSetExponentD    simdSetExponentD_arm_neon_asimd
/* integer datatype corresponding to double: SimdDInt32 */
#define SimdDInt32          int32x2_t
#define simdLoadDI(m)        vld1_s32(m)
#define simdSet1DI           vdup_n_s32
#define simdStoreDI(m, x)    vst1_s32(m, x)
#define simdLoadUDI(m)       vld1_s32(m)
#define simdStoreUDI(m, x)   vst1_s32(m, x)
#define simdSetZeroDI()      vdup_n_s32(0)
#define simdCvttD2I(x)       vmovn_s64(vcvtq_s64_f64(x))
#define simdCvtD2I(x)        vmovn_s64(vcvtnq_s64_f64(x))
#define simdCvtI2D(x)        vcvtq_f64_s64(vmovl_s32(x))
#define simdExtractDI(x, i)  vget_lane_s32(x, i)
/* Integer logical ops on SimdDInt32 */
#define simdSlliDI           vshl_n_s32
#define simdSrliDI           vshr_n_s32
#define simdAndDI            vand_s32
#define simdAndNotDI(a, b)    vbic_s32(b, a)
#define simdOrDI             vorr_s32
#define simdXorDI            veor_s32
/* Integer arithmetic ops on SimdDInt32 */
#define simdAddDI            vadd_s32
#define simdSubDI            vsub_s32
#define simdMulDI            vmul_s32
/* Boolean & comparison operations on SimdDouble */
#define SimdDBool           uint64x2_t
#define simdCmpEqD           vceqq_f64
#define simdCmpLtD           vcltq_f64
#define simdCmpLeD           vcleq_f64
#define simdAndDB            vandq_u64
#define simdOrDB             vorrq_u64
#define simdAnyTrueDB(x)     (vmaxvq_u32((uint32x4_t)(x)) != 0)
#define simdMaskD(a, sel)     (float64x2_t)(vandq_u64((uint64x2_t)(a), sel))
#define simdMaskNotD(a, sel)  (float64x2_t)(vbicq_u64((uint64x2_t)(a), sel))
#define simdBlendD(a, b, sel)     vbslq_f64(sel, b, a)
#define simdReduceD(a)       simdReduceD_arm_neon_asimd(a)
/* Boolean & comparison operations on SimdDInt32 */
#define SimdDIBool          uint32x2_t
#define simdCmpEqDI          vceq_s32
#define simdCmpLtDI          vclt_s32
#define simdAndDIB           vand_u32
#define simdOrDIB            vorr_u32
#define simdAnyTrueDIB(x)    (vmaxv_u32(x) != 0)
#define simdMaskDI(a, sel)      vand_s32(a, vreinterpret_s32_u32(sel))
#define simdMaskNotDI(a, sel)  vbic_s32(a, vreinterpret_s32_u32(sel))
#define simdBlendDI(a, b, sel)     vbsl_s32(sel, b, a)
/* Conversions between different booleans */
#define simdCvtDB2DIB(x)     vqmovn_u64(x)
#define simdCvtDIB2DB(x)     vorrq_u64(vmovl_u32(x), vshlq_n_u64(vmovl_u32(x), 32))

/* Float/double conversion */
#define simdCvtF2DD(f, d0, d1)  { *d0 = vcvt_f64_f32(vget_low_f32(f)); *d1 = vcvt_high_f64_f32(f); }
#define simdCvtDD2F(d0, d1)     vcvt_high_f32_f64(vcvt_f32_f64(d0), d1)

/****************************************************
 * DOUBLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static inline SimdDouble
simdGetExponentD_arm_neon_asimd(SimdDouble x)
{
    const float64x2_t expmask    = (float64x2_t)( vdupq_n_s64(0x7FF0000000000000LL) );
    int64x2_t         iexp;

    iexp = (int64x2_t)(simdAndD(x, expmask));
    iexp = vsubq_s64(vshrq_n_s64(iexp, 52), vdupq_n_s64(1023));
    return vcvtq_f64_s64(iexp);
}


static inline SimdDouble
simdGetMantissaD_arm_neon_asimd(SimdDouble x)
{
    const float64x2_t mantmask   = (float64x2_t)( vdupq_n_s64(0x000FFFFFFFFFFFFFLL) );
    const float64x2_t one        = vdupq_n_f64(1.0);

    /* Get mantissa */
    x = simdAndD(mantmask, x);
    /* Reset zero (but correctly biased) exponent */
    return simdOrD(x, one);
}


static inline SimdDouble
simdSetExponentD_arm_neon_asimd(SimdDouble x)
{
    int64x2_t  iexp = vcvtnq_s64_f64(x);

    iexp = vshlq_n_s64(vaddq_s64(iexp, vdupq_n_s64(1023)), 52);
    return (float64x2_t)(iexp);
}

static inline double
simdReduceD_arm_neon_asimd(SimdDouble a)
{
    a = vpaddq_f64(a, a);
    return vgetq_lane_f64(a, 0);
}

#endif /* GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD_DOUBLE_H */
