/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_ARM_NEON_ASIMD_H
#define GMX_SIMD_IMPL_ARM_NEON_ASIMD_H

#include <math.h>

#include <arm_neon.h>

/* ARM (AArch64) NEON Advanced SIMD instruction wrappers
 *
 * Please see documentation in gromacs/simd/simd.h for defines.
 */

/* Inherit single-precision and integer part from 32-bit arm */
#include "gromacs/simd/impl_arm_neon/impl_arm_neon.h"

/* Override some capability definitions from ARM 32-bit NEON - we now have double */
#define GMX_SIMD_HAVE_DOUBLE
#define GMX_SIMD_HAVE_DINT32
#define GMX_SIMD_HAVE_DINT32_EXTRACT
#define GMX_SIMD_HAVE_DINT32_LOGICAL
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS

/* Implementation details */
#define GMX_SIMD_DOUBLE_WIDTH        2
#define GMX_SIMD_DINT32_WIDTH        2

/* NEON ASIMD always has FMA support, so make sure we use that for single too. */
#undef  gmx_simd_fmadd_f
#define gmx_simd_fmadd_f(a, b, c)  vfmaq_f32(c, b, a)
#undef  gmx_simd_fmsub_f
#define gmx_simd_fmsub_f(a, b, c)  vnegq_f32(vfmsq_f32(c, b, a))
#undef  gmx_simd_fnmadd_f
#define gmx_simd_fnmadd_f(a, b, c) vfmsq_f32(c, b, a)
#undef  gmx_simd_fnmsub_f
#define gmx_simd_fnmsub_f(a, b, c) vnegq_f32(vfmaq_f32(c, b, a))

/* The rounding instructions were actually added already in ARMv8, but most
 * compilers did not add intrinsics for them. Make sure we use them for single
 * precision too when enabling NEON Advanced SIMD.
 */
#undef  gmx_simd_round_f
#define gmx_simd_round_f(x)        vrndnq_f32(x)
#undef  gmx_simd_trunc_f
#define gmx_simd_trunc_f(x)        vrndq_f32(x)

/* NEON Advanced SIMD has a real rounding conversion instruction */
#undef  gmx_simd_cvt_f2i
#define gmx_simd_cvt_f2i(x)        vcvtnq_s32_f32(x)

/* Since we redefine rounding/conversion-with-rounding, make
 * sure we use the new operations by redefining the routine
 * to set the exponent too.
 */
#undef  gmx_simd_set_exponent_f
#define gmx_simd_set_exponent_f    gmx_simd_set_exponent_f_arm_neon_asimd

/* We can do more efficient reduce with vector pairwise arithmetic */
#undef  gmx_simd_reduce_f
#define gmx_simd_reduce_f(a)       gmx_simd_reduce_f_arm_neon_asimd(a)

/* Pick the largest unsigned integer as a shortcut for any-true */
#undef  gmx_simd_anytrue_fb
#define gmx_simd_anytrue_fb(x)     (vmaxvq_u32(x) != 0)
#undef  gmx_simd_anytrue_fib
#define gmx_simd_anytrue_fib(x)    (vmaxvq_u32(x) != 0)

/* gcc-4.8 is missing the proper vreinterpretq casts
 * for 64-bit operands. However, since these datatypes
 * are opaque to the compiler we can safely cast one
 * to the other without any conversion happening.
 */

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_double_t          float64x2_t
#define gmx_simd_load_d            vld1q_f64
#define gmx_simd_load1_d           vld1q_dup_f64
#define gmx_simd_set1_d            vdupq_n_f64
#define gmx_simd_store_d           vst1q_f64
#define gmx_simd_loadu_d           vld1q_f64
#define gmx_simd_storeu_d          vst1q_f64
#define gmx_simd_setzero_d()       vdupq_n_f64(0.0)
#define gmx_simd_add_d             vaddq_f64
#define gmx_simd_sub_d             vsubq_f64
#define gmx_simd_mul_d             vmulq_f64
#define gmx_simd_fmadd_d(a, b, c)  vfmaq_f64(c, b, a)
#define gmx_simd_fmsub_d(a, b, c)  vnegq_f64(vfmsq_f64(c, b, a))
#define gmx_simd_fnmadd_d(a, b, c) vfmsq_f64(c, b, a)
#define gmx_simd_fnmsub_d(a, b, c) vnegq_f64(vfmaq_f64(c, b, a))
#define gmx_simd_and_d(a, b)        (float64x2_t)(vandq_s64((int64x2_t)(a), (int64x2_t)(b)))
#define gmx_simd_andnot_d(a, b)     (float64x2_t)(vbicq_s64((int64x2_t)(b), (int64x2_t)(a)))
#define gmx_simd_or_d(a, b)         (float64x2_t)(vorrq_s64((int64x2_t)(a), (int64x2_t)(b)))
#define gmx_simd_xor_d(a, b)        (float64x2_t)(veorq_s64((int64x2_t)(a), (int64x2_t)(b)))
#define gmx_simd_rsqrt_d            vrsqrteq_f64
#define gmx_simd_rsqrt_iter_d(lu, x) vmulq_f64(lu, vrsqrtsq_f64(vmulq_f64(lu, lu), x))
#define gmx_simd_rcp_d              vrecpeq_f64
#define gmx_simd_rcp_iter_d(lu, x)   vmulq_f64(lu, vrecpsq_f64(lu, x))
#define gmx_simd_fabs_d(x)         vabsq_f64(x)
#define gmx_simd_fneg_d(x)         vnegq_f64(x)
#define gmx_simd_max_d             vmaxq_f64
#define gmx_simd_min_d             vminq_f64
#define gmx_simd_round_d(x)        vrndnq_f64(x)
#define gmx_simd_trunc_d(x)        vrndq_f64(x)
#define gmx_simd_fraction_d(x)     vsubq_f64(x, gmx_simd_trunc_d(x))
#define gmx_simd_get_exponent_d    gmx_simd_get_exponent_d_arm_neon_asimd
#define gmx_simd_get_mantissa_d    gmx_simd_get_mantissa_d_arm_neon_asimd
#define gmx_simd_set_exponent_d    gmx_simd_set_exponent_d_arm_neon_asimd
/* integer datatype corresponding to double: gmx_simd_dint32_t */
#define gmx_simd_dint32_t          int32x2_t
#define gmx_simd_load_di(m)        vld1_s32(m)
#define gmx_simd_set1_di           vdup_n_s32
#define gmx_simd_store_di(m, x)    vst1_s32(m, x)
#define gmx_simd_loadu_di(m)       vld1_s32(m)
#define gmx_simd_storeu_di(m, x)   vst1_s32(m, x)
#define gmx_simd_setzero_di()      vdup_n_s32(0)
#define gmx_simd_cvtt_d2i(x)       vmovn_s64(vcvtq_s64_f64(x))
#define gmx_simd_cvt_d2i(x)        vmovn_s64(vcvtnq_s64_f64(x))
#define gmx_simd_cvt_i2d(x)        vcvtq_f64_s64(vmovl_s32(x))
#define gmx_simd_extract_di(x, i)  vget_lane_s32(x, i)
/* Integer logical ops on gmx_simd_dint32_t */
#define gmx_simd_slli_di           vshl_n_s32
#define gmx_simd_srli_di           vshr_n_s32
#define gmx_simd_and_di            vand_s32
#define gmx_simd_andnot_di(a, b)    vbic_s32(b, a)
#define gmx_simd_or_di             vorr_s32
#define gmx_simd_xor_di            veor_s32
/* Integer arithmetic ops on gmx_simd_dint32_t */
#define gmx_simd_add_di            vadd_s32
#define gmx_simd_sub_di            vsub_s32
#define gmx_simd_mul_di            vmul_s32
/* Boolean & comparison operations on gmx_simd_double_t */
#define gmx_simd_dbool_t           uint64x2_t
#define gmx_simd_cmpeq_d           vceqq_f64
#define gmx_simd_cmplt_d           vcltq_f64
#define gmx_simd_cmple_d           vcleq_f64
#define gmx_simd_and_db            vandq_u64
#define gmx_simd_or_db             vorrq_u64
#define gmx_simd_anytrue_db(x)     (vmaxvq_u32((uint32x4_t)(x)) != 0)
#define gmx_simd_blendzero_d(a, sel)     (float64x2_t)(vandq_u64((uint64x2_t)(a), sel))
#define gmx_simd_blendnotzero_d(a, sel)  (float64x2_t)(vbicq_u64((uint64x2_t)(a), sel))
#define gmx_simd_blendv_d(a, b, sel)     vbslq_f64(sel, b, a)
#define gmx_simd_reduce_d(a)       gmx_simd_reduce_d_arm_neon_asimd(a)
/* Boolean & comparison operations on gmx_simd_dint32_t */
#define gmx_simd_dibool_t          uint32x2_t
#define gmx_simd_cmpeq_di          vceq_s32
#define gmx_simd_cmplt_di          vclt_s32
#define gmx_simd_and_dib           vand_u32
#define gmx_simd_or_dib            vorr_u32
#define gmx_simd_anytrue_dib(x)    (vmaxv_u32(x) != 0)
#define gmx_simd_blendzero_di(a, sel)      vand_s32(a, vreinterpret_s32_u32(sel))
#define gmx_simd_blendnotzero_di(a, sel)  vbic_s32(a, vreinterpret_s32_u32(sel))
#define gmx_simd_blendv_di(a, b, sel)     vbsl_s32(sel, b, a)
/* Conversions between different booleans */
#define gmx_simd_cvt_db2dib(x)     vqmovn_u64(x)
#define gmx_simd_cvt_dib2db(x)     vorrq_u64(vmovl_u32(x), vshlq_n_u64(vmovl_u32(x), 32))

/* Float/double conversion */
#define gmx_simd_cvt_f2dd(f, d0, d1)  { *d0 = vcvt_f64_f32(vget_low_f32(f)); *d1 = vcvt_high_f64_f32(f); }
#define gmx_simd_cvt_dd2f(d0, d1)     vcvt_high_f32_f64(vcvt_f32_f64(d0), d1)

/****************************************************
 * SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static gmx_inline gmx_simd_float_t
gmx_simd_set_exponent_f_arm_neon_asimd(gmx_simd_float_t x)
{
    int32x4_t  iexp = vcvtnq_s32_f32(x);

    iexp = vshlq_n_s32(vaddq_s32(iexp, vdupq_n_s32(127)), 23);
    return vreinterpretq_f32_s32(iexp);
}

static gmx_inline float
gmx_simd_reduce_f_arm_neon_asimd(gmx_simd_float_t a)
{
    a = vpaddq_f32(a, a);
    a = vpaddq_f32(a, a);
    return vgetq_lane_f32(a, 0);
}


/****************************************************
 * DOUBLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static gmx_inline gmx_simd_double_t
gmx_simd_get_exponent_d_arm_neon_asimd(gmx_simd_double_t x)
{
    const float64x2_t expmask    = (float64x2_t)( vdupq_n_s64(0x7FF0000000000000LL) );
    int64x2_t         iexp;

    iexp = (int64x2_t)(gmx_simd_and_d(x, expmask));
    iexp = vsubq_s64(vshrq_n_s64(iexp, 52), vdupq_n_s64(1023));
    return vcvtq_f64_s64(iexp);
}


static gmx_inline gmx_simd_double_t
gmx_simd_get_mantissa_d_arm_neon_asimd(gmx_simd_double_t x)
{
    const float64x2_t mantmask   = (float64x2_t)( vdupq_n_s64(0x000FFFFFFFFFFFFFLL) );
    const float64x2_t one        = vdupq_n_f64(1.0);

    /* Get mantissa */
    x = gmx_simd_and_d(mantmask, x);
    /* Reset zero (but correctly biased) exponent */
    return gmx_simd_or_d(x, one);
}


static gmx_inline gmx_simd_double_t
gmx_simd_set_exponent_d_arm_neon_asimd(gmx_simd_double_t x)
{
    int64x2_t  iexp = vcvtnq_s64_f64(x);

    iexp = vshlq_n_s64(vaddq_s64(iexp, vdupq_n_s64(1023)), 52);
    return (float64x2_t)(iexp);
}

static gmx_inline double
gmx_simd_reduce_d_arm_neon_asimd(gmx_simd_double_t a)
{
    a = vpaddq_f64(a, a);
    return vgetq_lane_f64(a, 0);
}

#endif /* GMX_SIMD_IMPL_ARM_NEON_ASIMD_H */
