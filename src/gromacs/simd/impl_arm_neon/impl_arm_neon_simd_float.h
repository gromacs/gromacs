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
#define gmx_simd_float_t           float32x4_t
#define gmx_simd_load_f            vld1q_f32
#define gmx_simd_load1_f           vld1q_dup_f32
#define gmx_simd_set1_f            vdupq_n_f32
#define gmx_simd_store_f           vst1q_f32
#define gmx_simd_loadu_f           vld1q_f32
#define gmx_simd_storeu_f          vst1q_f32
#define gmx_simd_setzero_f()       vdupq_n_f32(0.0f)
#define gmx_simd_add_f             vaddq_f32
#define gmx_simd_sub_f             vsubq_f32
#define gmx_simd_mul_f             vmulq_f32
#ifdef __ARM_FEATURE_FMA
#    define gmx_simd_fmadd_f(a, b, c)  vfmaq_f32(c, b, a)
#    define gmx_simd_fmsub_f(a, b, c)  vnegq_f32(vfmsq_f32(c, b, a))
#    define gmx_simd_fnmadd_f(a, b, c) vfmaq_f32(c, b, a)
#    define gmx_simd_fnmsub_f(a, b, c) vnegq_f32(vfmaq_f32(c, b, a))
#else
#    define gmx_simd_fmadd_f(a, b, c)  vmlaq_f32(c, b, a)
#    define gmx_simd_fmsub_f(a, b, c)  vnegq_f32(vmlsq_f32(c, b, a))
#    define gmx_simd_fnmadd_f(a, b, c) vmlsq_f32(c, b, a)
#    define gmx_simd_fnmsub_f(a, b, c) vnegq_f32(vmlaq_f32(c, b, a))
#endif
#define gmx_simd_and_f(a, b)        vreinterpretq_f32_s32(vandq_s32(vreinterpretq_s32_f32(a), vreinterpretq_s32_f32(b)))
#define gmx_simd_andnot_f(a, b)     vreinterpretq_f32_s32(vbicq_s32(vreinterpretq_s32_f32(b), vreinterpretq_s32_f32(a)))
#define gmx_simd_or_f(a, b)         vreinterpretq_f32_s32(vorrq_s32(vreinterpretq_s32_f32(a), vreinterpretq_s32_f32(b)))
#define gmx_simd_xor_f(a, b)        vreinterpretq_f32_s32(veorq_s32(vreinterpretq_s32_f32(a), vreinterpretq_s32_f32(b)))
#define gmx_simd_rsqrt_f            vrsqrteq_f32
#define gmx_simd_rsqrt_iter_f(lu, x) vmulq_f32(lu, vrsqrtsq_f32(vmulq_f32(lu, lu), x))
#define gmx_simd_rcp_f              vrecpeq_f32
#define gmx_simd_rcp_iter_f(lu, x)   vmulq_f32(lu, vrecpsq_f32(lu, x))
#define gmx_simd_fabs_f(x)         vabsq_f32(x)
#define gmx_simd_fneg_f(x)         vnegq_f32(x)
#define gmx_simd_max_f             vmaxq_f32
#define gmx_simd_min_f             vminq_f32
#define gmx_simd_round_f(x)        gmx_simd_cvt_i2f(gmx_simd_cvt_f2i(x))
#define gmx_simd_trunc_f(x)        gmx_simd_cvt_i2f(gmx_simd_cvtt_f2i(x))
#define gmx_simd_fraction_f(x)     vsubq_f32(x, gmx_simd_trunc_f(x))
#define gmx_simd_get_exponent_f    gmx_simd_get_exponent_f_arm_neon
#define gmx_simd_get_mantissa_f    gmx_simd_get_mantissa_f_arm_neon
#define gmx_simd_set_exponent_f    gmx_simd_set_exponent_f_arm_neon
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_fint32_t         int32x4_t
#define gmx_simd_load_fi(m)        vld1q_s32(m)
#define gmx_simd_set1_fi           vdupq_n_s32
#define gmx_simd_store_fi(m, x)    vst1q_s32(m, x)
#define gmx_simd_loadu_fi(m)       vld1q_s32(m)
#define gmx_simd_storeu_fi(m, x)   vst1q_s32(m, x)
#define gmx_simd_setzero_fi()      vdupq_n_s32(0)
#define gmx_simd_cvtt_f2i          vcvtq_s32_f32
#define gmx_simd_cvt_f2i(x)        vcvtq_s32_f32(gmx_simd_add_f(gmx_simd_or_f(gmx_simd_and_f(vdupq_n_f32(-0.0f), x), vdupq_n_f32(0.5f)), x))
#define gmx_simd_cvt_i2f           vcvtq_f32_s32
#define gmx_simd_extract_fi(x, i)  vgetq_lane_s32(x, i)
/* Integer logical ops on gmx_simd_fint32_t */
#define gmx_simd_slli_fi           vshlq_n_s32
#define gmx_simd_srli_fi           vshrq_n_s32
#define gmx_simd_and_fi            vandq_s32
#define gmx_simd_andnot_fi(a, b)   vbicq_s32(b, a)
#define gmx_simd_or_fi             vorrq_s32
#define gmx_simd_xor_fi            veorq_s32
/* Integer arithmetic ops on gmx_simd_fint32_t */
#define gmx_simd_add_fi            vaddq_s32
#define gmx_simd_sub_fi            vsubq_s32
#define gmx_simd_mul_fi            vmulq_s32
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_fbool_t           uint32x4_t
#define gmx_simd_cmpeq_f           vceqq_f32
#define gmx_simd_cmplt_f           vcltq_f32
#define gmx_simd_cmple_f           vcleq_f32
#define gmx_simd_and_fb            vandq_u32
#define gmx_simd_or_fb             vorrq_u32
#define gmx_simd_anytrue_fb        gmx_simd_anytrue_fb_arm_neon
#define gmx_simd_blendzero_f(a, sel)     vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(a), sel))
#define gmx_simd_blendnotzero_f(a, sel)  vreinterpretq_f32_u32(vbicq_u32(vreinterpretq_u32_f32(a), sel))
#define gmx_simd_blendv_f(a, b, sel)     vbslq_f32(sel, b, a)
#define gmx_simd_reduce_f(a)       gmx_simd_reduce_f_arm_neon(a)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_fibool_t          uint32x4_t
#define gmx_simd_cmpeq_fi          vceqq_s32
#define gmx_simd_cmplt_fi          vcltq_s32
#define gmx_simd_and_fib           vandq_u32
#define gmx_simd_or_fib            vorrq_u32
#define gmx_simd_anytrue_fib       gmx_simd_anytrue_fb
#define gmx_simd_blendzero_fi(a, sel)     vandq_s32(a, vreinterpretq_s32_u32(sel))
#define gmx_simd_blendnotzero_fi(a, sel)  vbicq_s32(a, vreinterpretq_s32_u32(sel))
#define gmx_simd_blendv_fi(a, b, sel)     vbslq_s32(sel, b, a)
/* Conversions between different booleans */
#define gmx_simd_cvt_fb2fib(x)     (x)
#define gmx_simd_cvt_fib2fb(x)     (x)

/****************************************************
 * SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static gmx_inline gmx_simd_float_t
gmx_simd_get_exponent_f_arm_neon(gmx_simd_float_t x)
{
    const float32x4_t expmask    = vreinterpretq_f32_s32( vdupq_n_s32(0x7F800000) );
    int32x4_t         iexp;

    iexp = vreinterpretq_s32_f32(gmx_simd_and_f(x, expmask));
    iexp = vsubq_s32(vshrq_n_s32(iexp, 23), vdupq_n_s32(127));
    return vcvtq_f32_s32(iexp);
}


static gmx_inline gmx_simd_float_t
gmx_simd_get_mantissa_f_arm_neon(gmx_simd_float_t x)
{
    const float32x4_t mantmask   = vreinterpretq_f32_s32( vdupq_n_s32(0x007FFFFF) );
    const float32x4_t one        = vdupq_n_f32(1.0f);

    /* Get mantissa */
    x = gmx_simd_and_f(mantmask, x);
    /* Reset zero (but correctly biased) exponent */
    return gmx_simd_or_f(x, one);
}


static gmx_inline gmx_simd_float_t
gmx_simd_set_exponent_f_arm_neon(gmx_simd_float_t x)
{
    int32x4_t  iexp = gmx_simd_cvt_f2i(x);

    iexp = vshlq_n_s32(vaddq_s32(iexp, vdupq_n_s32(127)), 23);
    return vreinterpretq_f32_s32(iexp);
}

static gmx_inline float
gmx_simd_reduce_f_arm_neon(gmx_simd_float_t a)
{
    float32x4_t b = vextq_f32(a, a, 2);

    a = vaddq_f32(a, b);
    b = vextq_f32(a, a, 1);
    a = vaddq_f32(a, b);
    return vgetq_lane_f32(a, 0);
}

static gmx_inline int
gmx_simd_anytrue_fb_arm_neon(gmx_simd_fbool_t a)
{
    uint32x4_t b = vextq_u32(a, a, 2);

    a = gmx_simd_or_fb(a, b);
    b = vextq_u32(a, a, 1);
    a = gmx_simd_or_fb(a, b);
    return (vgetq_lane_u32(a, 0) != 0);
}

#endif /* GMX_SIMD_IMPL_ARM_NEON_SIMD_FLOAT_H */
