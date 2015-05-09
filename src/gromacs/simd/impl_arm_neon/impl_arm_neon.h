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

#ifndef GMX_SIMD_IMPL_ARM_NEON_H
#define GMX_SIMD_IMPL_ARM_NEON_H

#include <math.h>

#include <arm_neon.h>

/* ARM 32-bit NEON SIMD instruction wrappers
 *
 * Please see documentation in gromacs/simd/simd.h for defines.
 */

#define GMX_SIMD_V2

/* Capability definitions for ARM 32-bit NEON */
#define GMX_SIMD_HAVE_FLOAT
#undef  GMX_SIMD_HAVE_DOUBLE
#define GMX_SIMD_HAVE_HARDWARE
#define GMX_SIMD_HAVE_LOADU
#define GMX_SIMD_HAVE_STOREU
#define GMX_SIMD_HAVE_LOGICAL
#define GMX_SIMD_HAVE_FMA
#undef  GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD_HAVE_FINT32
#define GMX_SIMD_HAVE_FINT32_EXTRACT
#define GMX_SIMD_HAVE_FINT32_LOGICAL
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS
#define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT
#undef  GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT  /* No need for half-simd, width is 4 */

#undef  GMX_SIMD_HAVE_DINT32
#undef  GMX_SIMD_HAVE_DINT32_EXTRACT
#undef  GMX_SIMD_HAVE_DINT32_LOGICAL
#undef  GMX_SIMD_HAVE_DINT32_ARITHMETICS
#undef  GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_DOUBLE
#undef  GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE

#define GMX_SIMD4_HAVE_FLOAT
#undef  GMX_SIMD4_HAVE_DOUBLE

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH         4
#undef  GMX_SIMD_DOUBLE_WIDTH
#define GMX_SIMD_FINT32_WIDTH        4
#undef  GMX_SIMD_DINT32_WIDTH
#define GMX_SIMD_RSQRT_BITS          8
#define GMX_SIMD_RCP_BITS            8

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
#define gmx_simd_mul_mask_f(a, b, m)       gmx_simd_blendzero_f(gmx_simd_mul_f(a, b), m)
#define gmx_simd_fmadd_mask_f(a, b, c, m)  gmx_simd_blendzero_f(gmx_simd_fmadd_f(a, b, c), m)
#ifdef NDEBUG
#    define gmx_simd_rcp_mask_f(a, m)      gmx_simd_blendzero_f(gmx_simd_rcp_f(a), m)
#    define gmx_simd_rsqrt_mask_f(a, m)    gmx_simd_blendzero_f(gmx_simd_rsqrt_f(a), m)
#else
/* For masked rcp/rsqrt we need to make sure we do not use the masked-out arguments if FP exceptions are enabled */
#    define gmx_simd_rcp_mask_f(a, m)      gmx_simd_blendzero_f(gmx_simd_rcp_f(gmx_simd_blendv_f(gmx_simd_set1_f(1.0f), a, m)), m)
#    define gmx_simd_rsqrt_mask_f(a, m)    gmx_simd_blendzero_f(gmx_simd_rsqrt_f(gmx_simd_blendv_f(gmx_simd_set1_f(1.0f), a, m)), m)
#endif
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
#define gmx_simd_cmpnz_f(a)        vtstq_u32(vreinterpretq_u32_f32(a), vreinterpretq_u32_f32(a))
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
/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_f             gmx_simd_gather_load_transpose_f_arm_neon
#define gmx_simd_best_pair_alignment_f               gmx_simd_best_pair_alignment_f_arm_neon
#define gmx_simd_gather_loadu_transpose_f            gmx_simd_gather_loadu_transpose_f_arm_neon
#define gmx_simd_transpose_scatter_storeu_f          gmx_simd_transpose_scatter_storeu_f_arm_neon
#define gmx_simd_transpose_scatter_incru_f           gmx_simd_transpose_scatter_incru_f_arm_neon
#define gmx_simd_transpose_scatter_decru_f           gmx_simd_transpose_scatter_decru_f_arm_neon
#define gmx_simd_expand_scalars_to_triplets_f        gmx_simd_expand_scalars_to_triplets_f_arm_neon
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_arm_neon
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_arm_neon
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_arm_neon

/****************************************************
 *     NO DOUBLE PRECISION SIMD AVAILABLE           *
 ****************************************************/


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

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

#define GMX_NEON_TRANSPOSE4(v0, v1, v2, v3)                            \
    {                                                                   \
        float32x4x2_t gmx_neon_t0, gmx_neon_t1;                          \
        float32x4x2_t gmx_neon_t2, gmx_neon_t3;                          \
        gmx_neon_t0 = vuzpq_f32(v0, v2);                                \
        gmx_neon_t1 = vuzpq_f32(v1, v3);                                \
        gmx_neon_t2 = vtrnq_f32(gmx_neon_t0.val[0], gmx_neon_t1.val[0]); \
        gmx_neon_t3 = vtrnq_f32(gmx_neon_t0.val[1], gmx_neon_t1.val[1]); \
        v0          = gmx_neon_t2.val[0];                                        \
        v1          = gmx_neon_t3.val[0];                                        \
        v2          = gmx_neon_t2.val[1];                                        \
        v3          = gmx_neon_t3.val[1];                                        \
    }

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_arm_neon(const float *        base,
                                          const gmx_int32_t    offset[],
                                          gmx_simd_float_t    &v0,
                                          gmx_simd_float_t    &v1,
                                          gmx_simd_float_t    &v2,
                                          gmx_simd_float_t    &v3)
{
    /* Unfortunately we cannot use the beautiful Neon structured load
     * instructions since the data comes from four memory locations.
     */
    v0  = gmx_simd_load_f( base + align * offset[0] );
    v1  = gmx_simd_load_f( base + align * offset[1] );
    v2  = gmx_simd_load_f( base + align * offset[2] );
    v3  = gmx_simd_load_f( base + align * offset[3] );
    GMX_NEON_TRANSPOSE4(v0, v1, v2, v3);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_arm_neon(const float *        base,
                                          const gmx_int32_t    offset[],
                                          gmx_simd_float_t    &v0,
                                          gmx_simd_float_t    &v1)
{
    float32x4x2_t tmp;
    v0  = vcombine_f32(vld1_f32( base + align * offset[0] ),
                       vld1_f32( base + align * offset[2] ));
    v1  = vcombine_f32(vld1_f32( base + align * offset[1] ),
                       vld1_f32( base + align * offset[3] ));
    tmp = vtrnq_f32(v0, v1);
    v0  = tmp.val[0];
    v1  = tmp.val[1];
}

static const int gmx_simd_best_pair_alignment_f_arm_neon = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_f_arm_neon(const float *        base,
                                           const gmx_int32_t    offset[],
                                           gmx_simd_float_t    &v0,
                                           gmx_simd_float_t    &v1,
                                           gmx_simd_float_t    &v2)
{
    float32x4x2_t tmp;
    v0  = vcombine_f32(vld1_f32( base + align * offset[0] ),
                       vld1_f32( base + align * offset[2] ));
    v1  = vcombine_f32(vld1_f32( base + align * offset[1] ),
                       vld1_f32( base + align * offset[3] ));
    tmp = vtrnq_f32(v0, v1);
    v0  = tmp.val[0];
    v1  = tmp.val[1];
    v2  = vld1q_lane_f32( base + align * offset[0] + 2, v2, 0);
    v2  = vld1q_lane_f32( base + align * offset[1] + 2, v2, 1);
    v2  = vld1q_lane_f32( base + align * offset[2] + 2, v2, 2);
    v2  = vld1q_lane_f32( base + align * offset[3] + 2, v2, 3);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_f_arm_neon(float *             base,
                                             const gmx_int32_t   offset[],
                                             gmx_simd_float_t    v0,
                                             gmx_simd_float_t    v1,
                                             gmx_simd_float_t    v2)
{
    float32x4x2_t tmp;
    tmp = vtrnq_f32(v0, v1);

    vst1_f32( base + align * offset[0], vget_low_f32(tmp.val[0]) );
    vst1_f32( base + align * offset[1], vget_low_f32(tmp.val[1]) );
    vst1_f32( base + align * offset[2], vget_high_f32(tmp.val[0]) );
    vst1_f32( base + align * offset[3], vget_high_f32(tmp.val[1]) );

    vst1q_lane_f32( base + align * offset[0] + 2, v2, 0);
    vst1q_lane_f32( base + align * offset[1] + 2, v2, 1);
    vst1q_lane_f32( base + align * offset[2] + 2, v2, 2);
    vst1q_lane_f32( base + align * offset[3] + 2, v2, 3);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_f_arm_neon(float *            base,
                                            const gmx_int32_t  offset[],
                                            gmx_simd_float_t   v0,
                                            gmx_simd_float_t   v1,
                                            gmx_simd_float_t   v2)
{
    float32x4x2_t tmp;
    float32x2_t   t0, t1, t2, t3;
    tmp = vtrnq_f32(v0, v1);

    t0 = vget_low_f32(tmp.val[0]);
    t1 = vget_low_f32(tmp.val[1]);
    t2 = vget_high_f32(tmp.val[0]);
    t3 = vget_high_f32(tmp.val[1]);

    t0 = vadd_f32(t0, vld1_f32(base + align * offset[0]));
    t1 = vadd_f32(t1, vld1_f32(base + align * offset[1]));
    t2 = vadd_f32(t2, vld1_f32(base + align * offset[2]));
    t3 = vadd_f32(t3, vld1_f32(base + align * offset[3]));

    vst1_f32(base + align * offset[0], t0);
    vst1_f32(base + align * offset[1], t1);
    vst1_f32(base + align * offset[2], t2);
    vst1_f32(base + align * offset[3], t3);

    base[ align * offset[0] + 2] += vgetq_lane_f32(v2, 0);
    base[ align * offset[1] + 2] += vgetq_lane_f32(v2, 1);
    base[ align * offset[2] + 2] += vgetq_lane_f32(v2, 2);
    base[ align * offset[3] + 2] += vgetq_lane_f32(v2, 3);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_f_arm_neon(float *            base,
                                            const gmx_int32_t  offset[],
                                            gmx_simd_float_t   v0,
                                            gmx_simd_float_t   v1,
                                            gmx_simd_float_t   v2)
{
    float32x4x2_t tmp;
    float32x2_t   t0, t1, t2, t3;
    tmp = vtrnq_f32(v0, v1);

    t0 = vget_low_f32(tmp.val[0]);
    t1 = vget_low_f32(tmp.val[1]);
    t2 = vget_high_f32(tmp.val[0]);
    t3 = vget_high_f32(tmp.val[1]);

    t0 = vsub_f32(vld1_f32(base + align * offset[0]), t0);
    t1 = vsub_f32(vld1_f32(base + align * offset[1]), t1);
    t2 = vsub_f32(vld1_f32(base + align * offset[2]), t2);
    t3 = vsub_f32(vld1_f32(base + align * offset[3]), t3);

    vst1_f32(base + align * offset[0], t0);
    vst1_f32(base + align * offset[1], t1);
    vst1_f32(base + align * offset[2], t2);
    vst1_f32(base + align * offset[3], t3);

    base[ align * offset[0] + 2] -= vgetq_lane_f32(v2, 0);
    base[ align * offset[1] + 2] -= vgetq_lane_f32(v2, 1);
    base[ align * offset[2] + 2] -= vgetq_lane_f32(v2, 2);
    base[ align * offset[3] + 2] -= vgetq_lane_f32(v2, 3);
}

static gmx_inline void
gmx_simd_expand_scalars_to_triplets_f_arm_neon(gmx_simd_float_t    scalar,
                                               gmx_simd_float_t   &triplets0,
                                               gmx_simd_float_t   &triplets1,
                                               gmx_simd_float_t   &triplets2)
{
    float32x2_t lo, hi;
    float32x4_t t0, t1, t2, t3;

    lo = vget_low_f32(scalar);
    hi = vget_high_f32(scalar);

    t0 = vdupq_lane_f32(lo, 0);
    t1 = vdupq_lane_f32(lo, 1);
    t2 = vdupq_lane_f32(hi, 0);
    t3 = vdupq_lane_f32(hi, 1);

    triplets0 = vextq_f32(t0, t1, 1);
    triplets1 = vextq_f32(t1, t2, 2);
    triplets2 = vextq_f32(t2, t3, 3);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_arm_neon(const float *       base,
                                                    gmx_simd_fint32_t   offset,
                                                    gmx_simd_float_t   &v0,
                                                    gmx_simd_float_t   &v1,
                                                    gmx_simd_float_t   &v2,
                                                    gmx_simd_float_t   &v3)
{
    int ioffset[4] __attribute__((aligned(16)));
    gmx_simd_store_fi(ioffset, offset);
    gmx_simd_gather_load_transpose_f_arm_neon<align>(base, ioffset, v0, v1, v2, v3);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_arm_neon(const float *       base,
                                                    gmx_simd_fint32_t   offset,
                                                    gmx_simd_float_t   &v0,
                                                    gmx_simd_float_t   &v1)
{
    int ioffset[4] __attribute__((aligned(16)));
    gmx_simd_store_fi(ioffset, offset);
    gmx_simd_gather_load_transpose_f_arm_neon<align>(base, ioffset, v0, v1);
}


template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_f_arm_neon(const float *       base,
                                                     gmx_simd_fint32_t   offset,
                                                     gmx_simd_float_t   &v0,
                                                     gmx_simd_float_t   &v1)
{
    gmx_simd_gather_load_bysimdint_transpose_f_arm_neon<align>(base, offset, v0, v1);
}


static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_f_arm_neon(float *           m,
                                             gmx_simd_float_t  v0,
                                             gmx_simd_float_t  v1,
                                             gmx_simd_float_t  v2,
                                             gmx_simd_float_t  v3)
{
    GMX_NEON_TRANSPOSE4(v0, v1, v2, v3);
    v0 = gmx_simd_add_f(v0, v1);
    v2 = gmx_simd_add_f(v2, v3);
    v0 = gmx_simd_add_f(v0, v2);
    v2 = gmx_simd_add_f(v0, gmx_simd_load_f(m));
    gmx_simd_store_f(m, v2);

    return gmx_simd_reduce_f_arm_neon(v0);
}
#endif /* __cplusplus */


/* ARM 32-bit Neon is already 4-wide in single, so just reuse float type for SIMD4 */
#define gmx_simd4_float_t                gmx_simd_float_t
#define gmx_simd4_load_f                 gmx_simd_load_f
#define gmx_simd4_load1_f                gmx_simd_load1_f
#define gmx_simd4_set1_f                 gmx_simd_set1_f
#define gmx_simd4_store_f                gmx_simd_store_f
#define gmx_simd4_loadu_f                gmx_simd_loadu_f
#define gmx_simd4_storeu_f               gmx_simd_storeu_f
#define gmx_simd4_setzero_f              gmx_simd_setzero_f
#define gmx_simd4_add_f                  gmx_simd_add_f
#define gmx_simd4_sub_f                  gmx_simd_sub_f
#define gmx_simd4_mul_f                  gmx_simd_mul_f
#define gmx_simd4_fmadd_f                gmx_simd_fmadd_f
#define gmx_simd4_fmsub_f                gmx_simd_fmsub_f
#define gmx_simd4_fnmadd_f               gmx_simd_fnmadd_f
#define gmx_simd4_fnmsub_f               gmx_simd_fnmsub_f
#define gmx_simd4_and_f                  gmx_simd_and_f
#define gmx_simd4_andnot_f               gmx_simd_andnot_f
#define gmx_simd4_or_f                   gmx_simd_or_f
#define gmx_simd4_xor_f                  gmx_simd_xor_f
#define gmx_simd4_rsqrt_f                gmx_simd_rsqrt_f
#define gmx_simd4_fabs_f                 gmx_simd_fabs_f
#define gmx_simd4_fneg_f                 gmx_simd_fneg_f
#define gmx_simd4_max_f                  gmx_simd_max_f
#define gmx_simd4_min_f                  gmx_simd_min_f
#define gmx_simd4_round_f                gmx_simd_round_f
#define gmx_simd4_trunc_f                gmx_simd_trunc_f
#define gmx_simd4_dotproduct3_f          gmx_simd4_dotproduct3_f_arm_neon
#define gmx_simd4_transpose_f            gmx_simd4_transpose_f_arm_neon
#define gmx_simd4_fbool_t                gmx_simd_fbool_t
#define gmx_simd4_cmpeq_f                gmx_simd_cmpeq_f
#define gmx_simd4_cmplt_f                gmx_simd_cmplt_f
#define gmx_simd4_cmple_f                gmx_simd_cmple_f
#define gmx_simd4_and_fb                 gmx_simd_and_fb
#define gmx_simd4_or_fb                  gmx_simd_or_fb
#define gmx_simd4_anytrue_fb             gmx_simd_anytrue_fb
#define gmx_simd4_blendzero_f            gmx_simd_blendzero_f
#define gmx_simd4_blendnotzero_f         gmx_simd_blendnotzero_f
#define gmx_simd4_blendv_f               gmx_simd_blendv_f
#define gmx_simd4_reduce_f               gmx_simd_reduce_f

/* SIMD4 Dotproduct helper function */
static gmx_inline float
gmx_simd4_dotproduct3_f_arm_neon(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    c = gmx_simd_mul_f(a, b);
    /* set 4th element to 0, then add all of them */
    c = vsetq_lane_f32(0.0f, c, 3);
    return gmx_simd_reduce_f_arm_neon(c);
}

#ifdef __cplusplus
static gmx_inline void gmx_simdcall
gmx_simd4_transpose_f_arm_neon(gmx_simd4_float_t &v0, gmx_simd4_float_t &v1,
                               gmx_simd4_float_t &v2, gmx_simd4_float_t &v3)
{
    GMX_NEON_TRANSPOSE4(v0, v1, v2, v3);
}
#endif

/* Function to check whether SIMD operations have resulted in overflow.
 * For now, this is unfortunately a dummy for this architecture.
 */
static int
gmx_simd_check_and_reset_overflow(void)
{
    return 0;
}

#endif /* GMX_SIMD_IMPL_ARM_NEON_H */
