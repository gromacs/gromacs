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

#ifndef GMX_SIMD_IMPL_ARM_NEON_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_ARM_NEON_UTIL_FLOAT_H

#include "config.h"

#include <assert.h>

#include <arm_neon.h>

#include "impl_arm_neon_common.h"

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
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

#define GMX_NEON_TRANSPOSE4(v0, v1, v2, v3)                          \
    {                                                                    \
        float32x4x2_t gmx_neon_t0, gmx_neon_t1;                          \
        float32x4x2_t gmx_neon_t2, gmx_neon_t3;                          \
        gmx_neon_t0 = vuzpq_f32(v0, v2);                                 \
        gmx_neon_t1 = vuzpq_f32(v1, v3);                                 \
        gmx_neon_t2 = vtrnq_f32(gmx_neon_t0.val[0], gmx_neon_t1.val[0]); \
        gmx_neon_t3 = vtrnq_f32(gmx_neon_t0.val[1], gmx_neon_t1.val[1]); \
        v0          = gmx_neon_t2.val[0];                                \
        v1          = gmx_neon_t3.val[0];                                \
        v2          = gmx_neon_t2.val[1];                                \
        v3          = gmx_neon_t3.val[1];                                \
    }

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_transpose_f_arm_neon(const float *        base,
                                          const gmx_int32_t    offset[],
                                          gmx_simd_float_t    &v0,
                                          gmx_simd_float_t    &v1,
                                          gmx_simd_float_t    &v2,
                                          gmx_simd_float_t    &v3)
{
    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);

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
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_transpose_f_arm_neon(const float *        base,
                                          const gmx_int32_t    offset[],
                                          gmx_simd_float_t    &v0,
                                          gmx_simd_float_t    &v1)
{
    float32x4x2_t tmp;

    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 8 == 0);
    assert(align % 2 == 0);

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
static gmx_inline void gmx_simdcall
gmx_simd_gather_loadu_transpose_f_arm_neon(const float *        base,
                                           const gmx_int32_t    offset[],
                                           gmx_simd_float_t    &v0,
                                           gmx_simd_float_t    &v1,
                                           gmx_simd_float_t    &v2)
{
    float32x4x2_t tmp;

    assert((size_t)offset % 16 == 0);

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
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_storeu_f_arm_neon(float *             base,
                                             const gmx_int32_t   offset[],
                                             gmx_simd_float_t    v0,
                                             gmx_simd_float_t    v1,
                                             gmx_simd_float_t    v2)
{
    float32x4x2_t tmp;

    assert((size_t)offset % 16 == 0);

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
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_incru_f_arm_neon(float *            base,
                                            const gmx_int32_t  offset[],
                                            gmx_simd_float_t   v0,
                                            gmx_simd_float_t   v1,
                                            gmx_simd_float_t   v2)
{
    float32x4x2_t tmp;
    float32x2_t   t0, t1, t2, t3;

    assert((size_t)offset % 16 == 0);

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
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_decru_f_arm_neon(float *            base,
                                            const gmx_int32_t  offset[],
                                            gmx_simd_float_t   v0,
                                            gmx_simd_float_t   v1,
                                            gmx_simd_float_t   v2)
{
    float32x4x2_t tmp;
    float32x2_t   t0, t1, t2, t3;

    assert((size_t)offset % 16 == 0);

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

static gmx_inline void gmx_simdcall
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
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_f_arm_neon(const float *       base,
                                                    gmx_simd_fint32_t   offset,
                                                    gmx_simd_float_t   &v0,
                                                    gmx_simd_float_t   &v1,
                                                    gmx_simd_float_t   &v2,
                                                    gmx_simd_float_t   &v3)
{
    int ioffset[4] __attribute__((aligned(16)));

    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);

    gmx_simd_store_fi(ioffset, offset);
    gmx_simd_gather_load_transpose_f_arm_neon<align>(base, ioffset, v0, v1, v2, v3);
}


template <int align>
static gmx_inline void gmx_simdcall
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
static gmx_inline void gmx_simdcall
gmx_simd_gather_loadu_bysimdint_transpose_f_arm_neon(const float  *       base,
                                                     gmx_simd_fint32_t    offset,
                                                     gmx_simd_float_t    &v0,
                                                     gmx_simd_float_t    &v1)
{
    float32x4x2_t tmp;
    int           ioffset[4] __attribute__((aligned(16)));

    gmx_simd_store_fi(ioffset, offset);
    v0  = vcombine_f32(vld1_f32( base + align * ioffset[0] ),
                       vld1_f32( base + align * ioffset[2] ));
    v1  = vcombine_f32(vld1_f32( base + align * ioffset[1] ),
                       vld1_f32( base + align * ioffset[3] ));
    tmp = vtrnq_f32(v0, v1);
    v0  = tmp.val[0];
    v1  = tmp.val[1];
}


static gmx_inline float gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_f_arm_neon(float *           m,
                                             gmx_simd_float_t  v0,
                                             gmx_simd_float_t  v1,
                                             gmx_simd_float_t  v2,
                                             gmx_simd_float_t  v3)
{
    assert((size_t)m % 16 == 0);

    GMX_NEON_TRANSPOSE4(v0, v1, v2, v3);
    v0 = gmx_simd_add_f(v0, v1);
    v2 = gmx_simd_add_f(v2, v3);
    v0 = gmx_simd_add_f(v0, v2);
    v2 = gmx_simd_add_f(v0, gmx_simd_load_f(m));
    gmx_simd_store_f(m, v2);

    return gmx_simd_reduce_f_arm_neon(v0);
}
#endif /* __cplusplus */


#endif /* GMX_SIMD_IMPL_ARM_NEON_UTIL_FLOAT_H */
