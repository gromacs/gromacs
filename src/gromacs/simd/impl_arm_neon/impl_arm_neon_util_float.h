/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <arm_neon.h>

#include "gromacs/utility/basedefinitions.h"

#include "impl_arm_neon_simd_float.h"


namespace gmx
{

template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const float *        base,
                    const std::int32_t   offset[],
                    SimdFloat *          v0,
                    SimdFloat *          v1,
                    SimdFloat *          v2,
                    SimdFloat *          v3)
{
    assert(std::size_t(offset) % 16 == 0);
    assert(std::size_t(base) % 16 == 0);
    assert(align % 4 == 0);

    // Unfortunately we cannot use the beautiful Neon structured load
    // instructions since the data comes from four different memory locations.
    float32x4x2_t  t0 = vuzpq_f32(vld1q_f32( base + align * offset[0] ), vld1q_f32( base + align * offset[2] ));
    float32x4x2_t  t1 = vuzpq_f32(vld1q_f32( base + align * offset[1] ), vld1q_f32( base + align * offset[3] ));
    float32x4x2_t  t2 = vtrnq_f32(t0.val[0], t1.val[0]);
    float32x4x2_t  t3 = vtrnq_f32(t0.val[1], t1.val[1]);
    v0->simdInternal_ = t2.val[0];
    v1->simdInternal_ = t3.val[0];
    v2->simdInternal_ = t2.val[1];
    v3->simdInternal_ = t3.val[1];
}

template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const float *        base,
                    const std::int32_t   offset[],
                    SimdFloat *          v0,
                    SimdFloat *          v1)
{
    assert(std::size_t(offset) % 16 == 0);
    assert(std::size_t(base) % 8 == 0);
    assert(align % 2 == 0);

    v0->simdInternal_  = vcombine_f32(vld1_f32( base + align * offset[0] ),
                                      vld1_f32( base + align * offset[2] ));
    v1->simdInternal_  = vcombine_f32(vld1_f32( base + align * offset[1] ),
                                      vld1_f32( base + align * offset[3] ));

    float32x4x2_t tmp  = vtrnq_f32(v0->simdInternal_, v1->simdInternal_);

    v0->simdInternal_  = tmp.val[0];
    v1->simdInternal_  = tmp.val[1];
}

static const int c_simdBestPairAlignmentFloat = 2;

template <int align>
static inline void gmx_simdcall
gatherLoadUTranspose(const float *        base,
                     const std::int32_t   offset[],
                     SimdFloat *          v0,
                     SimdFloat *          v1,
                     SimdFloat *          v2)
{
    assert(std::size_t(offset) % 16 == 0);

    float32x4x2_t  t0 = vuzpq_f32(vld1q_f32( base + align * offset[0] ), vld1q_f32( base + align * offset[2] ));
    float32x4x2_t  t1 = vuzpq_f32(vld1q_f32( base + align * offset[1] ), vld1q_f32( base + align * offset[3] ));
    float32x4x2_t  t2 = vtrnq_f32(t0.val[0], t1.val[0]);
    float32x4x2_t  t3 = vtrnq_f32(t0.val[1], t1.val[1]);
    v0->simdInternal_ = t2.val[0];
    v1->simdInternal_ = t3.val[0];
    v2->simdInternal_ = t2.val[1];
}


template <int align>
static inline void gmx_simdcall
transposeScatterStoreU(float *              base,
                       const std::int32_t   offset[],
                       SimdFloat            v0,
                       SimdFloat            v1,
                       SimdFloat            v2)
{
    assert(std::size_t(offset) % 16 == 0);

    float32x4x2_t tmp = vtrnq_f32(v0.simdInternal_, v1.simdInternal_);

    vst1_f32( base + align * offset[0], vget_low_f32(tmp.val[0]) );
    vst1_f32( base + align * offset[1], vget_low_f32(tmp.val[1]) );
    vst1_f32( base + align * offset[2], vget_high_f32(tmp.val[0]) );
    vst1_f32( base + align * offset[3], vget_high_f32(tmp.val[1]) );

    vst1q_lane_f32( base + align * offset[0] + 2, v2.simdInternal_, 0);
    vst1q_lane_f32( base + align * offset[1] + 2, v2.simdInternal_, 1);
    vst1q_lane_f32( base + align * offset[2] + 2, v2.simdInternal_, 2);
    vst1q_lane_f32( base + align * offset[3] + 2, v2.simdInternal_, 3);
}


template <int align>
static inline void gmx_simdcall
transposeScatterIncrU(float *              base,
                      const std::int32_t   offset[],
                      SimdFloat            v0,
                      SimdFloat            v1,
                      SimdFloat            v2)
{
    assert(std::size_t(offset) % 16 == 0);

    if (align < 4)
    {
        float32x2_t   t0, t1, t2, t3;
        float32x4x2_t tmp = vtrnq_f32(v0.simdInternal_, v1.simdInternal_);

        t0 = vget_low_f32(tmp.val[0]);
        t1 = vget_low_f32(tmp.val[1]);
        t2 = vget_high_f32(tmp.val[0]);
        t3 = vget_high_f32(tmp.val[1]);

        t0 = vadd_f32(t0, vld1_f32(base + align * offset[0]));
        vst1_f32(base + align * offset[0], t0);
        base[ align * offset[0] + 2] += vgetq_lane_f32(v2.simdInternal_, 0);

        t1 = vadd_f32(t1, vld1_f32(base + align * offset[1]));
        vst1_f32(base + align * offset[1], t1);
        base[ align * offset[1] + 2] += vgetq_lane_f32(v2.simdInternal_, 1);

        t2 = vadd_f32(t2, vld1_f32(base + align * offset[2]));
        vst1_f32(base + align * offset[2], t2);
        base[ align * offset[2] + 2] += vgetq_lane_f32(v2.simdInternal_, 2);

        t3 = vadd_f32(t3, vld1_f32(base + align * offset[3]));
        vst1_f32(base + align * offset[3], t3);
        base[ align * offset[3] + 2] += vgetq_lane_f32(v2.simdInternal_, 3);
    }
    else
    {
        // Extra elements means we can use full width-4 load/store operations
        float32x4x2_t  t0 = vuzpq_f32(v0.simdInternal_, v2.simdInternal_);
        float32x4x2_t  t1 = vuzpq_f32(v1.simdInternal_, vdupq_n_f32(0.0f));
        float32x4x2_t  t2 = vtrnq_f32(t0.val[0], t1.val[0]);
        float32x4x2_t  t3 = vtrnq_f32(t0.val[1], t1.val[1]);
        float32x4_t    t4 = t2.val[0];
        float32x4_t    t5 = t3.val[0];
        float32x4_t    t6 = t2.val[1];
        float32x4_t    t7 = t3.val[1];

        vst1q_f32(base + align * offset[0], vaddq_f32(t4, vld1q_f32(base + align * offset[0])));
        vst1q_f32(base + align * offset[1], vaddq_f32(t5, vld1q_f32(base + align * offset[1])));
        vst1q_f32(base + align * offset[2], vaddq_f32(t6, vld1q_f32(base + align * offset[2])));
        vst1q_f32(base + align * offset[3], vaddq_f32(t7, vld1q_f32(base + align * offset[3])));
    }
}

template <int align>
static inline void gmx_simdcall
transposeScatterDecrU(float *              base,
                      const std::int32_t   offset[],
                      SimdFloat            v0,
                      SimdFloat            v1,
                      SimdFloat            v2)
{
    assert(std::size_t(offset) % 16 == 0);

    if (align < 4)
    {
        float32x2_t   t0, t1, t2, t3;
        float32x4x2_t tmp = vtrnq_f32(v0.simdInternal_, v1.simdInternal_);

        t0 = vget_low_f32(tmp.val[0]);
        t1 = vget_low_f32(tmp.val[1]);
        t2 = vget_high_f32(tmp.val[0]);
        t3 = vget_high_f32(tmp.val[1]);

        t0 = vsub_f32(vld1_f32(base + align * offset[0]), t0);
        vst1_f32(base + align * offset[0], t0);
        base[ align * offset[0] + 2] -= vgetq_lane_f32(v2.simdInternal_, 0);

        t1 = vsub_f32(vld1_f32(base + align * offset[1]), t1);
        vst1_f32(base + align * offset[1], t1);
        base[ align * offset[1] + 2] -= vgetq_lane_f32(v2.simdInternal_, 1);

        t2 = vsub_f32(vld1_f32(base + align * offset[2]), t2);
        vst1_f32(base + align * offset[2], t2);
        base[ align * offset[2] + 2] -= vgetq_lane_f32(v2.simdInternal_, 2);

        t3 = vsub_f32(vld1_f32(base + align * offset[3]), t3);
        vst1_f32(base + align * offset[3], t3);
        base[ align * offset[3] + 2] -= vgetq_lane_f32(v2.simdInternal_, 3);
    }
    else
    {
        // Extra elements means we can use full width-4 load/store operations
        float32x4x2_t  t0 = vuzpq_f32(v0.simdInternal_, v2.simdInternal_);
        float32x4x2_t  t1 = vuzpq_f32(v1.simdInternal_, vdupq_n_f32(0.0f));
        float32x4x2_t  t2 = vtrnq_f32(t0.val[0], t1.val[0]);
        float32x4x2_t  t3 = vtrnq_f32(t0.val[1], t1.val[1]);
        float32x4_t    t4 = t2.val[0];
        float32x4_t    t5 = t3.val[0];
        float32x4_t    t6 = t2.val[1];
        float32x4_t    t7 = t3.val[1];

        vst1q_f32(base + align * offset[0], vsubq_f32(vld1q_f32(base + align * offset[0]), t4));
        vst1q_f32(base + align * offset[1], vsubq_f32(vld1q_f32(base + align * offset[1]), t5));
        vst1q_f32(base + align * offset[2], vsubq_f32(vld1q_f32(base + align * offset[2]), t6));
        vst1q_f32(base + align * offset[3], vsubq_f32(vld1q_f32(base + align * offset[3]), t7));
    }
}

static inline void gmx_simdcall
expandScalarsToTriplets(SimdFloat    scalar,
                        SimdFloat *  triplets0,
                        SimdFloat *  triplets1,
                        SimdFloat *  triplets2)
{
    float32x2_t lo, hi;
    float32x4_t t0, t1, t2, t3;

    lo = vget_low_f32(scalar.simdInternal_);
    hi = vget_high_f32(scalar.simdInternal_);

    t0 = vdupq_lane_f32(lo, 0);
    t1 = vdupq_lane_f32(lo, 1);
    t2 = vdupq_lane_f32(hi, 0);
    t3 = vdupq_lane_f32(hi, 1);

    triplets0->simdInternal_ = vextq_f32(t0, t1, 1);
    triplets1->simdInternal_ = vextq_f32(t1, t2, 2);
    triplets2->simdInternal_ = vextq_f32(t2, t3, 3);
}


template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *  base,
                             SimdFInt32     offset,
                             SimdFloat *    v0,
                             SimdFloat *    v1,
                             SimdFloat *    v2,
                             SimdFloat *    v3)
{
    GMX_ALIGNED(int, GMX_SIMD_FINT32_WIDTH)  ioffset[GMX_SIMD_FINT32_WIDTH];

    assert(std::size_t(base) % 16 == 0);
    assert(align % 4 == 0);

    store(ioffset, offset);
    gatherLoadTranspose<align>(base, ioffset, v0, v1, v2, v3);
}

template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *   base,
                             SimdFInt32      offset,
                             SimdFloat *     v0,
                             SimdFloat *     v1)
{
    GMX_ALIGNED(int, GMX_SIMD_FINT32_WIDTH)  ioffset[GMX_SIMD_FINT32_WIDTH];

    store(ioffset, offset);
    gatherLoadTranspose<align>(base, ioffset, v0, v1);
}



template <int align>
static inline void gmx_simdcall
gatherLoadUBySimdIntTranspose(const float *  base,
                              SimdFInt32     offset,
                              SimdFloat *    v0,
                              SimdFloat *    v1)
{
    GMX_ALIGNED(int, GMX_SIMD_FINT32_WIDTH)  ioffset[GMX_SIMD_FINT32_WIDTH];

    store(ioffset, offset);
    v0->simdInternal_ = vcombine_f32(vld1_f32( base + align * ioffset[0] ),
                                     vld1_f32( base + align * ioffset[2] ));
    v1->simdInternal_ = vcombine_f32(vld1_f32( base + align * ioffset[1] ),
                                     vld1_f32( base + align * ioffset[3] ));
    float32x4x2_t tmp = vtrnq_f32(v0->simdInternal_, v1->simdInternal_ );
    v0->simdInternal_ = tmp.val[0];
    v1->simdInternal_ = tmp.val[1];
}

static inline float gmx_simdcall
reduceIncr4ReturnSum(float *    m,
                     SimdFloat  v0,
                     SimdFloat  v1,
                     SimdFloat  v2,
                     SimdFloat  v3)
{
    assert(std::size_t(m) % 16 == 0);

    float32x4x2_t  t0 = vuzpq_f32(v0.simdInternal_, v2.simdInternal_);
    float32x4x2_t  t1 = vuzpq_f32(v1.simdInternal_, v3.simdInternal_);
    float32x4x2_t  t2 = vtrnq_f32(t0.val[0], t1.val[0]);
    float32x4x2_t  t3 = vtrnq_f32(t0.val[1], t1.val[1]);
    v0.simdInternal_ = t2.val[0];
    v1.simdInternal_ = t3.val[0];
    v2.simdInternal_ = t2.val[1];
    v3.simdInternal_ = t3.val[1];

    v0 = v0 + v1;
    v2 = v2 + v3;
    v0 = v0 + v2;
    v2 = v0 + simdLoad(m);
    store(m, v2);

    return reduce(v0);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_NEON_UTIL_FLOAT_H
