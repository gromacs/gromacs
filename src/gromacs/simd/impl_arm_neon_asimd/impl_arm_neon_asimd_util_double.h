/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_ARM_NEON_ASIMD_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_ARM_NEON_ASIMD_UTIL_DOUBLE_H

#include "config.h"

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <arm_neon.h>

#include "gromacs/utility/basedefinitions.h"

#include "impl_arm_neon_asimd_simd_double.h"

namespace gmx
{

template <int align>
static gmx_inline void gmx_simdcall
gatherLoadTranspose(const double *        base,
                    const std::int32_t    offset[],
                    SimdDouble *          v0,
                    SimdDouble *          v1,
                    SimdDouble *          v2,
                    SimdDouble *          v3)
{
    float64x2_t t1, t2, t3, t4;

    assert(std::size_t(offset) % 8 == 0);
    assert(std::size_t(base) % 16 == 0);
    assert(align % 2 == 0);

    t1                = vld1q_f64(base + align * offset[0]);
    t2                = vld1q_f64(base + align * offset[1]);
    t3                = vld1q_f64(base + align * offset[0] + 2);
    t4                = vld1q_f64(base + align * offset[1] + 2);
    v0->simdInternal_ = vuzp1q_f64(t1, t2);
    v1->simdInternal_ = vuzp2q_f64(t1, t2);
    v2->simdInternal_ = vuzp1q_f64(t3, t4);
    v3->simdInternal_ = vuzp2q_f64(t3, t4);
}

template <int align>
static gmx_inline void gmx_simdcall
gatherLoadTranspose(const double *        base,
                    const std::int32_t    offset[],
                    SimdDouble *          v0,
                    SimdDouble *          v1)
{
    float64x2_t t1, t2;

    assert(std::size_t(offset) % 8 == 0);
    assert(std::size_t(base) % 16 == 0);
    assert(align % 2 == 0);

    t1                = vld1q_f64(base + align * offset[0]);
    t2                = vld1q_f64(base + align * offset[1]);
    v0->simdInternal_ = vuzp1q_f64(t1, t2);
    v1->simdInternal_ = vuzp2q_f64(t1, t2);
}

static const int c_simdBestPairAlignmentDouble = 2;

template <int align>
static gmx_inline void gmx_simdcall
gatherLoadUTranspose(const double *        base,
                     const std::int32_t    offset[],
                     SimdDouble *          v0,
                     SimdDouble *          v1,
                     SimdDouble *          v2)
{
    float64x2_t t1, t2;
    float64x1_t t3, t4;

    assert(std::size_t(offset) % 8 == 0);

    t1                = vld1q_f64(base + align * offset[0]);
    t2                = vld1q_f64(base + align * offset[1]);
    t3                = vld1_f64(base + align * offset[0] + 2);
    t4                = vld1_f64(base + align * offset[1] + 2);
    v0->simdInternal_ = vuzp1q_f64(t1, t2);
    v1->simdInternal_ = vuzp2q_f64(t1, t2);
    v2->simdInternal_ = vcombine_f64(t3, t4);
}

template <int align>
static gmx_inline void gmx_simdcall
transposeScatterStoreU(double *             base,
                       const std::int32_t   offset[],
                       SimdDouble           v0,
                       SimdDouble           v1,
                       SimdDouble           v2)
{
    float64x2_t t0, t1;

    assert(std::size_t(offset) % 8 == 0);

    t0  = vuzp1q_f64(v0.simdInternal_, v1.simdInternal_);
    t1  = vuzp2q_f64(v0.simdInternal_, v1.simdInternal_);
    vst1q_f64(base + align * offset[0], t0);
    vst1q_f64(base + align * offset[1], t1);
    vst1_f64(base + align * offset[0] + 2, vget_low_f64(v2.simdInternal_));
    vst1_f64(base + align * offset[1] + 2, vget_high_f64(v2.simdInternal_));
}

template <int align>
static gmx_inline void gmx_simdcall
transposeScatterIncrU(double *             base,
                      const std::int32_t   offset[],
                      SimdDouble           v0,
                      SimdDouble           v1,
                      SimdDouble           v2)
{
    float64x2_t t0, t1, t2;
    float64x1_t t3;

    assert(std::size_t(offset) % 8 == 0);

    t0  = vuzp1q_f64(v0.simdInternal_, v1.simdInternal_); // x0 y0
    t1  = vuzp2q_f64(v0.simdInternal_, v1.simdInternal_); // x1 y1

    t2 = vld1q_f64(base + align * offset[0]);
    t2 = vaddq_f64(t2, t0);
    vst1q_f64(base + align * offset[0], t2);

    t3 = vld1_f64(base + align * offset[0] + 2);
    t3 = vadd_f64(t3, vget_low_f64(v2.simdInternal_));
    vst1_f64(base + align * offset[0] + 2, t3);

    t2 = vld1q_f64(base + align * offset[1]);
    t2 = vaddq_f64(t2, t1);
    vst1q_f64(base + align * offset[1], t2);

    t3 = vld1_f64(base + align * offset[1] + 2);
    t3 = vadd_f64(t3, vget_high_f64(v2.simdInternal_));
    vst1_f64(base + align * offset[1] + 2, t3);
}

template <int align>
static gmx_inline void gmx_simdcall
transposeScatterDecrU(double *             base,
                      const std::int32_t   offset[],
                      SimdDouble           v0,
                      SimdDouble           v1,
                      SimdDouble           v2)
{
    float64x2_t t0, t1, t2;
    float64x1_t t3;

    assert(std::size_t(offset) % 8 == 0);

    t0  = vuzp1q_f64(v0.simdInternal_, v1.simdInternal_); // x0 y0
    t1  = vuzp2q_f64(v0.simdInternal_, v1.simdInternal_); // x1 y1

    t2 = vld1q_f64(base + align * offset[0]);
    t2 = vsubq_f64(t2, t0);
    vst1q_f64(base + align * offset[0], t2);

    t3 = vld1_f64(base + align * offset[0] + 2);
    t3 = vsub_f64(t3, vget_low_f64(v2.simdInternal_));
    vst1_f64(base + align * offset[0] + 2, t3);

    t2 = vld1q_f64(base + align * offset[1]);
    t2 = vsubq_f64(t2, t1);
    vst1q_f64(base + align * offset[1], t2);

    t3 = vld1_f64(base + align * offset[1] + 2);
    t3 = vsub_f64(t3, vget_high_f64(v2.simdInternal_));
    vst1_f64(base + align * offset[1] + 2, t3);
}

static gmx_inline void gmx_simdcall
expandScalarsToTriplets(SimdDouble    scalar,
                        SimdDouble *  triplets0,
                        SimdDouble *  triplets1,
                        SimdDouble *  triplets2)
{
    triplets0->simdInternal_ = vuzp1q_f64(scalar.simdInternal_, scalar.simdInternal_);
    triplets1->simdInternal_ = scalar.simdInternal_;
    triplets2->simdInternal_ = vuzp2q_f64(scalar.simdInternal_, scalar.simdInternal_);
}

template <int align>
static gmx_inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const double *  base,
                             SimdDInt32      offset,
                             SimdDouble *    v0,
                             SimdDouble *    v1,
                             SimdDouble *    v2,
                             SimdDouble *    v3)
{
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t     ioffset[GMX_SIMD_DINT32_WIDTH];

    assert(std::size_t(base) % 16 == 0);
    assert(align % 2 == 0);

    vst1_s32(ioffset, offset.simdInternal_);
    gatherLoadTranspose<align>(base, ioffset, v0, v1, v2, v3);
}


template <int align>
static gmx_inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const double *  base,
                             SimdDInt32      offset,
                             SimdDouble *    v0,
                             SimdDouble *    v1)
{
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t     ioffset[GMX_SIMD_DINT32_WIDTH];

    assert(std::size_t(base) % 16 == 0);
    assert(align % 2 == 0);

    vst1_s32(ioffset, offset.simdInternal_);
    gatherLoadTranspose<align>(base, ioffset, v0, v1);
}

template <int align>
static gmx_inline void gmx_simdcall
gatherLoadUBySimdIntTranspose(const double *  base,
                              SimdDInt32      offset,
                              SimdDouble *    v0,
                              SimdDouble *    v1)
{
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t     ioffset[GMX_SIMD_DINT32_WIDTH];

    vst1_s32(ioffset, offset.simdInternal_);

    float64x2_t t1, t2;

    t1                = vld1q_f64(base + align * ioffset[0]);
    t2                = vld1q_f64(base + align * ioffset[1]);
    v0->simdInternal_ = vuzp1q_f64(t1, t2);
    v1->simdInternal_ = vuzp2q_f64(t1, t2);
}


static gmx_inline double gmx_simdcall
reduceIncr4ReturnSum(double *    m,
                     SimdDouble  v0,
                     SimdDouble  v1,
                     SimdDouble  v2,
                     SimdDouble  v3)
{
    float64x2_t t1, t2, t3, t4;

    assert(std::size_t(m) % 8 == 0);

    t1 = vuzp1q_f64(v0.simdInternal_, v1.simdInternal_);
    t2 = vuzp2q_f64(v0.simdInternal_, v1.simdInternal_);
    t3 = vuzp1q_f64(v2.simdInternal_, v3.simdInternal_);
    t4 = vuzp2q_f64(v2.simdInternal_, v3.simdInternal_);

    t1 = vaddq_f64(t1, t2);
    t3 = vaddq_f64(t3, t4);

    t2 = vaddq_f64(t1, vld1q_f64(m));
    t4 = vaddq_f64(t3, vld1q_f64(m + 2));
    vst1q_f64(m, t2);
    vst1q_f64(m + 2, t4);

    t1 = vaddq_f64(t1, t3);
    t2 = vpaddq_f64(t1, t1);

    return vgetq_lane_f64(t2, 0);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_NEON_ASIMD_UTIL_DOUBLE_H
