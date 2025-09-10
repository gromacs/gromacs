/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*
 * armv8+sve support to GROMACS was contributed by the Research Organization for
 * Information Science and Technology (RIST).
 * Copyright (c) 2020 Research Organization for Information Science and Technology (RIST).
 */

#ifndef GMX_SIMD_IMPL_ARM_SVE_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_ARM_SVE_SIMD4_FLOAT_H

#include "config.h"

#include <arm_neon.h>

#include <cassert>
#include <cstddef>
#include <cstdint>

#include "gromacs/math/utilities.h"

namespace gmx
{

class Simd4Float
{
public:
    Simd4Float() {}

    Simd4Float(float f) : simdInternal_(vdupq_n_f32(f)) {}

    // Internal utility constructor to simplify return statements
    Simd4Float(float32x4_t simd) : simdInternal_(simd) {}

    float32x4_t simdInternal_;
};

class Simd4FBool
{
public:
    Simd4FBool() {}

    Simd4FBool(bool b) : simdInternal_(vdupq_n_u32(b ? 0xFFFFFFFF : 0)) {}

    // Internal utility constructor to simplify return statements
    Simd4FBool(uint32x4_t simd) : simdInternal_(simd) {}

    uint32x4_t simdInternal_;
};

static inline Simd4Float gmx_simdcall load4(const float* m)
{
    assert(std::size_t(m) % (GMX_SIMD4_WIDTH * sizeof(float)) == 0);
    return { vld1q_f32(m) };
}

static inline void gmx_simdcall store4(float* m, Simd4Float a)
{
    assert(std::size_t(m) % (GMX_SIMD4_WIDTH * sizeof(float)) == 0);
    vst1q_f32(m, a.simdInternal_);
}

static inline Simd4Float gmx_simdcall load4U(const float* m)
{
    return { vld1q_f32(m) };
}

static inline void gmx_simdcall store4U(float* m, Simd4Float a)
{
    vst1q_f32(m, a.simdInternal_);
}

static inline Simd4Float gmx_simdcall simd4SetZeroF()
{
    return { vdupq_n_f32(0.0f) };
}

static inline Simd4Float gmx_simdcall operator&(Simd4Float a, Simd4Float b)
{
    return { vreinterpretq_f32_s32(vandq_s32(vreinterpretq_s32_f32(a.simdInternal_),
                                             vreinterpretq_s32_f32(b.simdInternal_))) };
}

static inline Simd4Float gmx_simdcall andNot(Simd4Float a, Simd4Float b)
{
    return { vreinterpretq_f32_s32(vbicq_s32(vreinterpretq_s32_f32(b.simdInternal_),
                                             vreinterpretq_s32_f32(a.simdInternal_))) };
}

static inline Simd4Float gmx_simdcall operator|(Simd4Float a, Simd4Float b)
{
    return { vreinterpretq_f32_s32(vorrq_s32(vreinterpretq_s32_f32(a.simdInternal_),
                                             vreinterpretq_s32_f32(b.simdInternal_))) };
}

static inline Simd4Float gmx_simdcall operator^(Simd4Float a, Simd4Float b)
{
    return { vreinterpretq_f32_s32(veorq_s32(vreinterpretq_s32_f32(a.simdInternal_),
                                             vreinterpretq_s32_f32(b.simdInternal_))) };
}

static inline Simd4Float gmx_simdcall operator+(Simd4Float a, Simd4Float b)
{
    return { vaddq_f32(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall operator-(Simd4Float a, Simd4Float b)
{
    return { vsubq_f32(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall operator-(Simd4Float a)
{
    return { vnegq_f32(a.simdInternal_) };
}

static inline Simd4Float gmx_simdcall operator*(Simd4Float a, Simd4Float b)
{
    return { vmulq_f32(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { vfmaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { vnegq_f32(vfmsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)) };
}

static inline Simd4Float gmx_simdcall fnma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { vfmsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fnms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { vnegq_f32(vfmaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)) };
}

static inline Simd4Float gmx_simdcall rsqrt(Simd4Float x)
{
    return { vrsqrteq_f32(x.simdInternal_) };
}

static inline Simd4Float gmx_simdcall abs(Simd4Float x)
{
    return { vabsq_f32(x.simdInternal_) };
}

static inline Simd4Float gmx_simdcall max(Simd4Float a, Simd4Float b)
{
    return { vmaxq_f32(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall min(Simd4Float a, Simd4Float b)
{
    return { vminq_f32(a.simdInternal_, b.simdInternal_) };
}

// Round and trunc operations are defined at the end of this file, since they
// need to use float-to-integer and integer-to-float conversions.

static inline float gmx_simdcall reduce(Simd4Float a)
{
    float32x4_t b = a.simdInternal_;
    b             = vpaddq_f32(b, b);
    b             = vpaddq_f32(b, b);
    return vgetq_lane_f32(b, 0);
}

static inline Simd4FBool gmx_simdcall operator==(Simd4Float a, Simd4Float b)
{
    return { vceqq_f32(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4FBool gmx_simdcall operator!=(Simd4Float a, Simd4Float b)
{
    return { vmvnq_u32(vceqq_f32(a.simdInternal_, b.simdInternal_)) };
}

static inline Simd4FBool gmx_simdcall operator<(Simd4Float a, Simd4Float b)
{
    return { vcltq_f32(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4FBool gmx_simdcall operator<=(Simd4Float a, Simd4Float b)
{
    return { vcleq_f32(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4FBool gmx_simdcall operator&&(Simd4FBool a, Simd4FBool b)
{
    return { vandq_u32(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4FBool gmx_simdcall operator||(Simd4FBool a, Simd4FBool b)
{
    return { vorrq_u32(a.simdInternal_, b.simdInternal_) };
}

static inline bool gmx_simdcall anyTrue(Simd4FBool a)
{
    return (vmaxvq_u32(a.simdInternal_) != 0);
}

static inline Simd4Float gmx_simdcall selectByMask(Simd4Float a, Simd4FBool m)
{
    return { vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(a.simdInternal_), m.simdInternal_)) };
}

static inline Simd4Float gmx_simdcall selectByNotMask(Simd4Float a, Simd4FBool m)
{
    return { vreinterpretq_f32_u32(vbicq_u32(vreinterpretq_u32_f32(a.simdInternal_), m.simdInternal_)) };
}

static inline Simd4Float gmx_simdcall blend(Simd4Float a, Simd4Float b, Simd4FBool sel)
{
    return { vbslq_f32(sel.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline Simd4Float gmx_simdcall round(Simd4Float x)
{
    return { vrndnq_f32(x.simdInternal_) };
}

static inline Simd4Float gmx_simdcall trunc(Simd4Float x)
{
    return { vrndq_f32(x.simdInternal_) };
}

static inline float gmx_simdcall dotProduct(Simd4Float a, Simd4Float b)
{
    Simd4Float c;

    c = a * b;
    /* set 4th element to 0, then add all of them */
    c.simdInternal_ = vsetq_lane_f32(0.0f, c.simdInternal_, 3);
    return reduce(c);
}

static inline void gmx_simdcall transpose(Simd4Float* v0, Simd4Float* v1, Simd4Float* v2, Simd4Float* v3)
{
    float32x4x2_t t0  = vuzpq_f32(v0->simdInternal_, v2->simdInternal_);
    float32x4x2_t t1  = vuzpq_f32(v1->simdInternal_, v3->simdInternal_);
    float32x4x2_t t2  = vtrnq_f32(t0.val[0], t1.val[0]);
    float32x4x2_t t3  = vtrnq_f32(t0.val[1], t1.val[1]);
    v0->simdInternal_ = t2.val[0];
    v1->simdInternal_ = t3.val[0];
    v2->simdInternal_ = t2.val[1];
    v3->simdInternal_ = t3.val[1];
}

} // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_SVE_SIMD4_FLOAT_H
