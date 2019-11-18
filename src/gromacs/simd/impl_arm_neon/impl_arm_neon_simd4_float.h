/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2019, by the GROMACS development team, led by
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
#ifndef GMX_SIMD_IMPL_ARM_NEON_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_ARM_NEON_SIMD4_FLOAT_H

#include "config.h"

#include <cassert>
#include <cstddef>

#include <arm_neon.h>

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

    //! \brief Construct from scalar bool
    Simd4FBool(bool b) : simdInternal_(vdupq_n_u32(b ? 0xFFFFFFFF : 0)) {}

    // Internal utility constructor to simplify return statements
    Simd4FBool(uint32x4_t simd) : simdInternal_(simd) {}

    uint32x4_t simdInternal_;
};

static inline Simd4Float gmx_simdcall load4(const float* m)
{
    assert(size_t(m) % 16 == 0);
    return { vld1q_f32(m) };
}

static inline void gmx_simdcall store4(float* m, Simd4Float a)
{
    assert(size_t(m) % 16 == 0);
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
    return { vdupq_n_f32(0.0F) };
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

static inline Simd4Float gmx_simdcall operator-(Simd4Float x)
{
    return { vnegq_f32(x.simdInternal_) };
}

static inline Simd4Float gmx_simdcall operator*(Simd4Float a, Simd4Float b)
{
    return { vmulq_f32(a.simdInternal_, b.simdInternal_) };
}

// Override for Neon-Asimd
#if GMX_SIMD_ARM_NEON
static inline Simd4Float gmx_simdcall fma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
#    ifdef __ARM_FEATURE_FMA
        vfmaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)
#    else
        vmlaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)
#    endif
    };
}

static inline Simd4Float gmx_simdcall fms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
#    ifdef __ARM_FEATURE_FMA
        vnegq_f32(vfmsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_))
#    else
        vnegq_f32(vmlsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_))
#    endif
    };
}

static inline Simd4Float gmx_simdcall fnma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
#    ifdef __ARM_FEATURE_FMA
        vfmsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)
#    else
        vmlsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)
#    endif
    };
}

static inline Simd4Float gmx_simdcall fnms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
#    ifdef __ARM_FEATURE_FMA
        vnegq_f32(vfmaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_))
#    else
        vnegq_f32(vmlaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_))
#    endif
    };
}
#endif

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

// Override for Neon-Asimd
#if GMX_SIMD_ARM_NEON
static inline Simd4Float gmx_simdcall round(Simd4Float x)
{
    // Convert x to nearest integer
    float32x4_t signBitOfX = vreinterpretq_f32_u32(
            vandq_u32(vdupq_n_u32(0x80000000), vreinterpretq_u32_f32(x.simdInternal_)));
    float32x4_t half = vdupq_n_f32(0.5F);
    float32x4_t corr = vreinterpretq_f32_u32(
            vorrq_u32(vreinterpretq_u32_f32(half), vreinterpretq_u32_f32(signBitOfX)));

    int32x4_t integerX = vcvtq_s32_f32(vaddq_f32(x.simdInternal_, corr));

    // Convert back to float

    return { vcvtq_f32_s32(integerX) };
}

static inline Simd4Float gmx_simdcall trunc(Simd4Float x)
{
    return { vcvtq_f32_s32(vcvtq_s32_f32(x.simdInternal_)) };
}
#endif

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

// Override for Neon-Asimd
#if GMX_SIMD_ARM_NEON
static inline bool gmx_simdcall anyTrue(Simd4FBool a)
{
    uint32x4_t x = a.simdInternal_;
    uint32x4_t y = vextq_u32(x, x, 2);

    x = vorrq_u32(x, y);
    y = vextq_u32(x, x, 1);
    x = vorrq_u32(x, y);
    return (vgetq_lane_u32(x, 0) != 0);
}
#endif

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

// Override for Neon-Asimd
#if GMX_SIMD_ARM_NEON
static inline float gmx_simdcall reduce(Simd4Float a)
{
    float32x4_t x = a.simdInternal_;
    float32x4_t y = vextq_f32(x, x, 2);

    x = vaddq_f32(x, y);
    y = vextq_f32(x, x, 1);
    x = vaddq_f32(x, y);
    return vgetq_lane_f32(x, 0);
}

static inline float gmx_simdcall dotProduct(Simd4Float a, Simd4Float b)
{
    Simd4Float c;

    c = a * b;
    /* set 4th element to 0, then add all of them */
    c.simdInternal_ = vsetq_lane_f32(0.0F, c.simdInternal_, 3);
    return reduce(c);
}
#endif

} // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_NEON_SIMD4_FLOAT_H
