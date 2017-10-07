/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "config.h"

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <arm_neon.h>

#include "gromacs/math/utilities.h"

namespace gmx
{

class SimdFloat
{
    public:
        SimdFloat() {}

        SimdFloat(float f) : simdInternal_(vdupq_n_f32(f)) {}

        // Internal utility constructor to simplify return statements
        SimdFloat(float32x4_t simd) : simdInternal_(simd) {}

        float32x4_t  simdInternal_;
};

class SimdFInt32
{
    public:
        SimdFInt32() {}

        SimdFInt32(std::int32_t i) : simdInternal_(vdupq_n_s32(i)) {}

        // Internal utility constructor to simplify return statements
        SimdFInt32(int32x4_t simd) : simdInternal_(simd) {}

        int32x4_t  simdInternal_;
};

class SimdFBool
{
    public:
        SimdFBool() {}

        SimdFBool(bool b) : simdInternal_(vdupq_n_u32( b ? 0xFFFFFFFF : 0)) {}

        // Internal utility constructor to simplify return statements
        SimdFBool(uint32x4_t simd) : simdInternal_(simd) {}

        uint32x4_t  simdInternal_;
};

class SimdFIBool
{
    public:
        SimdFIBool() {}

        SimdFIBool(bool b) : simdInternal_(vdupq_n_u32( b ? 0xFFFFFFFF : 0)) {}

        // Internal utility constructor to simplify return statements
        SimdFIBool(uint32x4_t simd) : simdInternal_(simd) {}

        uint32x4_t  simdInternal_;
};

static inline SimdFloat gmx_simdcall
simdLoad(const float *m, SimdFloatTag = {})
{
    assert(std::size_t(m) % 16 == 0);
    return {
               vld1q_f32(m)
    };
}

static inline void gmx_simdcall
store(float *m, SimdFloat a)
{
    assert(std::size_t(m) % 16 == 0);
    vst1q_f32(m, a.simdInternal_);
}

static inline SimdFloat gmx_simdcall
simdLoadU(const float *m, SimdFloatTag = {})
{
    return {
               vld1q_f32(m)
    };
}

static inline void gmx_simdcall
storeU(float *m, SimdFloat a)
{
    vst1q_f32(m, a.simdInternal_);
}

static inline SimdFloat gmx_simdcall
setZeroF()
{
    return {
               vdupq_n_f32(0.0f)
    };
}

static inline SimdFInt32 gmx_simdcall
simdLoad(const std::int32_t * m, SimdFInt32Tag)
{
    assert(std::size_t(m) % 16 == 0);
    return {
               vld1q_s32(m)
    };
}

static inline void gmx_simdcall
store(std::int32_t * m, SimdFInt32 a)
{
    assert(std::size_t(m) % 16 == 0);
    vst1q_s32(m, a.simdInternal_);
}

static inline SimdFInt32 gmx_simdcall
simdLoadU(const std::int32_t *m, SimdFInt32Tag)
{
    return {
               vld1q_s32(m)
    };
}

static inline void gmx_simdcall
storeU(std::int32_t * m, SimdFInt32 a)
{
    vst1q_s32(m, a.simdInternal_);
}

static inline SimdFInt32 gmx_simdcall
setZeroFI()
{
    return {
               vdupq_n_s32(0)
    };
}

template<int index> gmx_simdcall
static inline std::int32_t
extract(SimdFInt32 a)
{
    return vgetq_lane_s32(a.simdInternal_, index);
}

static inline SimdFloat gmx_simdcall
operator&(SimdFloat a, SimdFloat b)
{
    return {
               vreinterpretq_f32_s32(vandq_s32(vreinterpretq_s32_f32(a.simdInternal_),
                                               vreinterpretq_s32_f32(b.simdInternal_)))
    };
}

static inline SimdFloat gmx_simdcall
andNot(SimdFloat a, SimdFloat b)
{
    return {
               vreinterpretq_f32_s32(vbicq_s32(vreinterpretq_s32_f32(b.simdInternal_),
                                               vreinterpretq_s32_f32(a.simdInternal_)))
    };
}

static inline SimdFloat gmx_simdcall
operator|(SimdFloat a, SimdFloat b)
{
    return {
               vreinterpretq_f32_s32(vorrq_s32(vreinterpretq_s32_f32(a.simdInternal_),
                                               vreinterpretq_s32_f32(b.simdInternal_)))
    };
}

static inline SimdFloat gmx_simdcall
operator^(SimdFloat a, SimdFloat b)
{
    return {
               vreinterpretq_f32_s32(veorq_s32(vreinterpretq_s32_f32(a.simdInternal_),
                                               vreinterpretq_s32_f32(b.simdInternal_)))
    };
}

static inline SimdFloat gmx_simdcall
operator+(SimdFloat a, SimdFloat b)
{
    return {
               vaddq_f32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator-(SimdFloat a, SimdFloat b)
{
    return {
               vsubq_f32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator-(SimdFloat x)
{
    return {
               vnegq_f32(x.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator*(SimdFloat a, SimdFloat b)
{
    return {
               vmulq_f32(a.simdInternal_, b.simdInternal_)
    };
}

// Override for Neon-Asimd
#if GMX_SIMD_ARM_NEON
static inline SimdFloat gmx_simdcall
fma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
#ifdef __ARM_FEATURE_FMA
               vfmaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)
#else
               vmlaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)
#endif
    };
}

static inline SimdFloat gmx_simdcall
fms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
#ifdef __ARM_FEATURE_FMA
               vnegq_f32(vfmsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_))
#else
               vnegq_f32(vmlsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_))
#endif
    };
}

static inline SimdFloat gmx_simdcall
fnma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
#ifdef __ARM_FEATURE_FMA
               vfmsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)
#else
               vmlsq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_)
#endif
    };
}

static inline SimdFloat gmx_simdcall
fnms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
#ifdef __ARM_FEATURE_FMA
               vnegq_f32(vfmaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_))
#else
               vnegq_f32(vmlaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_))
#endif
    };
}
#endif

static inline SimdFloat gmx_simdcall
rsqrt(SimdFloat x)
{
    return {
               vrsqrteq_f32(x.simdInternal_)
    };
}

// The SIMD implementation seems to overflow when we square lu for
// values close to FLOAT_MAX, so we fall back on the version in
// simd_math.h, which is probably slightly slower.
#if GMX_SIMD_HAVE_NATIVE_RSQRT_ITER_FLOAT
static inline SimdFloat gmx_simdcall
rsqrtIter(SimdFloat lu, SimdFloat x)
{
    return {
               vmulq_f32(lu.simdInternal_, vrsqrtsq_f32(vmulq_f32(lu.simdInternal_, lu.simdInternal_), x.simdInternal_))
    };
}
#endif

static inline SimdFloat gmx_simdcall
rcp(SimdFloat x)
{
    return {
               vrecpeq_f32(x.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
rcpIter(SimdFloat lu, SimdFloat x)
{
    return {
               vmulq_f32(lu.simdInternal_, vrecpsq_f32(lu.simdInternal_, x.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
maskAdd(SimdFloat a, SimdFloat b, SimdFBool m)
{
    b.simdInternal_ = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(b.simdInternal_),
                                                      m.simdInternal_));

    return {
               vaddq_f32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
maskzMul(SimdFloat a, SimdFloat b, SimdFBool m)
{
    SimdFloat tmp = a * b;

    return {
               vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(tmp.simdInternal_),
                                               m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
maskzFma(SimdFloat a, SimdFloat b, SimdFloat c, SimdFBool m)
{
#ifdef __ARM_FEATURE_FMA
    float32x4_t tmp = vfmaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_);
#else
    float32x4_t tmp = vmlaq_f32(c.simdInternal_, b.simdInternal_, a.simdInternal_);
#endif

    return {
               vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(tmp),
                                               m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
maskzRsqrt(SimdFloat x, SimdFBool m)
{
    // The result will always be correct since we mask the result with m, but
    // for debug builds we also want to make sure not to generate FP exceptions
#ifndef NDEBUG
    x.simdInternal_ = vbslq_f32(m.simdInternal_, x.simdInternal_, vdupq_n_f32(1.0f));
#endif
    return {
               vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(vrsqrteq_f32(x.simdInternal_)),
                                               m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
maskzRcp(SimdFloat x, SimdFBool m)
{
    // The result will always be correct since we mask the result with m, but
    // for debug builds we also want to make sure not to generate FP exceptions
#ifndef NDEBUG
    x.simdInternal_ = vbslq_f32(m.simdInternal_, x.simdInternal_, vdupq_n_f32(1.0f));
#endif
    return {
               vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(vrecpeq_f32(x.simdInternal_)),
                                               m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
abs(SimdFloat x)
{
    return {
               vabsq_f32( x.simdInternal_ )
    };
}

static inline SimdFloat gmx_simdcall
max(SimdFloat a, SimdFloat b)
{
    return {
               vmaxq_f32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
min(SimdFloat a, SimdFloat b)
{
    return {
               vminq_f32(a.simdInternal_, b.simdInternal_)
    };
}

// Round and trunc operations are defined at the end of this file, since they
// need to use float-to-integer and integer-to-float conversions.

static inline SimdFloat gmx_simdcall
frexp(SimdFloat value, SimdFInt32 * exponent)
{
    const int32x4_t    exponentMask   = vdupq_n_s32(0x7F800000);
    const int32x4_t    mantissaMask   = vdupq_n_s32(0x807FFFFF);
    const int32x4_t    exponentBias   = vdupq_n_s32(126); // add 1 to make our definition identical to frexp()
    const float32x4_t  half           = vdupq_n_f32(0.5f);
    int32x4_t          iExponent;

    iExponent               = vandq_s32(vreinterpretq_s32_f32(value.simdInternal_), exponentMask);
    iExponent               = vsubq_s32(vshrq_n_s32(iExponent, 23), exponentBias);
    exponent->simdInternal_ = iExponent;

    return {
               vreinterpretq_f32_s32(vorrq_s32(vandq_s32(vreinterpretq_s32_f32(value.simdInternal_),
                                                         mantissaMask),
                                               vreinterpretq_s32_f32(half)))
    };
}

template <MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall
ldexp(SimdFloat value, SimdFInt32 exponent)
{
    const int32x4_t exponentBias = vdupq_n_s32(127);
    int32x4_t       iExponent    = vaddq_s32(exponent.simdInternal_, exponentBias);

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponent = vmaxq_s32(iExponent, vdupq_n_s32(0));
    }

    iExponent = vshlq_n_s32( iExponent, 23);

    return {
               vmulq_f32(value.simdInternal_, vreinterpretq_f32_s32(iExponent))
    };
}

// Override for Neon-Asimd
#if GMX_SIMD_ARM_NEON
static inline float gmx_simdcall
reduce(SimdFloat a)
{
    float32x4_t x = a.simdInternal_;
    float32x4_t y = vextq_f32(x, x, 2);

    x = vaddq_f32(x, y);
    y = vextq_f32(x, x, 1);
    x = vaddq_f32(x, y);
    return vgetq_lane_f32(x, 0);
}
#endif

static inline SimdFBool gmx_simdcall
operator==(SimdFloat a, SimdFloat b)
{
    return {
               vceqq_f32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
operator!=(SimdFloat a, SimdFloat b)
{
    return {
               vmvnq_u32(vceqq_f32(a.simdInternal_, b.simdInternal_))
    };
}

static inline SimdFBool gmx_simdcall
operator<(SimdFloat a, SimdFloat b)
{
    return {
               vcltq_f32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
operator<=(SimdFloat a, SimdFloat b)
{
    return {
               vcleq_f32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
testBits(SimdFloat a)
{
    uint32x4_t tmp = vreinterpretq_u32_f32(a.simdInternal_);

    return {
               vtstq_u32(tmp, tmp)
    };
}

static inline SimdFBool gmx_simdcall
operator&&(SimdFBool a, SimdFBool b)
{

    return {
               vandq_u32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
operator||(SimdFBool a, SimdFBool b)
{
    return {
               vorrq_u32(a.simdInternal_, b.simdInternal_)
    };
}

// Override for Neon-Asimd
#if GMX_SIMD_ARM_NEON
static inline bool gmx_simdcall
anyTrue(SimdFBool a)
{
    uint32x4_t x = a.simdInternal_;
    uint32x4_t y = vextq_u32(x, x, 2);

    x = vorrq_u32(x, y);
    y = vextq_u32(x, x, 1);
    x = vorrq_u32(x, y);
    return (vgetq_lane_u32(x, 0) != 0);
}
#endif

static inline SimdFloat gmx_simdcall
selectByMask(SimdFloat a, SimdFBool m)
{
    return {
               vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(a.simdInternal_),
                                               m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
selectByNotMask(SimdFloat a, SimdFBool m)
{
    return {
               vreinterpretq_f32_u32(vbicq_u32(vreinterpretq_u32_f32(a.simdInternal_),
                                               m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
blend(SimdFloat a, SimdFloat b, SimdFBool sel)
{
    return {
               vbslq_f32(sel.simdInternal_, b.simdInternal_, a.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator<<(SimdFInt32 a, int n)
{
    return {
               vshlq_s32(a.simdInternal_, vdupq_n_s32(n >= 32 ? 32 : n))
    };
}

static inline SimdFInt32 gmx_simdcall
operator>>(SimdFInt32 a, int n)
{
    return {
               vshlq_s32(a.simdInternal_, vdupq_n_s32(n >= 32 ? -32 : -n))
    };
}

static inline SimdFInt32 gmx_simdcall
operator&(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vandq_s32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
andNot(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vbicq_s32(b.simdInternal_, a.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator|(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vorrq_s32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator^(SimdFInt32 a, SimdFInt32 b)
{
    return {
               veorq_s32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator+(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vaddq_s32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator-(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vsubq_s32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator*(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vmulq_s32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator==(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vceqq_s32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
testBits(SimdFInt32 a)
{
    return {
               vtstq_s32(a.simdInternal_, a.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator<(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vcltq_s32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator&&(SimdFIBool a, SimdFIBool b)
{
    return {
               vandq_u32(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator||(SimdFIBool a, SimdFIBool b)
{
    return {
               vorrq_u32(a.simdInternal_, b.simdInternal_)
    };
}

// Override for Neon-Asimd
#if GMX_SIMD_ARM_NEON
static inline bool gmx_simdcall
anyTrue(SimdFIBool a)
{
    uint32x4_t x = a.simdInternal_;
    uint32x4_t y = vextq_u32(x, x, 2);

    x = vorrq_u32(x, y);
    y = vextq_u32(x, x, 1);
    x = vorrq_u32(x, y);
    return (vgetq_lane_u32(x, 0) != 0);
}
#endif

static inline SimdFInt32 gmx_simdcall
selectByMask(SimdFInt32 a, SimdFIBool m)
{
    return {
               vandq_s32(a.simdInternal_, vreinterpretq_s32_u32(m.simdInternal_))
    };
}

static inline SimdFInt32 gmx_simdcall
selectByNotMask(SimdFInt32 a, SimdFIBool m)
{
    return {
               vbicq_s32(a.simdInternal_, vreinterpretq_s32_u32(m.simdInternal_))
    };
}

static inline SimdFInt32 gmx_simdcall
blend(SimdFInt32 a, SimdFInt32 b, SimdFIBool sel)
{
    return {
               vbslq_s32(sel.simdInternal_, b.simdInternal_, a.simdInternal_)
    };
}

// Override for Neon-Asimd
#if GMX_SIMD_ARM_NEON
static inline SimdFInt32 gmx_simdcall
cvtR2I(SimdFloat a)
{
    float32x4_t signBitOfA = vreinterpretq_f32_u32(vandq_u32(vdupq_n_u32(0x80000000), vreinterpretq_u32_f32(a.simdInternal_)));
    float32x4_t half       = vdupq_n_f32(0.5f);
    float32x4_t corr       = vreinterpretq_f32_u32(vorrq_u32(vreinterpretq_u32_f32(half), vreinterpretq_u32_f32(signBitOfA)));

    return {
               vcvtq_s32_f32(vaddq_f32(a.simdInternal_, corr))
    };
}
#endif

static inline SimdFInt32 gmx_simdcall
cvttR2I(SimdFloat a)
{
    return {
               vcvtq_s32_f32(a.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
cvtI2R(SimdFInt32 a)
{
    return {
               vcvtq_f32_s32(a.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
cvtB2IB(SimdFBool a)
{
    return {
               a.simdInternal_
    };
}

static inline SimdFBool gmx_simdcall
cvtIB2B(SimdFIBool a)
{
    return {
               a.simdInternal_
    };
}

// Override for Neon-Asimd
#if GMX_SIMD_ARM_NEON
static inline SimdFloat gmx_simdcall
round(SimdFloat x)
{
    return cvtI2R(cvtR2I(x));
}

static inline SimdFloat gmx_simdcall
trunc(SimdFloat x)
{
    return cvtI2R(cvttR2I(x));
}
#endif

}      // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_NEON_SIMD_FLOAT_H
