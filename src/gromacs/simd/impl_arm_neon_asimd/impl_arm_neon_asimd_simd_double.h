/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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

#ifndef GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD_DOUBLE_H

#include "config.h"

#include <arm_neon.h>

#include <cassert>

#include "gromacs/math/utilities.h"

#include "impl_arm_neon_asimd_simd_float.h"

namespace gmx
{

class SimdDouble
{
public:
    SimdDouble() {}

    SimdDouble(double d) : simdInternal_(vdupq_n_f64(d)) {}

    // Internal utility constructor to simplify return statements
    SimdDouble(float64x2_t simd) : simdInternal_(simd) {}

    float64x2_t simdInternal_;
};

class SimdDInt32
{
public:
    SimdDInt32() {}

    SimdDInt32(std::int32_t i) : simdInternal_(vdup_n_s32(i)) {}

    // Internal utility constructor to simplify return statements
    SimdDInt32(int32x2_t simd) : simdInternal_(simd) {}

    int32x2_t simdInternal_;
};

class SimdDBool
{
public:
    SimdDBool() {}

    SimdDBool(bool b) : simdInternal_(vdupq_n_u64(b ? 0xFFFFFFFFFFFFFFFF : 0)) {}

    // Internal utility constructor to simplify return statements
    SimdDBool(uint64x2_t simd) : simdInternal_(simd) {}

    uint64x2_t simdInternal_;
};

class SimdDIBool
{
public:
    SimdDIBool() {}

    SimdDIBool(bool b) : simdInternal_(vdup_n_u32(b ? 0xFFFFFFFF : 0)) {}

    // Internal utility constructor to simplify return statements
    SimdDIBool(uint32x2_t simd) : simdInternal_(simd) {}

    uint32x2_t simdInternal_;
};

static inline SimdDouble gmx_simdcall simdLoad(const double* m, SimdDoubleTag = {})
{
    assert(std::size_t(m) % 16 == 0);
    return { vld1q_f64(m) };
}

static inline void gmx_simdcall store(double* m, SimdDouble a)
{
    assert(std::size_t(m) % 16 == 0);
    vst1q_f64(m, a.simdInternal_);
}

static inline SimdDouble gmx_simdcall simdLoadU(const double* m, SimdDoubleTag = {})
{
    return { vld1q_f64(m) };
}

static inline void gmx_simdcall storeU(double* m, SimdDouble a)
{
    vst1q_f64(m, a.simdInternal_);
}

static inline SimdDouble gmx_simdcall setZeroD()
{
    return { vdupq_n_f64(0.0) };
}

static inline SimdDInt32 gmx_simdcall simdLoad(const std::int32_t* m, SimdDInt32Tag)
{
    assert(std::size_t(m) % 8 == 0);
    return { vld1_s32(m) };
}

static inline void gmx_simdcall store(std::int32_t* m, SimdDInt32 a)
{
    assert(std::size_t(m) % 8 == 0);
    vst1_s32(m, a.simdInternal_);
}

static inline SimdDInt32 gmx_simdcall simdLoadU(const std::int32_t* m, SimdDInt32Tag)
{
    return { vld1_s32(m) };
}

static inline void gmx_simdcall storeU(std::int32_t* m, SimdDInt32 a)
{
    vst1_s32(m, a.simdInternal_);
}

static inline SimdDInt32 gmx_simdcall setZeroDI()
{
    return { vdup_n_s32(0) };
}

template<int index>
gmx_simdcall static inline std::int32_t extract(SimdDInt32 a)
{
    return vget_lane_s32(a.simdInternal_, index);
}

static inline SimdDouble gmx_simdcall operator&(SimdDouble a, SimdDouble b)
{
    return { vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(a.simdInternal_),
                                             vreinterpretq_u64_f64(b.simdInternal_))) };
}

static inline SimdDouble gmx_simdcall andNot(SimdDouble a, SimdDouble b)
{
    return { vreinterpretq_f64_u64(vbicq_u64(vreinterpretq_u64_f64(b.simdInternal_),
                                             vreinterpretq_u64_f64(a.simdInternal_))) };
}

static inline SimdDouble gmx_simdcall operator|(SimdDouble a, SimdDouble b)
{
    return { vreinterpretq_f64_u64(vorrq_u64(vreinterpretq_u64_f64(a.simdInternal_),
                                             vreinterpretq_u64_f64(b.simdInternal_))) };
}

static inline SimdDouble gmx_simdcall operator^(SimdDouble a, SimdDouble b)
{
    return { vreinterpretq_f64_u64(veorq_u64(vreinterpretq_u64_f64(a.simdInternal_),
                                             vreinterpretq_u64_f64(b.simdInternal_))) };
}

static inline SimdDouble gmx_simdcall operator+(SimdDouble a, SimdDouble b)
{
    return { vaddq_f64(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDouble gmx_simdcall operator-(SimdDouble a, SimdDouble b)
{
    return { vsubq_f64(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDouble gmx_simdcall operator-(SimdDouble x)
{
    return { vnegq_f64(x.simdInternal_) };
}

static inline SimdDouble gmx_simdcall operator*(SimdDouble a, SimdDouble b)
{
    return { vmulq_f64(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDouble gmx_simdcall fma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return { vfmaq_f64(c.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline SimdDouble gmx_simdcall fms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return { vnegq_f64(vfmsq_f64(c.simdInternal_, b.simdInternal_, a.simdInternal_)) };
}

static inline SimdDouble gmx_simdcall fnma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return { vfmsq_f64(c.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline SimdDouble gmx_simdcall fnms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return { vnegq_f64(vfmaq_f64(c.simdInternal_, b.simdInternal_, a.simdInternal_)) };
}

static inline SimdDouble gmx_simdcall rsqrt(SimdDouble x)
{
    return { vrsqrteq_f64(x.simdInternal_) };
}

static inline SimdDouble gmx_simdcall rsqrtIter(SimdDouble lu, SimdDouble x)
{
    return { vmulq_f64(lu.simdInternal_,
                       vrsqrtsq_f64(vmulq_f64(lu.simdInternal_, lu.simdInternal_), x.simdInternal_)) };
}

static inline SimdDouble gmx_simdcall rcp(SimdDouble x)
{
    return { vrecpeq_f64(x.simdInternal_) };
}

static inline SimdDouble gmx_simdcall rcpIter(SimdDouble lu, SimdDouble x)
{
    return { vmulq_f64(lu.simdInternal_, vrecpsq_f64(lu.simdInternal_, x.simdInternal_)) };
}

static inline SimdDouble gmx_simdcall maskAdd(SimdDouble a, SimdDouble b, SimdDBool m)
{
    float64x2_t addend =
            vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(b.simdInternal_), m.simdInternal_));

    return { vaddq_f64(a.simdInternal_, addend) };
}

static inline SimdDouble gmx_simdcall maskzMul(SimdDouble a, SimdDouble b, SimdDBool m)
{
    float64x2_t prod = vmulq_f64(a.simdInternal_, b.simdInternal_);
    return { vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(prod), m.simdInternal_)) };
}

static inline SimdDouble gmx_simdcall maskzFma(SimdDouble a, SimdDouble b, SimdDouble c, SimdDBool m)
{
    float64x2_t prod = vfmaq_f64(c.simdInternal_, b.simdInternal_, a.simdInternal_);

    return { vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(prod), m.simdInternal_)) };
}

static inline SimdDouble gmx_simdcall maskzRsqrt(SimdDouble x, SimdDBool m)
{
    // The result will always be correct since we mask the result with m, but
    // for debug builds we also want to make sure not to generate FP exceptions
#ifndef NDEBUG
    x.simdInternal_ = vbslq_f64(m.simdInternal_, x.simdInternal_, vdupq_n_f64(1.0));
#endif
    return { vreinterpretq_f64_u64(
            vandq_u64(vreinterpretq_u64_f64(vrsqrteq_f64(x.simdInternal_)), m.simdInternal_)) };
}

static inline SimdDouble gmx_simdcall maskzRcp(SimdDouble x, SimdDBool m)
{
    // The result will always be correct since we mask the result with m, but
    // for debug builds we also want to make sure not to generate FP exceptions
#ifndef NDEBUG
    x.simdInternal_ = vbslq_f64(m.simdInternal_, x.simdInternal_, vdupq_n_f64(1.0));
#endif
    return { vreinterpretq_f64_u64(
            vandq_u64(vreinterpretq_u64_f64(vrecpeq_f64(x.simdInternal_)), m.simdInternal_)) };
}

static inline SimdDouble gmx_simdcall abs(SimdDouble x)
{
    return { vabsq_f64(x.simdInternal_) };
}

static inline SimdDouble gmx_simdcall max(SimdDouble a, SimdDouble b)
{
    return { vmaxq_f64(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDouble gmx_simdcall min(SimdDouble a, SimdDouble b)
{
    return { vminq_f64(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDouble gmx_simdcall round(SimdDouble x)
{
    return { vrndnq_f64(x.simdInternal_) };
}

static inline SimdDouble gmx_simdcall trunc(SimdDouble x)
{
    return { vrndq_f64(x.simdInternal_) };
}

template<MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble frexp(SimdDouble value, SimdDInt32* exponent)
{
    const float64x2_t exponentMask = vreinterpretq_f64_u64(vdupq_n_u64(0x7FF0000000000000LL));
    const float64x2_t mantissaMask = vreinterpretq_f64_u64(vdupq_n_u64(0x800FFFFFFFFFFFFFLL));

    const int64x2_t exponentBias = vdupq_n_s64(1022); // add 1 to make our definition identical to frexp()
    const float64x2_t half = vdupq_n_f64(0.5);
    int64x2_t         iExponent;

    iExponent = vandq_s64(vreinterpretq_s64_f64(value.simdInternal_), vreinterpretq_s64_f64(exponentMask));
    iExponent = vsubq_s64(vshrq_n_s64(iExponent, 52), exponentBias);

    float64x2_t result = vreinterpretq_f64_u64(vorrq_u64(
            vandq_u64(vreinterpretq_u64_f64(value.simdInternal_), vreinterpretq_u64_f64(mantissaMask)),
            vreinterpretq_u64_f64(half)));

    if (opt == MathOptimization::Safe)
    {
        uint64x2_t valueIsZero = vceqq_f64(value.simdInternal_, vdupq_n_f64(0.0));
        iExponent              = vbicq_s64(iExponent, vreinterpretq_s64_u64(valueIsZero));
        result                 = vbslq_f64(valueIsZero, value.simdInternal_, result);
    }

    exponent->simdInternal_ = vmovn_s64(iExponent);

    return { result };
}

template<MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble ldexp(SimdDouble value, SimdDInt32 exponent)
{
    const int32x2_t exponentBias = vdup_n_s32(1023);
    int32x2_t       iExponent    = vadd_s32(exponent.simdInternal_, exponentBias);
    int64x2_t       iExponent64;

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponent = vmax_s32(iExponent, vdup_n_s32(0));
    }

    iExponent64 = vmovl_s32(iExponent);
    iExponent64 = vshlq_n_s64(iExponent64, 52);

    return { vmulq_f64(value.simdInternal_, vreinterpretq_f64_s64(iExponent64)) };
}

static inline double gmx_simdcall reduce(SimdDouble a)
{
    float64x2_t b = vpaddq_f64(a.simdInternal_, a.simdInternal_);
    return vgetq_lane_f64(b, 0);
}

static inline SimdDBool gmx_simdcall operator==(SimdDouble a, SimdDouble b)
{
    return { vceqq_f64(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDBool gmx_simdcall operator!=(SimdDouble a, SimdDouble b)
{
    return { vreinterpretq_u64_u32(
            vmvnq_u32(vreinterpretq_u32_u64(vceqq_f64(a.simdInternal_, b.simdInternal_)))) };
}

static inline SimdDBool gmx_simdcall operator<(SimdDouble a, SimdDouble b)
{
    return { vcltq_f64(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDBool gmx_simdcall operator<=(SimdDouble a, SimdDouble b)
{
    return { vcleq_f64(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDBool gmx_simdcall testBits(SimdDouble a)
{
    return { vtstq_s64(vreinterpretq_s64_f64(a.simdInternal_), vreinterpretq_s64_f64(a.simdInternal_)) };
}

static inline SimdDBool gmx_simdcall operator&&(SimdDBool a, SimdDBool b)
{
    return { vandq_u64(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDBool gmx_simdcall operator||(SimdDBool a, SimdDBool b)
{
    return { vorrq_u64(a.simdInternal_, b.simdInternal_) };
}

static inline bool gmx_simdcall anyTrue(SimdDBool a)
{
    return (vmaxvq_u32(vreinterpretq_u32_u64(a.simdInternal_)) != 0);
}

static inline SimdDouble gmx_simdcall selectByMask(SimdDouble a, SimdDBool m)
{
    return { vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(a.simdInternal_), m.simdInternal_)) };
}

static inline SimdDouble gmx_simdcall selectByNotMask(SimdDouble a, SimdDBool m)
{
    return { vreinterpretq_f64_u64(vbicq_u64(vreinterpretq_u64_f64(a.simdInternal_), m.simdInternal_)) };
}

static inline SimdDouble gmx_simdcall blend(SimdDouble a, SimdDouble b, SimdDBool sel)
{
    return { vbslq_f64(sel.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall operator&(SimdDInt32 a, SimdDInt32 b)
{
    return { vand_s32(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall andNot(SimdDInt32 a, SimdDInt32 b)
{
    return { vbic_s32(b.simdInternal_, a.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall operator|(SimdDInt32 a, SimdDInt32 b)
{
    return { vorr_s32(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall operator^(SimdDInt32 a, SimdDInt32 b)
{
    return { veor_s32(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall operator+(SimdDInt32 a, SimdDInt32 b)
{
    return { vadd_s32(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall operator-(SimdDInt32 a, SimdDInt32 b)
{
    return { vsub_s32(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall operator*(SimdDInt32 a, SimdDInt32 b)
{
    return { vmul_s32(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDIBool gmx_simdcall operator==(SimdDInt32 a, SimdDInt32 b)
{
    return { vceq_s32(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDIBool gmx_simdcall testBits(SimdDInt32 a)
{
    return { vtst_s32(a.simdInternal_, a.simdInternal_) };
}

static inline SimdDIBool gmx_simdcall operator<(SimdDInt32 a, SimdDInt32 b)
{
    return { vclt_s32(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDIBool gmx_simdcall operator&&(SimdDIBool a, SimdDIBool b)
{
    return { vand_u32(a.simdInternal_, b.simdInternal_) };
}

static inline SimdDIBool gmx_simdcall operator||(SimdDIBool a, SimdDIBool b)
{
    return { vorr_u32(a.simdInternal_, b.simdInternal_) };
}

static inline bool gmx_simdcall anyTrue(SimdDIBool a)
{
    return (vmaxv_u32(a.simdInternal_) != 0);
}

static inline SimdDInt32 gmx_simdcall selectByMask(SimdDInt32 a, SimdDIBool m)
{
    return { vand_s32(a.simdInternal_, vreinterpret_s32_u32(m.simdInternal_)) };
}

static inline SimdDInt32 gmx_simdcall selectByNotMask(SimdDInt32 a, SimdDIBool m)
{
    return { vbic_s32(a.simdInternal_, vreinterpret_s32_u32(m.simdInternal_)) };
}

static inline SimdDInt32 gmx_simdcall blend(SimdDInt32 a, SimdDInt32 b, SimdDIBool sel)
{
    return { vbsl_s32(sel.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall cvtR2I(SimdDouble a)
{
    return { vmovn_s64(vcvtnq_s64_f64(a.simdInternal_)) };
}

static inline SimdDInt32 gmx_simdcall cvttR2I(SimdDouble a)
{
    return { vmovn_s64(vcvtq_s64_f64(a.simdInternal_)) };
}

static inline SimdDouble gmx_simdcall cvtI2R(SimdDInt32 a)
{
    return { vcvtq_f64_s64(vmovl_s32(a.simdInternal_)) };
}

static inline SimdDIBool gmx_simdcall cvtB2IB(SimdDBool a)
{
    return { vqmovn_u64(a.simdInternal_) };
}

static inline SimdDBool gmx_simdcall cvtIB2B(SimdDIBool a)
{
    return { vorrq_u64(vmovl_u32(a.simdInternal_), vshlq_n_u64(vmovl_u32(a.simdInternal_), 32)) };
}

static inline void gmx_simdcall cvtF2DD(SimdFloat f, SimdDouble* d0, SimdDouble* d1)
{
    d0->simdInternal_ = vcvt_f64_f32(vget_low_f32(f.simdInternal_));
    d1->simdInternal_ = vcvt_high_f64_f32(f.simdInternal_);
}

static inline SimdFloat gmx_simdcall cvtDD2F(SimdDouble d0, SimdDouble d1)
{
    return { vcvt_high_f32_f64(vcvt_f32_f64(d0.simdInternal_), d1.simdInternal_) };
}

} // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_NEON_ASIMD_SIMD_DOUBLE_H
