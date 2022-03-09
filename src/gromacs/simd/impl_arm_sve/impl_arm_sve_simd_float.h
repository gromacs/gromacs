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

#ifndef GMX_SIMD_IMPL_ARM_SVE_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_ARM_SVE_SIMD_FLOAT_H

#include "config.h"

#include <arm_sve.h>

#include <cassert>
#include <cstddef>
#include <cstdint>

#include "gromacs/math/utilities.h"

namespace gmx
{


class SimdFloat
{
private:
    typedef svfloat32_t simdInternalType_
            __attribute__((arm_sve_vector_bits(GMX_SIMD_ARM_SVE_LENGTH_VALUE)));

public:
    SimdFloat() {}

    SimdFloat(const float f) { this->simdInternal_ = svdup_n_f32(f); }

    SimdFloat(svfloat32_t simd) : simdInternal_(simd) {}

    simdInternalType_ simdInternal_;
};

class SimdFInt32
{
private:
    typedef svint32_t simdInternalType_
            __attribute__((arm_sve_vector_bits(GMX_SIMD_ARM_SVE_LENGTH_VALUE)));

public:
    SimdFInt32() {}

    SimdFInt32(const int32_t i) { this->simdInternal_ = svdup_n_s32(i); }

    SimdFInt32(svint32_t simd) : simdInternal_(simd) {}

    simdInternalType_ simdInternal_;
};

class SimdFBool
{
private:
    typedef svbool_t simdInternalType_
            __attribute__((arm_sve_vector_bits(GMX_SIMD_ARM_SVE_LENGTH_VALUE)));

public:
    SimdFBool() {}

    SimdFBool(const bool b) { this->simdInternal_ = svdup_n_b32(b); }

    SimdFBool(svbool_t simd) : simdInternal_(simd) {}

    simdInternalType_ simdInternal_;
};

class SimdFIBool
{
private:
    typedef svbool_t simdInternalType_
            __attribute__((arm_sve_vector_bits(GMX_SIMD_ARM_SVE_LENGTH_VALUE)));

public:
    SimdFIBool() {}

    SimdFIBool(const bool b) { this->simdInternal_ = svdup_n_b32(b); }

    SimdFIBool(svbool_t simd) : simdInternal_(simd) {}

    simdInternalType_ simdInternal_;
};

static inline SimdFloat gmx_simdcall simdLoad(const float* m, SimdFloatTag = {})
{
    assert(0 == (std::size_t(m) % GMX_SIMD_ALIGNMENT));
    svbool_t pg = svptrue_b32();
    return { svld1_f32(pg, m) };
}

static inline SimdFloat gmx_simdcall simdLoad(SimdFloat* m, int offset, SimdFloatTag = {})
{
    assert(0 == (std::size_t(m) % GMX_SIMD_ALIGNMENT));
    svbool_t pg = svptrue_b32();
    return { svld1_f32(pg, reinterpret_cast<float*>(m) + offset * svcntw()) };
}

static inline SimdFloat gmx_simdcall simdLoadFloat(const float* m)
{
    assert(0 == (std::size_t(m) % GMX_SIMD_ALIGNMENT));
    svbool_t pg = svptrue_b32();
    return { svld1_f32(pg, m) };
}

static inline void gmx_simdcall store(float* m, SimdFloat a)
{
    assert(0 == (std::size_t(m) % GMX_SIMD_ALIGNMENT));
    svbool_t pg = svptrue_b32();
    svst1_f32(pg, m, a.simdInternal_);
}

static inline SimdFloat gmx_simdcall simdLoadU(const float* m, SimdFloatTag = {})
{
    svbool_t pg = svptrue_b32();
    return { svld1_f32(pg, m) };
}

static inline void gmx_simdcall storeU(float* m, SimdFloat a)
{
    svbool_t pg = svptrue_b32();
    svst1_f32(pg, m, a.simdInternal_);
}

static inline SimdFloat gmx_simdcall setZeroF()
{
    return { svdup_n_f32(0.0f) };
}

static inline void gmx_simdcall simdIncr(SimdFloat*& p, SimdFloatTag)
{
    p = reinterpret_cast<SimdFloat*>(reinterpret_cast<uint64_t>(p) + svcntw());
}

static inline SimdFInt32 gmx_simdcall simdLoad(const std::int32_t* m, SimdFInt32Tag)
{
    assert(0 == (std::size_t(m) % GMX_SIMD_ALIGNMENT));
    svbool_t pg = svptrue_b32();
    return { svld1_s32(pg, m) };
}

static inline void gmx_simdcall store(std::int32_t* m, SimdFInt32 a)
{
    assert(0 == (std::size_t(m) % GMX_SIMD_ALIGNMENT));
    svbool_t pg = svptrue_b32();
    svst1_s32(pg, m, a.simdInternal_);
}

static inline SimdFInt32 gmx_simdcall simdLoadU(const std::int32_t* m, SimdFInt32Tag)
{
    svbool_t pg = svptrue_b32();
    return { svld1_s32(pg, m) };
}

static inline void gmx_simdcall storeU(std::int32_t* m, SimdFInt32 a)
{
    svbool_t pg = svptrue_b32();
    svst1_s32(pg, m, a.simdInternal_);
}

static inline SimdFInt32 gmx_simdcall setZeroFI()
{
    return { svdup_n_s32(0) };
}

template<int index>
gmx_simdcall static inline std::int32_t extract(SimdFInt32 a)
{
    svbool_t pg = svwhilelt_b32(0, index);
    return svlasta_s32(pg, a.simdInternal_);
}

template<int index>
gmx_simdcall static inline float extract(SimdFloat a)
{
    svbool_t pg = svwhilelt_b32(0, index);
    return svlasta_f32(pg, a.simdInternal_);
}

static inline SimdFloat gmx_simdcall operator&(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svreinterpret_f32_s32(svand_s32_x(
            pg, svreinterpret_s32_f32(a.simdInternal_), svreinterpret_s32_f32(b.simdInternal_))) };
}

static inline SimdFloat gmx_simdcall andNot(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svreinterpret_f32_s32(svbic_s32_x(
            pg, svreinterpret_s32_f32(b.simdInternal_), svreinterpret_s32_f32(a.simdInternal_))) };
}

static inline SimdFloat gmx_simdcall operator|(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svreinterpret_f32_s32(svorr_s32_x(
            pg, svreinterpret_s32_f32(a.simdInternal_), svreinterpret_s32_f32(b.simdInternal_))) };
}

static inline SimdFloat gmx_simdcall operator^(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svreinterpret_f32_s32(sveor_s32_x(
            pg, svreinterpret_s32_f32(a.simdInternal_), svreinterpret_s32_f32(b.simdInternal_))) };
}

static inline SimdFloat gmx_simdcall operator+(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svadd_f32_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFloat gmx_simdcall operator-(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svsub_f32_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFloat gmx_simdcall operator-(SimdFloat a)
{
    svbool_t pg = svptrue_b32();
    return { svneg_f32_x(pg, a.simdInternal_) };
}

static inline SimdFloat gmx_simdcall operator*(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svmul_f32_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFloat gmx_simdcall fma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    svbool_t pg = svptrue_b32();
    return { svmad_f32_x(pg, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline SimdFloat gmx_simdcall fms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    svbool_t pg = svptrue_b32();
    return { svnmsb_f32_x(pg, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline SimdFloat gmx_simdcall fnma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    svbool_t pg = svptrue_b32();
    return { svmsb_f32_x(pg, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline SimdFloat gmx_simdcall fnms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    svbool_t pg = svptrue_b32();
    return { svnmad_f32_x(pg, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline SimdFloat gmx_simdcall rsqrt(SimdFloat x)
{
    return { svrsqrte_f32(x.simdInternal_) };
}

// The SIMD implementation seems to overflow when we square lu for
// values close to FLOAT_MAX, so we fall back on the version in
// simd_math.h, which is probably slightly slower.
#if GMX_SIMD_HAVE_NATIVE_RSQRT_ITER_FLOAT
static inline SimdFloat gmx_simdcall rsqrtIter(SimdFloat lu, SimdFloat x)
{
    svbool_t    pg = svptrue_b32();
    svfloat32_t tmp1, tmp2;
    tmp1 = svmul_f32_x(pg, x.simdInternal_, lu.simdInternal_);
    tmp2 = svmul_n_f32_x(pg, lu.simdInternal_, -0.5f);
    tmp1 = svmad_n_f32_x(pg, tmp1, lu.simdInternal_, -3.0f);
    return { svmul_f32_x(pg, tmp1, tmp2) };
}

#endif

static inline SimdFloat gmx_simdcall rcp(SimdFloat x)
{
    return { svrecpe_f32(x.simdInternal_) };
}

static inline SimdFloat gmx_simdcall rcpIter(SimdFloat lu, SimdFloat x)
{
    svbool_t pg = svptrue_b32();
    return { svmul_f32_x(pg, lu.simdInternal_, svrecps_f32(lu.simdInternal_, x.simdInternal_)) };
}

static inline SimdFloat gmx_simdcall maskAdd(SimdFloat a, SimdFloat b, SimdFBool m)
{
    return { svadd_f32_m(m.simdInternal_, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFloat gmx_simdcall maskzMul(SimdFloat a, SimdFloat b, SimdFBool m)
{
    return { svmul_f32_z(m.simdInternal_, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFloat gmx_simdcall maskzFma(SimdFloat a, SimdFloat b, SimdFloat c, SimdFBool m)
{
    return { svmad_f32_z(m.simdInternal_, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline SimdFloat gmx_simdcall maskzRsqrt(SimdFloat x, SimdFBool m)
{
    // The result will always be correct since we mask the result with m, but
    // for debug builds we also want to make sure not to generate FP exceptions
#ifndef NDEBUG
    x.simdInternal_ = svsel_f32(m.simdInternal_, x.simdInternal_, svdup_n_f32(1.0f));
#endif
    return { svsel_f32(m.simdInternal_, svrsqrte_f32(x.simdInternal_), svdup_n_f32(0.0f)) };
}

static inline SimdFloat gmx_simdcall maskzRcp(SimdFloat x, SimdFBool m)
{
    // The result will always be correct since we mask the result with m, but
    // for debug builds we also want to make sure not to generate FP exceptions
#ifndef NDEBUG
    x.simdInternal_ = svsel_f32(m.simdInternal_, x.simdInternal_, svdup_n_f32(1.0f));
#endif
    return { svsel_f32(m.simdInternal_, svrecpe_f32(x.simdInternal_), svdup_n_f32(0.0f)) };
}

static inline SimdFloat gmx_simdcall abs(SimdFloat x)
{
    svbool_t pg = svptrue_b32();
    return { svabs_f32_x(pg, x.simdInternal_) };
}

static inline SimdFloat gmx_simdcall max(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svmax_f32_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFloat gmx_simdcall min(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svmin_f32_x(pg, a.simdInternal_, b.simdInternal_) };
}

// Round and trunc operations are defined at the end of this file, since they
// need to use float-to-integer and integer-to-float conversions.

template<MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall frexp(SimdFloat value, SimdFInt32* exponent)
{
    svbool_t        pg           = svptrue_b32();
    const svint32_t exponentMask = svdup_n_s32(0x7F800000);
    const svint32_t mantissaMask = svdup_n_s32(0x807FFFFF);
    const svint32_t exponentBias = svdup_n_s32(126); // add 1 to make our definition identical to frexp()
    const svfloat32_t half       = svdup_n_f32(0.5f);
    svint32_t         iExponent;

    iExponent = svand_s32_x(pg, svreinterpret_s32_f32(value.simdInternal_), exponentMask);
    iExponent = svsub_s32_x(
            pg, svreinterpret_s32_u32(svlsr_n_u32_x(pg, svreinterpret_u32_s32(iExponent), 23)), exponentBias);

    svfloat32_t result = svreinterpret_f32_s32(
            svorr_s32_x(pg,
                        svand_s32_x(pg, svreinterpret_s32_f32(value.simdInternal_), mantissaMask),
                        svreinterpret_s32_f32(half)));

    if (opt == MathOptimization::Safe)
    {
        svbool_t valueIsZero = svcmpeq_n_f32(pg, value.simdInternal_, 0.0F);
        iExponent            = svsel_s32(valueIsZero, svdup_n_s32(0), iExponent);
        result               = svsel_f32(valueIsZero, value.simdInternal_, result);
    }

    exponent->simdInternal_ = iExponent;
    return { result };
}

template<MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall ldexp(SimdFloat value, SimdFInt32 exponent)
{
    svbool_t        pg           = svptrue_b32();
    const svint32_t exponentBias = svdup_n_s32(127);
    svint32_t       iExponent    = svadd_s32_x(pg, exponent.simdInternal_, exponentBias);

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponent = svmax_n_s32_x(pg, iExponent, 0);
    }

    iExponent = svlsl_n_s32_x(pg, iExponent, 23);

    return { svmul_f32_x(pg, value.simdInternal_, svreinterpret_f32_s32(iExponent)) };
}

static inline float gmx_simdcall reduce(SimdFloat a)
{
    svbool_t pg = svptrue_b32();
    return svaddv_f32(pg, a.simdInternal_);
}

static inline SimdFBool gmx_simdcall operator==(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svcmpeq_f32(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFBool gmx_simdcall operator!=(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svcmpne_f32(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFBool gmx_simdcall operator<(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svcmplt_f32(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFBool gmx_simdcall operator<=(SimdFloat a, SimdFloat b)
{
    svbool_t pg = svptrue_b32();
    return { svcmple_f32(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFBool gmx_simdcall testBits(SimdFloat a)
{
    svbool_t pg = svptrue_b32();
    return { svcmpne_n_s32(pg, svreinterpret_s32_f32(a.simdInternal_), 0) };
}

static inline SimdFBool gmx_simdcall operator&&(SimdFBool a, SimdFBool b)
{
    svbool_t pg = svptrue_b32();
    return { svand_b_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFBool gmx_simdcall operator||(SimdFBool a, SimdFBool b)
{
    svbool_t pg = svptrue_b32();
    return { svorr_b_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline bool gmx_simdcall anyTrue(SimdFBool a)
{
    svbool_t pg = svptrue_b32();
    return svptest_any(pg, a.simdInternal_);
}

static inline bool gmx_simdcall extractFirst(SimdFBool a)
{
    svbool_t pg = svptrue_b32();
    return svptest_first(pg, a.simdInternal_);
}

static inline SimdFloat gmx_simdcall selectByMask(SimdFloat a, SimdFBool m)
{
    return { svsel_f32(m.simdInternal_, a.simdInternal_, svdup_n_f32(0.0f)) };
}

static inline SimdFloat gmx_simdcall selectByNotMask(SimdFloat a, SimdFBool m)
{
    return { svsel_f32(m.simdInternal_, svdup_n_f32(0.0f), a.simdInternal_) };
}

static inline SimdFloat gmx_simdcall blend(SimdFloat a, SimdFloat b, SimdFBool sel)
{
    return { svsel_f32(sel.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline SimdFInt32 gmx_simdcall operator&(SimdFInt32 a, SimdFInt32 b)
{
    svbool_t pg = svptrue_b32();
    return { svand_s32_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFInt32 gmx_simdcall andNot(SimdFInt32 a, SimdFInt32 b)
{
    svbool_t pg = svptrue_b32();
    return { svbic_s32_x(pg, b.simdInternal_, a.simdInternal_) };
}

static inline SimdFInt32 gmx_simdcall operator|(SimdFInt32 a, SimdFInt32 b)
{
    svbool_t pg = svptrue_b32();
    return { svorr_s32_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFInt32 gmx_simdcall operator^(SimdFInt32 a, SimdFInt32 b)
{
    svbool_t pg = svptrue_b32();
    return { sveor_s32_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFInt32 gmx_simdcall operator+(SimdFInt32 a, SimdFInt32 b)
{
    svbool_t pg = svptrue_b32();
    return { svadd_s32_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFInt32 gmx_simdcall operator-(SimdFInt32 a, SimdFInt32 b)
{
    svbool_t pg = svptrue_b32();
    return { svsub_s32_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFInt32 gmx_simdcall operator*(SimdFInt32 a, SimdFInt32 b)
{
    svbool_t pg = svptrue_b32();
    return { svmul_s32_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFIBool gmx_simdcall operator==(SimdFInt32 a, SimdFInt32 b)
{
    svbool_t pg = svptrue_b32();
    return { svcmpeq_s32(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFIBool gmx_simdcall testBits(SimdFInt32 a)
{
    svbool_t pg = svptrue_b32();
    return { svcmpne_n_s32(pg, a.simdInternal_, (int32_t)0) };
}

static inline SimdFIBool gmx_simdcall operator<(SimdFInt32 a, SimdFInt32 b)
{
    svbool_t pg = svptrue_b32();
    return { svcmplt_s32(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFIBool gmx_simdcall operator&&(SimdFIBool a, SimdFIBool b)
{
    svbool_t pg = svptrue_b32();
    return { svand_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdFIBool gmx_simdcall operator||(SimdFIBool a, SimdFIBool b)
{
    svbool_t pg = svptrue_b32();
    return { svorr_b_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline bool gmx_simdcall anyTrue(SimdFIBool a)
{
    svbool_t pg = svptrue_b32();
    return svptest_any(pg, a.simdInternal_);
}

static inline SimdFInt32 gmx_simdcall selectByMask(SimdFInt32 a, SimdFIBool m)
{
    return { svsel_s32(m.simdInternal_, a.simdInternal_, svdup_n_s32(0)) };
}

static inline SimdFInt32 gmx_simdcall selectByNotMask(SimdFInt32 a, SimdFIBool m)
{
    return { svsel_s32(m.simdInternal_, svdup_n_s32(0), a.simdInternal_) };
}

static inline SimdFInt32 gmx_simdcall blend(SimdFInt32 a, SimdFInt32 b, SimdFIBool sel)
{
    return { svsel_s32(sel.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline SimdFInt32 gmx_simdcall cvtR2I(SimdFloat a)
{
    svbool_t pg = svptrue_b32();
    return { svcvt_s32_x(pg, svrinta_f32_x(pg, a.simdInternal_)) };
}

static inline SimdFInt32 gmx_simdcall cvttR2I(SimdFloat a)
{
    svbool_t pg = svptrue_b32();
    return { svcvt_s32_x(pg, a.simdInternal_) };
}

static inline SimdFloat gmx_simdcall cvtI2R(SimdFInt32 a)
{
    svbool_t pg = svptrue_b32();
    return { svcvt_f32_x(pg, a.simdInternal_) };
}

static inline SimdFIBool gmx_simdcall cvtB2IB(SimdFBool a)
{
    return { a.simdInternal_ };
}

static inline SimdFBool gmx_simdcall cvtIB2B(SimdFIBool a)
{
    return { a.simdInternal_ };
}

static inline SimdFloat gmx_simdcall round(SimdFloat x)
{
    svbool_t pg = svptrue_b32();
    return { svrinta_f32_x(pg, x.simdInternal_) };
}

static inline SimdFloat gmx_simdcall trunc(SimdFloat x)
{
    return cvtI2R(cvttR2I(x));
}

} // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_SVE_SIMD_FLOAT_H
