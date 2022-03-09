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

#ifndef GMX_SIMD_IMPL_ARM_SVE_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_ARM_SVE_SIMD_DOUBLE_H

#include "config.h"

#include <arm_sve.h>

#include <cassert>
#include <cstddef>
#include <cstdint>

#include "gromacs/math/utilities.h"

#include "impl_arm_sve_simd_float.h"

namespace gmx
{

class SimdDouble
{
private:
    typedef svfloat64_t simdInternalType_
            __attribute__((arm_sve_vector_bits(GMX_SIMD_ARM_SVE_LENGTH_VALUE)));

public:
    SimdDouble() {}

    SimdDouble(const double d) { this->simdInternal_ = svdup_n_f64(d); }

    SimdDouble(svfloat64_t simd) : simdInternal_(simd) {}

    simdInternalType_ simdInternal_;
};

class SimdDInt32
{
private:
    typedef svint64_t simdInternalType_
            __attribute__((arm_sve_vector_bits(GMX_SIMD_ARM_SVE_LENGTH_VALUE)));

public:
    SimdDInt32() {}

    SimdDInt32(const int32_t i) { this->simdInternal_ = svdup_n_s64(i); }

    SimdDInt32(svint64_t simd) : simdInternal_(simd) {}

    simdInternalType_ simdInternal_;
};

class SimdDBool
{
private:
    typedef svbool_t simdInternalType_
            __attribute__((arm_sve_vector_bits(GMX_SIMD_ARM_SVE_LENGTH_VALUE)));

public:
    SimdDBool() {}

    SimdDBool(const bool b) { this->simdInternal_ = svdup_n_b64(b); }

    SimdDBool(svbool_t simd) : simdInternal_(simd) {}

    simdInternalType_ simdInternal_;
};

class SimdDIBool
{
private:
    typedef svbool_t simdInternalType_
            __attribute__((arm_sve_vector_bits(GMX_SIMD_ARM_SVE_LENGTH_VALUE)));

public:
    SimdDIBool() {}

    SimdDIBool(const bool b) { this->simdInternal_ = svdup_n_b64(b); }

    SimdDIBool(svbool_t simd) : simdInternal_(simd) {}

    simdInternalType_ simdInternal_;
};

static inline SimdDouble gmx_simdcall simdLoad(const double* m, SimdDoubleTag = {})
{
    assert(0 == (std::size_t(m) % GMX_SIMD_ALIGNMENT));
    svbool_t pg = svptrue_b64();
    return { svld1_f64(pg, m) };
}

static inline SimdDouble gmx_simdcall simdLoad(SimdDouble* m, int offset, SimdDoubleTag = {})
{
    assert(0 == (std::size_t(m) % GMX_SIMD_ALIGNMENT));
    svbool_t pg = svptrue_b64();
    return { svld1_f64(pg, reinterpret_cast<double*>(m) + offset * svcntd()) };
}

static inline SimdDouble gmx_simdcall simdLoadDouble(const double* m)
{
    assert(0 == (std::size_t(m) % GMX_SIMD_ALIGNMENT));
    svbool_t pg = svptrue_b64();
    return { svld1_f64(pg, m) };
}

static inline void gmx_simdcall store(double* m, SimdDouble a)
{
    assert(0 == (std::size_t(m) % GMX_SIMD_ALIGNMENT));
    svbool_t pg = svptrue_b64();
    svst1_f64(pg, m, a.simdInternal_);
}

static inline SimdDouble gmx_simdcall simdLoadU(const double* m, SimdDoubleTag = {})
{
    svbool_t pg = svptrue_b64();
    return { svld1_f64(pg, m) };
}

static inline void gmx_simdcall storeU(double* m, SimdDouble a)
{
    svbool_t pg = svptrue_b64();
    svst1_f64(pg, m, a.simdInternal_);
}

static inline SimdDouble gmx_simdcall setZeroD()
{
    return { svdup_n_f64(0.0) };
}

static inline SimdDInt32 gmx_simdcall simdLoad(const std::int32_t* m, SimdDInt32Tag)
{
    assert(0 == (std::size_t(m) % GMX_SIMD_ALIGNMENT));
    svbool_t pg = SVE_SIMD_FLOAT_HALF_DOUBLE_MASK;
    return { svunpklo_s64(svld1_s32(pg, m)) };
}

static inline void gmx_simdcall store(std::int32_t* m, SimdDInt32 a)
{
    assert(0 == (std::size_t(m) % GMX_SIMD_ALIGNMENT));
    svbool_t pg = SVE_SIMD_FLOAT_HALF_DOUBLE_MASK;
    svst1_s32(pg, m, svuzp1(svreinterpret_s32_s64(a.simdInternal_), svreinterpret_s32_s64(a.simdInternal_)));
}

static inline SimdDInt32 gmx_simdcall simdLoadU(const std::int32_t* m, SimdDInt32Tag)
{
    svbool_t pg = SVE_SIMD_FLOAT_HALF_DOUBLE_MASK;
    return { svunpklo_s64(svld1_s32(pg, m)) };
}

static inline void gmx_simdcall storeU(std::int32_t* m, SimdDInt32 a)
{
    svbool_t pg = SVE_SIMD_FLOAT_HALF_DOUBLE_MASK;
    svst1_s32(pg, m, svuzp1(svreinterpret_s32_s64(a.simdInternal_), svreinterpret_s32_s64(a.simdInternal_)));
}

static inline SimdDInt32 gmx_simdcall setZeroDI()
{
    return { svdup_n_s64(0) };
}

template<int index>
gmx_simdcall static inline std::int32_t extract(SimdDInt32 a)
{
    svbool_t pg = svwhilelt_b64(0, index);
    return svlasta_s64(pg, a.simdInternal_);
}

template<int index>
gmx_simdcall static inline double extract(SimdDouble a)
{
    svbool_t pg = svwhilelt_b64(0, index);
    return svlasta_f64(pg, a.simdInternal_);
}

static inline SimdDouble gmx_simdcall operator&(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svreinterpret_f64_s64(svand_s64_x(
            pg, svreinterpret_s64_f64(a.simdInternal_), svreinterpret_s64_f64(b.simdInternal_))) };
}

static inline SimdDouble gmx_simdcall andNot(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svreinterpret_f64_s64(svbic_s64_x(
            pg, svreinterpret_s64_f64(b.simdInternal_), svreinterpret_s64_f64(a.simdInternal_))) };
}

static inline SimdDouble gmx_simdcall operator|(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svreinterpret_f64_s64(svorr_s64_x(
            pg, svreinterpret_s64_f64(a.simdInternal_), svreinterpret_s64_f64(b.simdInternal_))) };
}

static inline SimdDouble gmx_simdcall operator^(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svreinterpret_f64_s64(sveor_s64_x(
            pg, svreinterpret_s64_f64(a.simdInternal_), svreinterpret_s64_f64(b.simdInternal_))) };
}

static inline SimdDouble gmx_simdcall operator+(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svadd_f64_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDouble gmx_simdcall operator-(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svsub_f64_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDouble gmx_simdcall operator-(SimdDouble a)
{
    svbool_t pg = svptrue_b64();
    return { svneg_f64_x(pg, a.simdInternal_) };
}

static inline SimdDouble gmx_simdcall operator*(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svmul_f64_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDouble gmx_simdcall fma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    svbool_t pg = svptrue_b64();
    return { svmad_f64_x(pg, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline SimdDouble gmx_simdcall fms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    svbool_t pg = svptrue_b64();
    return { svnmsb_f64_x(pg, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline SimdDouble gmx_simdcall fnma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    svbool_t pg = svptrue_b64();
    return { svmsb_f64_x(pg, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline SimdDouble gmx_simdcall fnms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    svbool_t pg = svptrue_b64();
    return { svnmad_f64_x(pg, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline SimdDouble gmx_simdcall rsqrt(SimdDouble x)
{
    return { svrsqrte_f64(x.simdInternal_) };
}

// The SIMD implementation seems to overflow when we square lu for
// values close to FLOAT_MAX, so we fall back on the version in
// simd_math.h, which is probably slightly slower.
#if GMX_SIMD_HAVE_NATIVE_RSQRT_ITER_DOUBLE
static inline SimdDouble gmx_simdcall rsqrtIter(SimdDouble lu, SimdDouble x)
{
    return { vmulq_f64(lu.simdInternal_,
                       vrsqrtsq_f32(vmulq_f32(lu.simdInternal_, lu.simdInternal_), x.simdInternal_)) };
}
#endif

static inline SimdDouble gmx_simdcall rcp(SimdDouble x)
{
    return { svrecpe_f64(x.simdInternal_) };
}

static inline SimdDouble gmx_simdcall rcpIter(SimdDouble lu, SimdDouble x)
{
    svbool_t pg = svptrue_b64();
    return { svmul_f64_x(pg, lu.simdInternal_, svrecps_f64(lu.simdInternal_, x.simdInternal_)) };
}

static inline SimdDouble gmx_simdcall maskAdd(SimdDouble a, SimdDouble b, SimdDBool m)
{
    return { svadd_f64_m(m.simdInternal_, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDouble gmx_simdcall maskzMul(SimdDouble a, SimdDouble b, SimdDBool m)
{
    return { svmul_f64_z(m.simdInternal_, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDouble gmx_simdcall maskzFma(SimdDouble a, SimdDouble b, SimdDouble c, SimdDBool m)
{
    return { svmad_f64_z(m.simdInternal_, a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline SimdDouble gmx_simdcall maskzRsqrt(SimdDouble x, SimdDBool m)
{
    // The result will always be correct since we mask the result with m, but
    // for debug builds we also want to make sure not to generate FP exceptions
#ifndef NDEBUG
    x.simdInternal_ = svsel_f64(m.simdInternal_, x.simdInternal_, svdup_n_f64(1.0));
#endif
    return { svsel_f64(m.simdInternal_, svrsqrte_f64(x.simdInternal_), svdup_n_f64(0.0)) };
}

static inline SimdDouble gmx_simdcall maskzRcp(SimdDouble x, SimdDBool m)
{
    // The result will always be correct since we mask the result with m, but
    // for debug builds we also want to make sure not to generate FP exceptions
#ifndef NDEBUG
    x.simdInternal_ = svsel_f64(m.simdInternal_, x.simdInternal_, svdup_n_f64(1.0));
#endif
    return { svsel_f64(m.simdInternal_, svrecpe_f64(x.simdInternal_), svdup_n_f64(0.0)) };
}

static inline SimdDouble gmx_simdcall abs(SimdDouble x)
{
    svbool_t pg = svptrue_b64();
    return { svabs_f64_x(pg, x.simdInternal_) };
}

static inline SimdDouble gmx_simdcall max(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svmax_f64_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDouble gmx_simdcall min(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svmin_f64_x(pg, a.simdInternal_, b.simdInternal_) };
}

// Round and trunc operations are defined at the end of this file, since they
// need to use double-to-integer and integer-to-double conversions.

template<MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble gmx_simdcall frexp(SimdDouble value, SimdDInt32* exponent)
{
    svbool_t        pg           = svptrue_b64();
    const svint64_t exponentMask = svdup_n_s64(0x7FF0000000000000LL);
    const svint64_t mantissaMask = svdup_n_s64(0x800FFFFFFFFFFFFFLL);
    const svint64_t exponentBias = svdup_n_s64(1022LL); // add 1 to make our definition identical to frexp()
    const svfloat64_t half = svdup_n_f64(0.5);
    svint64_t         iExponent;

    iExponent = svand_s64_x(pg, svreinterpret_s64_f64(value.simdInternal_), exponentMask);
    iExponent = svsub_s64_x(
            pg, svreinterpret_s64_u64(svlsr_n_u64_x(pg, svreinterpret_u64_s64(iExponent), 52)), exponentBias);


    svfloat64_t result = svreinterpret_f64_s64(
            svorr_s64_x(pg,
                        svand_s64_x(pg, svreinterpret_s64_f64(value.simdInternal_), mantissaMask),
                        svreinterpret_s64_f64(half)));

    if (opt == MathOptimization::Safe)
    {
        svbool_t valueIsZero = svcmpeq_n_f64(pg, value.simdInternal_, 0.0);
        iExponent            = svsel_s64(valueIsZero, svdup_n_s64(0), iExponent);
        result               = svsel_f64(valueIsZero, value.simdInternal_, result);
    }

    exponent->simdInternal_ = iExponent;
    return { result };
}

template<MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble gmx_simdcall ldexp(SimdDouble value, SimdDInt32 exponent)
{
    svbool_t        pg           = svptrue_b64();
    const svint64_t exponentBias = svdup_n_s64(1023);
    svint64_t       iExponent    = svadd_s64_x(pg, exponent.simdInternal_, exponentBias);

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponent = svmax_n_s64_x(pg, iExponent, 0);
    }

    iExponent = svlsl_n_s64_x(pg, iExponent, 52);

    return { svmul_f64_x(pg, value.simdInternal_, svreinterpret_f64_s64(iExponent)) };
}

static inline double gmx_simdcall reduce(SimdDouble a)
{
    svbool_t pg = svptrue_b64();
    return svadda_f64(pg, 0.0, a.simdInternal_);
}

static inline SimdDBool gmx_simdcall operator==(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svcmpeq_f64(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDBool gmx_simdcall operator!=(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svcmpne_f64(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDBool gmx_simdcall operator<(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svcmplt_f64(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDBool gmx_simdcall operator<=(SimdDouble a, SimdDouble b)
{
    svbool_t pg = svptrue_b64();
    return { svcmple_f64(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDBool gmx_simdcall testBits(SimdDouble a)
{
    svbool_t pg = svptrue_b64();
    return { svcmpne_n_s64(pg, svreinterpret_s64_f64(a.simdInternal_), 0) };
}

static inline SimdDBool gmx_simdcall operator&&(SimdDBool a, SimdDBool b)
{
    svbool_t pg = svptrue_b64();
    return { svand_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDBool gmx_simdcall operator||(SimdDBool a, SimdDBool b)
{
    svbool_t pg = svptrue_b64();
    return { svorr_b_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline bool gmx_simdcall anyTrue(SimdDBool a)
{
    svbool_t pg = svptrue_b64();
    return svptest_any(pg, a.simdInternal_);
}

static inline bool gmx_simdcall extractFirst(SimdDBool a)
{
    svbool_t pg = svptrue_b64();
    return svptest_first(pg, a.simdInternal_);
}

static inline SimdDouble gmx_simdcall selectByMask(SimdDouble a, SimdDBool m)
{
    return { svsel_f64(m.simdInternal_, a.simdInternal_, svdup_n_f64(0.0)) };
}

static inline SimdDouble gmx_simdcall selectByNotMask(SimdDouble a, SimdDBool m)
{
    return { svsel_f64(m.simdInternal_, svdup_n_f64(0.0), a.simdInternal_) };
}

static inline SimdDouble gmx_simdcall blend(SimdDouble a, SimdDouble b, SimdDBool sel)
{
    return { svsel_f64(sel.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall operator&(SimdDInt32 a, SimdDInt32 b)
{
    svbool_t pg = svptrue_b64();
    return { svand_s64_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall andNot(SimdDInt32 a, SimdDInt32 b)
{
    svbool_t pg = svptrue_b64();
    return { svbic_s64_x(pg, b.simdInternal_, a.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall operator|(SimdDInt32 a, SimdDInt32 b)
{
    svbool_t pg = svptrue_b64();
    return { svorr_s64_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall operator^(SimdDInt32 a, SimdDInt32 b)
{
    svbool_t pg = svptrue_b64();
    return { sveor_s64_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall operator+(SimdDInt32 a, SimdDInt32 b)
{
    svbool_t pg = svptrue_b64();
    return { svadd_s64_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall operator-(SimdDInt32 a, SimdDInt32 b)
{
    svbool_t pg = svptrue_b64();
    return { svsub_s64_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall operator*(SimdDInt32 a, SimdDInt32 b)
{
    svbool_t pg = svptrue_b64();
    return { svmul_s64_x(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDIBool gmx_simdcall operator==(SimdDInt32 a, SimdDInt32 b)
{
    svbool_t pg = svptrue_b64();
    return { svcmpeq_s64(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDIBool gmx_simdcall testBits(SimdDInt32 a)
{
    svbool_t pg = svptrue_b64();
    return { svcmpne_n_s64(pg, a.simdInternal_, (int64_t)0) };
}

static inline SimdDIBool gmx_simdcall operator<(SimdDInt32 a, SimdDInt32 b)
{
    svbool_t pg = svptrue_b64();
    return { svcmplt_s64(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDIBool gmx_simdcall operator&&(SimdDIBool a, SimdDIBool b)
{
    svbool_t pg = svptrue_b64();
    return { svand_b_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline SimdDIBool gmx_simdcall operator||(SimdDIBool a, SimdDIBool b)
{
    svbool_t pg = svptrue_b64();
    return { svorr_b_z(pg, a.simdInternal_, b.simdInternal_) };
}

static inline bool gmx_simdcall anyTrue(SimdDIBool a)
{
    svbool_t pg = svptrue_b64();
    return svptest_any(pg, a.simdInternal_);
}

static inline SimdDInt32 gmx_simdcall selectByMask(SimdDInt32 a, SimdDIBool m)
{
    return { svsel_s64(m.simdInternal_, a.simdInternal_, svdup_n_s64(0)) };
}

static inline SimdDInt32 gmx_simdcall selectByNotMask(SimdDInt32 a, SimdDIBool m)
{
    return { svsel_s64(m.simdInternal_, svdup_n_s64(0), a.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall blend(SimdDInt32 a, SimdDInt32 b, SimdDIBool sel)
{
    return { svsel_s64(sel.simdInternal_, b.simdInternal_, a.simdInternal_) };
}

static inline SimdDInt32 gmx_simdcall cvtR2I(SimdDouble a)
{
    svbool_t pg = svptrue_b64();
    return { svcvt_s64_x(pg, svrinta_f64_x(pg, a.simdInternal_)) };
}

static inline SimdDInt32 gmx_simdcall cvttR2I(SimdDouble a)
{
    // FIXME ???
    svbool_t pg = svptrue_b64();
    return { svcvt_s64_x(pg, a.simdInternal_) };
}

static inline SimdDouble gmx_simdcall cvtI2R(SimdDInt32 a)
{
    svbool_t pg = svptrue_b64();
    return { svcvt_f64_x(pg, a.simdInternal_) };
}

static inline SimdDIBool gmx_simdcall cvtB2IB(SimdDBool a)
{
    return { a.simdInternal_ };
}

static inline SimdDBool gmx_simdcall cvtIB2B(SimdDIBool a)
{
    return { a.simdInternal_ };
}

static inline SimdDouble gmx_simdcall round(SimdDouble x)
{
    svbool_t pg = svptrue_b64();
    return { svrinta_f64_x(pg, x.simdInternal_) };
}

static inline SimdDouble gmx_simdcall trunc(SimdDouble x)
{
    return cvtI2R(cvttR2I(x));
}

static inline void gmx_simdcall cvtF2DD(SimdFloat gmx_unused f,
                                        SimdDouble gmx_unused* d0,
                                        SimdDouble gmx_unused* d1)
{
    assert(GMX_SIMD_FLOAT_WIDTH == 2 * GMX_SIMD_DOUBLE_WIDTH);
    svbool_t pg       = svptrue_b32();
    d0->simdInternal_ = svcvt_f64_f32_x(pg, svzip1(f.simdInternal_, f.simdInternal_));
    d1->simdInternal_ = svcvt_f64_f32_x(pg, svzip2(f.simdInternal_, f.simdInternal_));
}

static inline SimdFloat gmx_simdcall cvtDD2F(SimdDouble gmx_unused d0, SimdDouble gmx_unused d1)
{
    svbool_t pg = svptrue_b64();
    assert(GMX_SIMD_FLOAT_WIDTH == 2 * GMX_SIMD_DOUBLE_WIDTH);
    return { svuzp1_f32(svcvt_f32_f64_x(pg, d0.simdInternal_), svcvt_f32_f64_x(pg, d1.simdInternal_)) };
}

} // namespace gmx

#endif // GMX_SIMD_IMPL_ARM_SVE_SIMD_DOUBLE_H
