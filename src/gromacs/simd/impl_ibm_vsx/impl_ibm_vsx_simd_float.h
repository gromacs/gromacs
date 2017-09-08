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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD_FLOAT_H

#include "config.h"

#include "gromacs/math/utilities.h"
#include "gromacs/utility/basedefinitions.h"

#include "impl_ibm_vsx_definitions.h"

namespace gmx
{

class SimdFloat
{
    public:
        SimdFloat() {}

        // gcc-4.9 does not recognize that we use the parameter
        SimdFloat(float gmx_unused f) : simdInternal_(vec_splats(f)) {}

        // Internal utility constructor to simplify return statements
        SimdFloat(__vector float simd) : simdInternal_(simd) {}

        __vector float  simdInternal_;
};

class SimdFInt32
{
    public:
        SimdFInt32() {}

        // gcc-4.9 does not recognize that we use the parameter
        SimdFInt32(std::int32_t gmx_unused i) : simdInternal_(vec_splats(i)) {}

        // Internal utility constructor to simplify return statements
        SimdFInt32(__vector signed int simd) : simdInternal_(simd) {}

        __vector signed int  simdInternal_;
};

class SimdFBool
{
    public:
        SimdFBool() {}

        SimdFBool(bool b) : simdInternal_(reinterpret_cast<__vector vsxBool int>(vec_splats( b ? 0xFFFFFFFF : 0))) {}

        // Internal utility constructor to simplify return statements
        SimdFBool(__vector vsxBool int simd) : simdInternal_(simd) {}

        __vector vsxBool int  simdInternal_;
};

class SimdFIBool
{
    public:
        SimdFIBool() {}

        SimdFIBool(bool b) : simdInternal_(reinterpret_cast<__vector vsxBool int>(vec_splats( b ? 0xFFFFFFFF : 0))) {}

        // Internal utility constructor to simplify return statements
        SimdFIBool(__vector vsxBool int simd) : simdInternal_(simd) {}

        __vector vsxBool int  simdInternal_;
};

// Note that the interfaces we use here have been a mess in xlc;
// currently version 13.1.5 is required.

static inline SimdFloat gmx_simdcall
simdLoad(const float *m)
{
    return {
               *reinterpret_cast<const __vector float *>(m)
    };
}

static inline void gmx_simdcall
store(float *m, SimdFloat a)
{
    *reinterpret_cast<__vector float *>(m) = a.simdInternal_;
}

static inline SimdFloat gmx_simdcall
simdLoadU(const float *m)
{
    return {
               *reinterpret_cast<const __vector float *>(m)
    };
}

static inline void gmx_simdcall
storeU(float *m, SimdFloat a)
{
    *reinterpret_cast<__vector float *>(m) = a.simdInternal_;
}

static inline SimdFloat gmx_simdcall
setZeroF()
{
    return {
               vec_splats(0.0f)
    };
}

static inline SimdFInt32 gmx_simdcall
simdLoadFI(const std::int32_t * m)
{
    return {
               *reinterpret_cast<const __vector int *>(m)
    };
}

static inline void gmx_simdcall
store(std::int32_t * m, SimdFInt32 a)
{
    *reinterpret_cast<__vector int *>(m) = a.simdInternal_;
}

static inline SimdFInt32 gmx_simdcall
simdLoadUFI(const std::int32_t *m)
{
    return {
               *reinterpret_cast<const __vector int *>(m)
    };
}

static inline void gmx_simdcall
storeU(std::int32_t * m, SimdFInt32 a)
{
    *reinterpret_cast<__vector int *>(m) = a.simdInternal_;
}

static inline SimdFInt32 gmx_simdcall
setZeroFI()
{
    return {
               vec_splats(static_cast<int>(0))
    };
}

// gcc-4.9 does not detect that vec_extract() uses its argument
template<int index>
static inline std::int32_t gmx_simdcall
extract(SimdFInt32 gmx_unused a)
{
    return vec_extract(a.simdInternal_, index);
}

static inline SimdFloat gmx_simdcall
operator&(SimdFloat a, SimdFloat b)
{
    return {
               vec_and(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
andNot(SimdFloat a, SimdFloat b)
{
    return {
               vec_andc(b.simdInternal_, a.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator|(SimdFloat a, SimdFloat b)
{
    return {
               vec_or(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator^(SimdFloat a, SimdFloat b)
{
    return {
               vec_xor(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator+(SimdFloat a, SimdFloat b)
{
    return {
               vec_add(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator-(SimdFloat a, SimdFloat b)
{
    return {
               vec_sub(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
operator-(SimdFloat x)
{
    return {
               -x.simdInternal_
    };
}

static inline SimdFloat gmx_simdcall
operator*(SimdFloat a, SimdFloat b)
{
    return {
               vec_mul(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               vec_madd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               vec_msub(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fnma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               vec_nmsub(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
fnms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return {
               vec_nmadd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
rsqrt(SimdFloat x)
{
    return {
               vec_rsqrte(x.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
rcp(SimdFloat x)
{
    return {
               vec_re(x.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
maskAdd(SimdFloat a, SimdFloat b, SimdFBool m)
{
    return {
               vec_add(a.simdInternal_, vec_and(b.simdInternal_, reinterpret_cast<__vector float>(m.simdInternal_)))
    };
}

static inline SimdFloat gmx_simdcall
maskzMul(SimdFloat a, SimdFloat b, SimdFBool m)
{
    SimdFloat prod = a * b;

    return {
               vec_and(prod.simdInternal_, reinterpret_cast<__vector float>(m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
maskzFma(SimdFloat a, SimdFloat b, SimdFloat c, SimdFBool m)
{
    SimdFloat prod = fma(a, b, c);

    return {
               vec_and(prod.simdInternal_, reinterpret_cast<__vector float>(m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
maskzRsqrt(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = vec_sel(vec_splats(1.0f), x.simdInternal_, m.simdInternal_);
#endif
    return {
               vec_and(vec_rsqrte(x.simdInternal_), reinterpret_cast<__vector float>(m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
maskzRcp(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = vec_sel(vec_splats(1.0f), x.simdInternal_, m.simdInternal_);
#endif
    return {
               vec_and(vec_re(x.simdInternal_), reinterpret_cast<__vector float>(m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
abs(SimdFloat x)
{
    return {
               vec_abs( x.simdInternal_ )
    };
}

static inline SimdFloat gmx_simdcall
max(SimdFloat a, SimdFloat b)
{
    return {
               vec_max(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
min(SimdFloat a, SimdFloat b)
{
    return {
               vec_min(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
round(SimdFloat x)
{
    return {
               vec_round( x.simdInternal_ )
    };
}

static inline SimdFloat gmx_simdcall
trunc(SimdFloat x)
{
    return {
               vec_trunc( x.simdInternal_ )
    };
}

static inline SimdFloat gmx_simdcall
frexp(SimdFloat value, SimdFInt32 * exponent)
{
    const __vector float      exponentMask   = reinterpret_cast<__vector float>(vec_splats(0x7F800000U));
    const __vector signed int exponentBias   = vec_splats(126);
    const __vector float      half           = vec_splats(0.5f);
    __vector signed int       iExponent;

    iExponent               = reinterpret_cast<__vector signed int>(vec_and(value.simdInternal_, exponentMask));
    iExponent               = vec_sub(vec_sr(iExponent, vec_splats(23U)), exponentBias);
    exponent->simdInternal_ = iExponent;

    return {
               vec_or( vec_andc(value.simdInternal_, exponentMask), half)
    };
}

template <MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall
ldexp(SimdFloat value, SimdFInt32 exponent)
{
    const __vector signed int exponentBias   = vec_splats(127);
    __vector signed int       iExponent;

    iExponent  = vec_add(exponent.simdInternal_, exponentBias);

    if (opt == MathOptimization::Safe)
    {
        // Make sure biased argument is not negative
        iExponent = vec_max(iExponent, vec_splat_s32(0));
    }

    iExponent = vec_sl( iExponent, vec_splats(23U));

    return {
               vec_mul(value.simdInternal_, reinterpret_cast<__vector float>(iExponent))
    };
}

static inline float gmx_simdcall
reduce(SimdFloat x)
{
    const __vector unsigned char perm1 = { 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7 };
    const __vector unsigned char perm2 = { 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3 };

    x.simdInternal_ = vec_add(x.simdInternal_, vec_perm(x.simdInternal_, x.simdInternal_, perm1));
    x.simdInternal_ = vec_add(x.simdInternal_, vec_perm(x.simdInternal_, x.simdInternal_, perm2));
    return vec_extract(x.simdInternal_, 0);
}

static inline SimdFBool gmx_simdcall
operator==(SimdFloat a, SimdFloat b)
{
    return {
               vec_cmpeq(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
operator!=(SimdFloat a, SimdFloat b)
{
    return {
               vec_or(vec_cmpgt(a.simdInternal_, b.simdInternal_),
                      vec_cmplt(a.simdInternal_, b.simdInternal_))
    };
}

static inline SimdFBool gmx_simdcall
operator<(SimdFloat a, SimdFloat b)
{
    return {
               vec_cmplt(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
operator<=(SimdFloat a, SimdFloat b)
{
    return {
               vec_cmple(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
testBits(SimdFloat a)
{
    return {
               vec_cmpgt( reinterpret_cast<__vector unsigned int>(a.simdInternal_), vec_splats(0U))
    };
}

static inline SimdFBool gmx_simdcall
operator&&(SimdFBool a, SimdFBool b)
{
    return {
               vec_and(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFBool gmx_simdcall
operator||(SimdFBool a, SimdFBool b)
{
    return {
               vec_or(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdFBool a)
{
    return vec_any_ne(a.simdInternal_, reinterpret_cast<__vector vsxBool int>(vec_splats(0)));
}

static inline SimdFloat gmx_simdcall
selectByMask(SimdFloat a, SimdFBool m)
{
    return {
               vec_and(a.simdInternal_, reinterpret_cast<__vector float>(m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
selectByNotMask(SimdFloat a, SimdFBool m)
{
    return {
               vec_andc(a.simdInternal_, reinterpret_cast<__vector float>(m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
blend(SimdFloat a, SimdFloat b, SimdFBool sel)
{
    return {
               vec_sel(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator<<(SimdFInt32 a, int n)
{
    return {
               vec_sl(a.simdInternal_, vec_splats(static_cast<unsigned int>(n)))
    };
}

static inline SimdFInt32 gmx_simdcall
operator>>(SimdFInt32 a, int n)
{
    return {
               vec_sr(a.simdInternal_, vec_splats(static_cast<unsigned int>(n)))
    };
}

static inline SimdFInt32 gmx_simdcall
operator&(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vec_and(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
andNot(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vec_andc(b.simdInternal_, a.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator|(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vec_or(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator^(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vec_xor(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator+(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vec_add(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator-(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vec_sub(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
operator*(SimdFInt32 a, SimdFInt32 b)
{
    return {
               a.simdInternal_ * b.simdInternal_
    };
}

static inline SimdFIBool gmx_simdcall
operator==(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vec_cmpeq(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
testBits(SimdFInt32 a)
{
    return {
               vec_cmpgt( reinterpret_cast<__vector unsigned int>(a.simdInternal_), vec_splats(0U))
    };
}

static inline SimdFIBool gmx_simdcall
operator<(SimdFInt32 a, SimdFInt32 b)
{
    return {
               vec_cmplt(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator&&(SimdFIBool a, SimdFIBool b)
{
    return {
               vec_and(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdFIBool gmx_simdcall
operator||(SimdFIBool a, SimdFIBool b)
{
    return {
               vec_or(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdFIBool a)
{
    return vec_any_ne(a.simdInternal_, reinterpret_cast<__vector vsxBool int>(vec_splats(0)));
}

static inline SimdFInt32 gmx_simdcall
selectByMask(SimdFInt32 a, SimdFIBool m)
{
    return {
               vec_and(a.simdInternal_, reinterpret_cast<__vector signed int>(m.simdInternal_))
    };
}

static inline SimdFInt32 gmx_simdcall
selectByNotMask(SimdFInt32 a, SimdFIBool m)
{
    return {
               vec_andc(a.simdInternal_, reinterpret_cast<__vector signed int>(m.simdInternal_))
    };
}

static inline SimdFInt32 gmx_simdcall
blend(SimdFInt32 a, SimdFInt32 b, SimdFIBool sel)
{
    return {
               vec_sel(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
cvtR2I(SimdFloat a)
{
    return {
               vec_cts(vec_round(a.simdInternal_), 0)
    };
}

static inline SimdFInt32 gmx_simdcall
cvttR2I(SimdFloat a)
{
    return {
               vec_cts(a.simdInternal_, 0)
    };
}

static inline SimdFloat gmx_simdcall
cvtI2R(SimdFInt32 a)
{
    return {
               vec_ctf(a.simdInternal_, 0)
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

static inline SimdFloat gmx_simdcall
copysign(SimdFloat x, SimdFloat y)
{
#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__)
    __vector float res;
    __asm__ ("xvcpsgnsp %x0,%x1,%x2" : "=wf" (res) : "wf" (y.simdInternal_), "wf" (x.simdInternal_));
    return {
               res
    };
#else
    return {
               vec_cpsgn(y.simdInternal_, x.simdInternal_)
    };
#endif
}

}      // namespace gmx

#endif // GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD_FLOAT_H
