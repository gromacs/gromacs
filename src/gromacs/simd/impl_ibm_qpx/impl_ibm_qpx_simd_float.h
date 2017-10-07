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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_FLOAT_H

#include "config.h"

// Assert is buggy on xlc with high optimization, so we skip it for QPX
#include <cstddef>
#include <cstdint>

#ifdef __clang__
#include <qpxmath.h>
#endif

#include "gromacs/math/utilities.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

class SimdFloat
{
    public:
        SimdFloat() {}

        SimdFloat(float f) : simdInternal_(vec_splats(f)) {}

        // Internal utility constructor to simplify return statements
        SimdFloat(vector4double simd) : simdInternal_(simd) {}

        vector4double  simdInternal_;
};

class SimdFInt32
{
    public:
        SimdFInt32() {}

        SimdFInt32(std::int32_t i)
        {
            GMX_ALIGNED(int, GMX_SIMD_FINT32_WIDTH) idata[GMX_SIMD_FINT32_WIDTH];
            idata[0]      = i;
            simdInternal_ = vec_splat(vec_ldia(0, idata), 0);
        }

        // Internal utility constructor to simplify return statements
        SimdFInt32(vector4double simd) : simdInternal_(simd) {}

        vector4double  simdInternal_;
};

class SimdFBool
{
    public:
        SimdFBool() {}

        SimdFBool(bool b) : simdInternal_(vec_splats(b ? 1.0 : -1.0)) {}

        // Internal utility constructor to simplify return statements
        SimdFBool(vector4double simd) : simdInternal_(simd) {}

        vector4double  simdInternal_;
};

static inline SimdFloat gmx_simdcall
simdLoad(const float *m, SimdFloatTag = {})
{
#ifdef NDEBUG
    return {
               vec_ld(0, const_cast<float *>(m))
    };
#else
    return {
               vec_lda(0, const_cast<float *>(m))
    };
#endif
}

static inline void gmx_simdcall
store(float *m, SimdFloat a)
{
#ifdef NDEBUG
    vec_st(a.simdInternal_, 0, m);
#else
    vec_sta(a.simdInternal_, 0, m);
#endif
}

static inline SimdFloat gmx_simdcall
setZeroF()
{
    return {
               vec_splats(0.0)
    };
}

static inline SimdFInt32 gmx_simdcall
simdLoad(const std::int32_t * m, SimdFInt32Tag)
{
#ifdef NDEBUG
    return {
               vec_ldia(0, const_cast<int *>(m))
    };
#else
    return {
               vec_ldiaa(0, const_cast<int *>(m))
    };
#endif
}

static inline void gmx_simdcall
store(std::int32_t * m, SimdFInt32 a)
{
    vec_st(a.simdInternal_, 0, m);
}

static inline SimdFInt32 gmx_simdcall
setZeroFI()
{
    return {
               vec_splats(0.0)
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
               vec_neg(x.simdInternal_)
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
               vec_add(a.simdInternal_, vec_sel(vec_splats(0.0), b.simdInternal_, m.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
maskzMul(SimdFloat a, SimdFloat b, SimdFBool m)
{
    return {
               vec_sel(vec_splats(0.0), vec_mul(a.simdInternal_, b.simdInternal_), m.simdInternal_)
    };
}

static inline SimdFloat
maskzFma(SimdFloat a, SimdFloat b, SimdFloat c, SimdFBool m)
{
    return {
               vec_sel(vec_splats(0.0), vec_madd(a.simdInternal_, b.simdInternal_, c.simdInternal_), m.simdInternal_)
    };
}

static inline SimdFloat
maskzRsqrt(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = vec_sel(vec_splats(1.0), x.simdInternal_, m.simdInternal_);
#endif
    return {
               vec_sel(vec_splats(0.0), vec_rsqrte(x.simdInternal_), m.simdInternal_)
    };
}

static inline SimdFloat
maskzRcp(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = vec_sel(vec_splats(1.0), x.simdInternal_, m.simdInternal_);
#endif
    return {
               vec_sel(vec_splats(0.0), vec_re(x.simdInternal_), m.simdInternal_)
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
               vec_sel(b.simdInternal_, a.simdInternal_, vec_sub(a.simdInternal_, b.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
min(SimdFloat a, SimdFloat b)
{
    return {
               vec_sel(b.simdInternal_, a.simdInternal_, vec_sub(b.simdInternal_, a.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
round(SimdFloat x)
{
    // Note: It is critical to use vec_cfid(vec_ctid(a)) for the implementation
    // here, since vec_round() does not adhere to the FP control
    // word rounding scheme. We rely on float-to-float and float-to-integer
    // rounding being the same for half-way values in a few algorithms.
    return {
               vec_cfid(vec_ctid(x.simdInternal_))
    };
}

static inline SimdFloat gmx_simdcall
trunc(SimdFloat x)
{
    return {
               vec_trunc(x.simdInternal_)
    };
}

static inline SimdFloat
frexp(SimdFloat value, SimdFInt32 * exponent)
{
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)  rdata[GMX_SIMD_FLOAT_WIDTH];
    GMX_ALIGNED(int, GMX_SIMD_FLOAT_WIDTH)    idata[GMX_SIMD_FLOAT_WIDTH];

    vec_st(value.simdInternal_, 0, rdata);

    for (std::size_t i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        rdata[i] = std::frexp(rdata[i], idata + i);
    }

    exponent->simdInternal_ = vec_ldia(0, idata);
    value.simdInternal_     = vec_ld(0, rdata);

    return value;
}

template <MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat
ldexp(SimdFloat value, SimdFInt32 exponent)
{
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)  rdata[GMX_SIMD_FLOAT_WIDTH];
    GMX_ALIGNED(int, GMX_SIMD_FLOAT_WIDTH)    idata[GMX_SIMD_FLOAT_WIDTH];

    vec_st(value.simdInternal_,    0, rdata);
    vec_st(exponent.simdInternal_, 0, idata);

    for (std::size_t i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        rdata[i] = std::ldexp(rdata[i], idata[i]);
    }

    value.simdInternal_     = vec_ld(0, rdata);

    return value;
}

static inline float gmx_simdcall
reduce(SimdFloat x)
{
    vector4double y = vec_sldw(x.simdInternal_, x.simdInternal_, 2);
    vector4double z;

    y = vec_add(y, x.simdInternal_);
    z = vec_sldw(y, y, 1);
    y = vec_add(y, z);
    return vec_extract(y, 0);
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
               vec_not(vec_cmpeq(a.simdInternal_, b.simdInternal_))
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
               vec_or(vec_cmplt(a.simdInternal_, b.simdInternal_), vec_cmpeq(a.simdInternal_, b.simdInternal_))
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
    vector4double b = vec_sldw(a.simdInternal_, a.simdInternal_, 2);

    a.simdInternal_ = vec_or(a.simdInternal_, b);
    b               = vec_sldw(a.simdInternal_, a.simdInternal_, 1);
    b               = vec_or(a.simdInternal_, b);
    return (vec_extract(b, 0) > 0);
}

static inline SimdFloat gmx_simdcall
selectByMask(SimdFloat a, SimdFBool m)
{
    return {
               vec_sel(vec_splats(0.0), a.simdInternal_, m.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
selectByNotMask(SimdFloat a, SimdFBool m)
{
    return {
               vec_sel(a.simdInternal_, vec_splats(0.0), m.simdInternal_)
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
cvtR2I(SimdFloat a)
{
    return {
               vec_ctiw(a.simdInternal_)
    };
}

static inline SimdFInt32 gmx_simdcall
cvttR2I(SimdFloat a)
{
    return {
               vec_ctiwz(a.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
cvtI2R(SimdFInt32 a)
{
    return {
               vec_cfid(a.simdInternal_)
    };
}

static inline SimdFloat gmx_simdcall
copysign(SimdFloat x, SimdFloat y)
{
    return {
               vec_cpsgn(y.simdInternal_, x.simdInternal_)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_FLOAT_H
