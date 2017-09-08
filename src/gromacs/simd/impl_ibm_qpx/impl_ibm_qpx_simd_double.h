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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_DOUBLE_H
#define GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_DOUBLE_H

#include "config.h"

// Assert is buggy on xlc with high optimization, so we skip it for QPX
#include <cmath>
#include <cstddef>
#include <cstdint>

#ifdef __clang__
#include <qpxmath.h>
#endif

#include "gromacs/math/utilities.h"
#include "gromacs/utility/basedefinitions.h"

#include "impl_ibm_qpx_simd_float.h"

namespace gmx
{

class SimdDouble
{
    public:
        SimdDouble() {}

        SimdDouble(double d) : simdInternal_(vec_splats(d)) {}

        // Internal utility constructor to simplify return statements
        SimdDouble(vector4double simd) : simdInternal_(simd) {}

        vector4double  simdInternal_;
};

class SimdDInt32
{
    public:
        SimdDInt32() {}

        SimdDInt32(std::int32_t i)
        {
            GMX_ALIGNED(int, GMX_SIMD_DINT32_WIDTH) idata[GMX_SIMD_DINT32_WIDTH];
            idata[0]      = i;
            simdInternal_ = vec_splat(vec_ldia(0, idata), 0);
        }

        // Internal utility constructor to simplify return statements
        SimdDInt32(vector4double simd) : simdInternal_(simd) {}

        vector4double  simdInternal_;
};

class SimdDBool
{
    public:
        SimdDBool() {}

        SimdDBool(bool b) : simdInternal_(vec_splats(b ? 1.0 : -1.0)) {}

        // Internal utility constructor to simplify return statements
        SimdDBool(vector4double simd) : simdInternal_(simd) {}

        vector4double  simdInternal_;
};

static inline SimdDouble gmx_simdcall
simdLoad(const double *m)
{
#ifdef NDEBUG
    return {
               vec_ld(0, const_cast<double *>(m))
    };
#else
    return {
               vec_lda(0, const_cast<double *>(m))
    };
#endif
}

static inline void gmx_simdcall
store(double *m, SimdDouble a)
{
#ifdef NDEBUG
    vec_st(a.simdInternal_, 0, m);
#else
    vec_sta(a.simdInternal_, 0, m);
#endif
}

static inline SimdDouble gmx_simdcall
setZeroD()
{
    return {
               vec_splats(0.0)
    };
}

static inline SimdDInt32 gmx_simdcall
simdLoadDI(const std::int32_t * m)
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
store(std::int32_t * m, SimdDInt32 a)
{
    vec_st(a.simdInternal_, 0, m);
}

static inline SimdDInt32 gmx_simdcall
setZeroDI()
{
    return {
               vec_splats(0.0)
    };
}

static inline SimdDouble gmx_simdcall
operator+(SimdDouble a, SimdDouble b)
{
    return {
               vec_add(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator-(SimdDouble a, SimdDouble b)
{
    return {
               vec_sub(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator-(SimdDouble x)
{
    return {
               vec_neg(x.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
operator*(SimdDouble a, SimdDouble b)
{
    return {
               vec_mul(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               vec_madd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               vec_msub(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fnma(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               vec_nmsub(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
fnms(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               vec_nmadd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
rsqrt(SimdDouble x)
{
    return {
               vec_rsqrte(x.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
rcp(SimdDouble x)
{
    return {
               vec_re(x.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
maskAdd(SimdDouble a, SimdDouble b, SimdDBool m)
{
    return {
               vec_add(a.simdInternal_, vec_sel(vec_splats(0.0), b.simdInternal_, m.simdInternal_))
    };
}

static inline SimdDouble gmx_simdcall
maskzMul(SimdDouble a, SimdDouble b, SimdDBool m)
{
    return {
               vec_sel(vec_splats(0.0), vec_mul(a.simdInternal_, b.simdInternal_), m.simdInternal_)
    };
}

static inline SimdDouble
maskzFma(SimdDouble a, SimdDouble b, SimdDouble c, SimdDBool m)
{
    return {
               vec_sel(vec_splats(0.0), vec_madd(a.simdInternal_, b.simdInternal_, c.simdInternal_), m.simdInternal_)
    };
}

static inline SimdDouble
maskzRsqrt(SimdDouble x, SimdDBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = vec_sel(vec_splats(1.0), x.simdInternal_, m.simdInternal_);
#endif
    return {
               vec_sel(vec_splats(0.0), vec_rsqrte(x.simdInternal_), m.simdInternal_)
    };
}

static inline SimdDouble
maskzRcp(SimdDouble x, SimdDBool m)
{
#ifndef NDEBUG
    x.simdInternal_ = vec_sel(vec_splats(1.0), x.simdInternal_, m.simdInternal_);
#endif
    return {
               vec_sel(vec_splats(0.0), vec_re(x.simdInternal_), m.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
abs(SimdDouble x)
{
    return {
               vec_abs( x.simdInternal_ )
    };
}

static inline SimdDouble gmx_simdcall
max(SimdDouble a, SimdDouble b)
{
    return {
               vec_sel(b.simdInternal_, a.simdInternal_, vec_sub(a.simdInternal_, b.simdInternal_))
    };
}

static inline SimdDouble gmx_simdcall
min(SimdDouble a, SimdDouble b)
{
    return {
               vec_sel(b.simdInternal_, a.simdInternal_, vec_sub(b.simdInternal_, a.simdInternal_))
    };
}

static inline SimdDouble gmx_simdcall
round(SimdDouble x)
{
    // Note: It is critical to use vec_cfid(vec_ctid(a)) for the implementation
    // here, since vec_round() does not adhere to the FP control
    // word rounding scheme. We rely on float-to-float and float-to-integer
    // rounding being the same for half-way values in a few algorithms.
    return {
               vec_cfid(vec_ctid(x.simdInternal_))
    };
}

static inline SimdDouble gmx_simdcall
trunc(SimdDouble x)
{
    return {
               vec_trunc(x.simdInternal_)
    };
}

static inline SimdDouble
frexp(SimdDouble value, SimdDInt32 * exponent)
{
    GMX_ALIGNED(double, GMX_SIMD_DOUBLE_WIDTH) rdata[GMX_SIMD_DOUBLE_WIDTH];
    GMX_ALIGNED(int, GMX_SIMD_DOUBLE_WIDTH)    idata[GMX_SIMD_DOUBLE_WIDTH];

    vec_st(value.simdInternal_, 0, rdata);

    for (std::size_t i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        rdata[i] = std::frexp(rdata[i], idata + i);
    }

    exponent->simdInternal_ = vec_ldia(0, idata);
    value.simdInternal_     = vec_ld(0, rdata);

    return value;
}

template <MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble
ldexp(SimdDouble value, SimdDInt32 exponent)
{
    GMX_ALIGNED(double, GMX_SIMD_DOUBLE_WIDTH) rdata[GMX_SIMD_DOUBLE_WIDTH];
    GMX_ALIGNED(int, GMX_SIMD_DOUBLE_WIDTH)    idata[GMX_SIMD_DOUBLE_WIDTH];

    vec_st(value.simdInternal_,    0, rdata);
    vec_st(exponent.simdInternal_, 0, idata);

    for (std::size_t i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        rdata[i] = std::ldexp(rdata[i], idata[i]);
    }

    value.simdInternal_     = vec_ld(0, rdata);

    return value;
}

static inline double gmx_simdcall
reduce(SimdDouble x)
{
    vector4double y = vec_sldw(x.simdInternal_, x.simdInternal_, 2);
    vector4double z;

    y = vec_add(y, x.simdInternal_);
    z = vec_sldw(y, y, 1);
    y = vec_add(y, z);
    return vec_extract(y, 0);
}

static inline SimdDBool gmx_simdcall
operator==(SimdDouble a, SimdDouble b)
{
    return {
               vec_cmpeq(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDBool gmx_simdcall
operator!=(SimdDouble a, SimdDouble b)
{
    return {
               vec_not(vec_cmpeq(a.simdInternal_, b.simdInternal_))
    };
}

static inline SimdDBool gmx_simdcall
operator<(SimdDouble a, SimdDouble b)
{
    return {
               vec_cmplt(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDBool gmx_simdcall
operator<=(SimdDouble a, SimdDouble b)
{
    return {
               vec_or(vec_cmplt(a.simdInternal_, b.simdInternal_), vec_cmpeq(a.simdInternal_, b.simdInternal_))
    };
}

static inline SimdDBool gmx_simdcall
operator&&(SimdDBool a, SimdDBool b)
{
    return {
               vec_and(a.simdInternal_, b.simdInternal_)
    };
}

static inline SimdDBool gmx_simdcall
operator||(SimdDBool a, SimdDBool b)
{
    return {
               vec_or(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(SimdDBool a)
{
    vector4double b = vec_sldw(a.simdInternal_, a.simdInternal_, 2);

    a.simdInternal_ = vec_or(a.simdInternal_, b);
    b               = vec_sldw(a.simdInternal_, a.simdInternal_, 1);
    b               = vec_or(a.simdInternal_, b);
    return (vec_extract(b, 0) > 0);
}

static inline SimdDouble gmx_simdcall
selectByMask(SimdDouble a, SimdDBool m)
{
    return {
               vec_sel(vec_splats(0.0), a.simdInternal_, m.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
selectByNotMask(SimdDouble a, SimdDBool m)
{
    return {
               vec_sel(a.simdInternal_, vec_splats(0.0), m.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
blend(SimdDouble a, SimdDouble b, SimdDBool sel)
{
    return {
               vec_sel(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
cvtR2I(SimdDouble a)
{
    return {
               vec_ctiw(a.simdInternal_)
    };
}

static inline SimdDInt32 gmx_simdcall
cvttR2I(SimdDouble a)
{
    return {
               vec_ctiwz(a.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
cvtI2R(SimdDInt32 a)
{
    return {
               vec_cfid(a.simdInternal_)
    };
}

static inline SimdDouble gmx_simdcall
cvtF2D(SimdFloat f)
{
    return {
               f.simdInternal_
    };
}

static inline SimdFloat gmx_simdcall
cvtD2F(SimdDouble d)
{
    return {
               d.simdInternal_
    };
}

static inline SimdDouble gmx_simdcall
copysign(SimdDouble x, SimdDouble y)
{
    return {
               vec_cpsgn(y.simdInternal_, x.simdInternal_)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_DOUBLE_H
