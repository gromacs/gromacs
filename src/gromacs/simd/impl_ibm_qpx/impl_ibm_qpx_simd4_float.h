/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD4_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD4_FLOAT_H

#include "config.h"

// Assert is buggy on xlc with high optimization, so we skip it for QPX
#include <cstddef>

#ifdef __clang__
#include <qpxmath.h>
#endif

namespace gmx
{

class Simd4Float
{
    public:
        Simd4Float() {}

        Simd4Float(float d) : simdInternal_(vec_splats(d)) {}

        // Internal utility constructor to simplify return statements
        Simd4Float(vector4double simd) : simdInternal_(simd) {}

        vector4double  simdInternal_;
};

class Simd4FBool
{
    public:
        Simd4FBool() {}

        //! \brief Construct from scalar bool
        Simd4FBool(bool b) : simdInternal_(vec_splats(b ? 1.0 : -1.0)) {}

        // Internal utility constructor to simplify return statements
        Simd4FBool(vector4double simd) : simdInternal_(simd) {}

        vector4double  simdInternal_;
};

static inline Simd4Float gmx_simdcall
load4(const float *m)
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
store4(float *m, Simd4Float a)
{
#ifdef NDEBUG
    vec_st(a.simdInternal_, 0, m);
#else
    vec_sta(a.simdInternal_, 0, m);
#endif
}

static inline Simd4Float gmx_simdcall
simd4SetZeroF()
{
    return {
               vec_splats(0.0)
    };
}

static inline Simd4Float gmx_simdcall
operator+(Simd4Float a, Simd4Float b)
{
    return {
               vec_add(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
operator-(Simd4Float a, Simd4Float b)
{
    return {
               vec_sub(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
operator-(Simd4Float x)
{
    return {
               vec_neg(x.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
operator*(Simd4Float a, Simd4Float b)
{
    return {
               vec_mul(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
fma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               vec_madd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
fms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               vec_msub(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
fnma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               vec_nmsub(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
fnms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               vec_nmadd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
rsqrt(Simd4Float x)
{
    return {
               vec_rsqrte(x.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
abs(Simd4Float x)
{
    return {
               vec_abs( x.simdInternal_ )
    };
}

static inline Simd4Float gmx_simdcall
max(Simd4Float a, Simd4Float b)
{
    return {
               vec_sel(b.simdInternal_, a.simdInternal_, vec_sub(a.simdInternal_, b.simdInternal_))
    };
}

static inline Simd4Float gmx_simdcall
min(Simd4Float a, Simd4Float b)
{
    return {
               vec_sel(b.simdInternal_, a.simdInternal_, vec_sub(b.simdInternal_, a.simdInternal_))
    };
}

static inline Simd4Float gmx_simdcall
round(Simd4Float x)
{
    // Note: It is critical to use vec_cfid(vec_ctid(a)) for the implementation
    // here, since vec_round() does not adhere to the FP control
    // word rounding scheme. We rely on float-to-float and float-to-integer
    // rounding being the same for half-way values in a few algorithms.

    return {
               vec_cfid(vec_ctid(x.simdInternal_))
    };
}

static inline Simd4Float gmx_simdcall
trunc(Simd4Float x)
{
    return {
               vec_trunc(x.simdInternal_)
    };
}

static inline float gmx_simdcall
dotProduct(Simd4Float a, Simd4Float b)
{
    vector4double dp_sh0 = vec_mul(a.simdInternal_, b.simdInternal_);
    vector4double dp_sh1 = vec_sldw(dp_sh0, dp_sh0, 1);
    vector4double dp_sh2 = vec_sldw(dp_sh0, dp_sh0, 2);
    vector4double dp     = vec_add(dp_sh2, vec_add(dp_sh0, dp_sh1));

    return vec_extract(dp, 0);
}

static inline void gmx_simdcall
transpose(Simd4Float * v0, Simd4Float * v1,
          Simd4Float * v2, Simd4Float * v3)
{
    vector4double t0 = vec_perm(v0->simdInternal_, v2->simdInternal_, vec_gpci(00415));
    vector4double t1 = vec_perm(v0->simdInternal_, v2->simdInternal_, vec_gpci(02637));
    vector4double t2 = vec_perm(v1->simdInternal_, v3->simdInternal_, vec_gpci(00415));
    vector4double t3 = vec_perm(v1->simdInternal_, v3->simdInternal_, vec_gpci(02637));
    v0->simdInternal_ = vec_perm(t0, t2, vec_gpci(00415));
    v1->simdInternal_ = vec_perm(t0, t2, vec_gpci(02637));
    v2->simdInternal_ = vec_perm(t1, t3, vec_gpci(00415));
    v3->simdInternal_ = vec_perm(t1, t3, vec_gpci(02637));
}

static inline Simd4FBool gmx_simdcall
operator==(Simd4Float a, Simd4Float b)
{
    return {
               vec_cmpeq(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4FBool gmx_simdcall
operator!=(Simd4Float a, Simd4Float b)
{
    return {
               vec_not(vec_cmpeq(a.simdInternal_, b.simdInternal_))
    };
}

static inline Simd4FBool gmx_simdcall
operator<(Simd4Float a, Simd4Float b)
{
    return {
               vec_cmplt(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4FBool gmx_simdcall
operator<=(Simd4Float a, Simd4Float b)
{
    return {
               vec_or(vec_cmplt(a.simdInternal_, b.simdInternal_), vec_cmpeq(a.simdInternal_, b.simdInternal_))
    };
}

static inline Simd4FBool gmx_simdcall
operator&&(Simd4FBool a, Simd4FBool b)
{
    return {
               vec_and(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4FBool gmx_simdcall
operator||(Simd4FBool a, Simd4FBool b)
{
    return {
               vec_or(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(Simd4FBool a)
{
    vector4double b = vec_sldw(a.simdInternal_, a.simdInternal_, 2);

    a.simdInternal_ = vec_or(a.simdInternal_, b);
    b               = vec_sldw(a.simdInternal_, a.simdInternal_, 1);
    b               = vec_or(a.simdInternal_, b);
    return (vec_extract(b, 0) > 0);
}

static inline Simd4Float gmx_simdcall
selectByMask(Simd4Float a, Simd4FBool m)
{
    return {
               vec_sel(vec_splats(0.0), a.simdInternal_, m.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
selectByNotMask(Simd4Float a, Simd4FBool m)
{
    return {
               vec_sel(a.simdInternal_, vec_splats(0.0), m.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
blend(Simd4Float a, Simd4Float b, Simd4FBool sel)
{
    return {
               vec_sel(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

static inline float gmx_simdcall
reduce(Simd4Float x)
{
    vector4double y = vec_sldw(x.simdInternal_, x.simdInternal_, 2);
    vector4double z;

    y = vec_add(y, x.simdInternal_);
    z = vec_sldw(y, y, 1);
    y = vec_add(y, z);
    return vec_extract(y, 0);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD4_FLOAT_H
