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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VMX_SIMD4_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VMX_SIMD4_FLOAT_H

#include "config.h"

#include "gromacs/utility/basedefinitions.h"

#include "impl_ibm_vmx_definitions.h"

namespace gmx
{

class Simd4Float
{
public:
    Simd4Float() {}

    Simd4Float(float f)
    {
        __vector unsigned char perm;

        simdInternal_ = vec_lde(0, const_cast<float*>(&f));
        perm          = vec_lvsl(0, const_cast<float*>(&f));
        simdInternal_ = vec_perm(simdInternal_, simdInternal_, perm);
        simdInternal_ = vec_splat(simdInternal_, 0);
    }

    // Internal utility constructor to simplify return statements
    Simd4Float(__vector float simd) : simdInternal_(simd) {}

    __vector float simdInternal_;
};

class Simd4FBool
{
public:
    Simd4FBool() {}

    Simd4FBool(bool b)
    {
        simdInternal_ = reinterpret_cast<__vector vmxBool int>(vec_splat_u32(b ? 0xFFFFFFFF : 0));
    }

    // Internal utility constructor to simplify return statements
    Simd4FBool(__vector vmxBool int simd) : simdInternal_(simd) {}

    __vector vmxBool int simdInternal_;
};

static inline Simd4Float gmx_simdcall load4(const float* m)
{
    return { vec_ld(0, const_cast<float*>(m)) };
}

static inline void gmx_simdcall store4(float* m, Simd4Float a)
{
    vec_st(a.simdInternal_, 0, const_cast<float*>(m));
}

static inline Simd4Float gmx_simdcall simd4SetZeroF()
{
    return { reinterpret_cast<__vector float>(vec_splat_u32(0)) };
}

static inline Simd4Float gmx_simdcall operator&(Simd4Float a, Simd4Float b)
{
    return { vec_and(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall andNot(Simd4Float a, Simd4Float b)
{
    return { vec_andc(b.simdInternal_, a.simdInternal_) };
}

static inline Simd4Float gmx_simdcall operator|(Simd4Float a, Simd4Float b)
{
    return { vec_or(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall operator^(Simd4Float a, Simd4Float b)
{
    return { vec_xor(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall operator+(Simd4Float a, Simd4Float b)
{
    return { vec_add(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall operator-(Simd4Float a, Simd4Float b)
{
    return { vec_sub(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall operator-(Simd4Float x)
{
    return { vec_xor(x.simdInternal_,
                     reinterpret_cast<__vector float>(vec_sl(vec_splat_u32(-1), vec_splat_u32(-1)))) };
}

static inline Simd4Float gmx_simdcall operator*(Simd4Float a, Simd4Float b)
{
    return { vec_madd(a.simdInternal_, b.simdInternal_,
                      reinterpret_cast<__vector float>(vec_splat_u32(0))) };
}

static inline Simd4Float gmx_simdcall fma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { vec_madd(a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { vec_madd(a.simdInternal_, b.simdInternal_, -c.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fnma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { vec_nmsub(a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fnms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { -vec_madd(a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Float gmx_simdcall rsqrt(Simd4Float x)
{
    return { vec_rsqrte(x.simdInternal_) };
}

static inline Simd4Float gmx_simdcall abs(Simd4Float x)
{
    return { vec_abs(x.simdInternal_) };
}

static inline Simd4Float gmx_simdcall max(Simd4Float a, Simd4Float b)
{
    return { vec_max(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall min(Simd4Float a, Simd4Float b)
{
    return { vec_min(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall round(Simd4Float x)
{
    return { vec_round(x.simdInternal_) };
}

static inline Simd4Float gmx_simdcall trunc(Simd4Float x)
{
    return { vec_trunc(x.simdInternal_) };
}

static inline float gmx_simdcall dotProduct(Simd4Float a, Simd4Float b)
{
    float res;

    __vector float c = vec_madd(a.simdInternal_, b.simdInternal_,
                                reinterpret_cast<__vector float>(vec_splat_u32(0)));
    // Keep only elements 0,1,2 by shifting in zero from right (xor of a vector with itself is 0)
    c = vec_sld(c, vec_xor(a.simdInternal_, a.simdInternal_), 4);
    // calculate sum
    c = vec_add(c, vec_sld(c, c, 8));
    c = vec_add(c, vec_sld(c, c, 4));
    vec_ste(c, 0, &res);
    return res;
}

static inline void gmx_simdcall transpose(Simd4Float* v0, Simd4Float* v1, Simd4Float* v2, Simd4Float* v3)
{
    __vector float t0 = vec_mergeh(v0->simdInternal_, v2->simdInternal_);
    __vector float t1 = vec_mergel(v0->simdInternal_, v2->simdInternal_);
    __vector float t2 = vec_mergeh(v1->simdInternal_, v3->simdInternal_);
    __vector float t3 = vec_mergel(v1->simdInternal_, v3->simdInternal_);
    v0->simdInternal_ = vec_mergeh(t0, t2);
    v1->simdInternal_ = vec_mergel(t0, t2);
    v2->simdInternal_ = vec_mergeh(t1, t3);
    v3->simdInternal_ = vec_mergel(t1, t3);
}

static inline Simd4FBool gmx_simdcall operator==(Simd4Float a, Simd4Float b)
{
    return { vec_cmpeq(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4FBool gmx_simdcall operator!=(Simd4Float a, Simd4Float b)
{
    return { vec_or(vec_cmpgt(a.simdInternal_, b.simdInternal_),
                    vec_cmplt(a.simdInternal_, b.simdInternal_)) };
}

static inline Simd4FBool gmx_simdcall operator<(Simd4Float a, Simd4Float b)
{
    return { vec_cmplt(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4FBool gmx_simdcall operator<=(Simd4Float a, Simd4Float b)
{
    return { vec_cmple(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4FBool gmx_simdcall operator&&(Simd4FBool a, Simd4FBool b)
{
    return { vec_and(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4FBool gmx_simdcall operator||(Simd4FBool a, Simd4FBool b)
{
    return { vec_or(a.simdInternal_, b.simdInternal_) };
}

static inline bool gmx_simdcall anyTrue(Simd4FBool a)
{
    return vec_any_ne(a.simdInternal_, reinterpret_cast<__vector vmxBool int>(vec_splat_u32(0)));
}

static inline Simd4Float gmx_simdcall selectByMask(Simd4Float a, Simd4FBool m)
{
    return { vec_and(a.simdInternal_, reinterpret_cast<__vector float>(m.simdInternal_)) };
}

static inline Simd4Float gmx_simdcall selectByNotMask(Simd4Float a, Simd4FBool m)
{
    return { vec_andc(a.simdInternal_, reinterpret_cast<__vector float>(m.simdInternal_)) };
}

static inline Simd4Float gmx_simdcall blend(Simd4Float a, Simd4Float b, Simd4FBool sel)
{
    return { vec_sel(a.simdInternal_, b.simdInternal_, sel.simdInternal_) };
}

static inline float gmx_simdcall reduce(Simd4Float a)
{
    __vector float c = a.simdInternal_;
    float          res;

    // calculate sum
    c = vec_add(c, vec_sld(c, c, 8));
    c = vec_add(c, vec_sld(c, c, 4));
    vec_ste(c, 0, &res);
    return res;
}

} // namespace gmx

#endif // GMX_SIMD_IMPLEMENTATION_IBM_VMX_SIMD4_FLOAT_H
