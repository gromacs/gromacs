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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD4_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD4_FLOAT_H

#include "config.h"

#include "gromacs/simd/impl_ibm_vsx/impl_ibm_vsx_definitions.h"
#include "gromacs/simd/impl_ibm_vsx/impl_ibm_vsx_simd_float.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

class Simd4Float
{
public:
    Simd4Float() {}

    Simd4Float(float f) : simdInternal_(vec_splats(f)) {}

    // Internal utility constructor to simplify return statements
    Simd4Float(__vector float simd) : simdInternal_(simd) {}

    __vector float simdInternal_;
};

class Simd4FBool
{
public:
    Simd4FBool() {}

    //! \brief Construct from scalar bool
    Simd4FBool(bool b) :
        simdInternal_(reinterpret_cast<__vector vsxBool int>(vec_splats(b ? 0xFFFFFFFF : 0)))
    {
    }

    // Internal utility constructor to simplify return statements
    Simd4FBool(__vector vsxBool int simd) : simdInternal_(simd) {}

    __vector vsxBool int simdInternal_;
};

// The VSX load & store operations are a bit of a mess. The interface is different
// for xlc version 12, xlc version 13, and gcc. Long-term IBM recommends
// simply using pointer dereferencing both for aligned and unaligned loads.
// That's nice, but unfortunately xlc still bugs out when the pointer is
// not aligned. Sticking to vec_xl/vec_xst isn't a solution either, since
// that appears to be buggy for some _aligned_ loads :-)
//
// For now, we use pointer dereferencing for all aligned load/stores, and
// for unaligned ones with gcc. On xlc we use vec_xlw4/vec_xstw4 for
// unaligned memory operations. The latest docs recommend using the overloaded
// vec_xl/vec_xst, but that is not supported on xlc version 12. We'll
// revisit things once xlc is a bit more stable - for now you probably want
// to stick to gcc...

static inline Simd4Float gmx_simdcall load4(const float* m)
{
    return { *reinterpret_cast<const __vector float*>(m) };
}

static inline void gmx_simdcall store4(float* m, Simd4Float a)
{
    *reinterpret_cast<__vector float*>(m) = a.simdInternal_;
}

static inline Simd4Float gmx_simdcall load4U(const float* m)
{
    return { vec_xl(0, m) };
}

static inline void gmx_simdcall store4U(float* m, Simd4Float a)
{
    vec_xst(a.simdInternal_, 0, m);
}

static inline Simd4Float gmx_simdcall simd4SetZeroF()
{
    return { vec_splats(0.0F) };
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
    return { -x.simdInternal_ };
}

static inline Simd4Float gmx_simdcall operator*(Simd4Float a, Simd4Float b)
{
    return { vec_mul(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { vec_madd(a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { vec_msub(a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fnma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { vec_nmsub(a.simdInternal_, b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fnms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { vec_nmadd(a.simdInternal_, b.simdInternal_, c.simdInternal_) };
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
    const __vector unsigned char perm1 = { 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7 };
    const __vector unsigned char perm2 = { 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3 };
    __vector float               c     = vec_mul(a.simdInternal_, b.simdInternal_);
    __vector float               sum;
    sum = vec_add(c, vec_perm(c, c, perm1));
    sum = vec_add(sum, vec_perm(c, c, perm2));
    return vec_extract(sum, 0);
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
    return vec_any_ne(a.simdInternal_, reinterpret_cast<__vector vsxBool int>(vec_splats(0)));
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

static inline float gmx_simdcall reduce(Simd4Float x)
{
    const __vector unsigned char perm1 = { 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7 };
    const __vector unsigned char perm2 = { 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3 };

    x.simdInternal_ = vec_add(x.simdInternal_, vec_perm(x.simdInternal_, x.simdInternal_, perm1));
    x.simdInternal_ = vec_add(x.simdInternal_, vec_perm(x.simdInternal_, x.simdInternal_, perm2));
    return vec_extract(x.simdInternal_, 0);
}

} // namespace gmx

#endif // GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD4_FLOAT_H
