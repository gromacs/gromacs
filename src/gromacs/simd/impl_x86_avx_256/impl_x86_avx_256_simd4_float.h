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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_256_SIMD4_FLOAT_H

#include "config.h"

#include <cassert>
#include <cstddef>

#include <immintrin.h>

namespace gmx
{

class Simd4Float
{
    public:
        Simd4Float() {}

        Simd4Float(float f) : simdInternal_(_mm_set1_ps(f)) {}

        // Internal utility constructor to simplify return statements
        Simd4Float(__m128 simd) : simdInternal_(simd) {}

        __m128  simdInternal_;
};

class Simd4FBool
{
    public:
        Simd4FBool() {}

        //! \brief Construct from scalar bool
        Simd4FBool(bool b) : simdInternal_(_mm_castsi128_ps(_mm_set1_epi32( b ? 0xFFFFFFFF : 0))) {}

        // Internal utility constructor to simplify return statements
        Simd4FBool(__m128 simd) : simdInternal_(simd) {}

        __m128  simdInternal_;
};

static inline Simd4Float gmx_simdcall
load4(const float *m)
{
    assert(std::size_t(m) % 16 == 0);
    return {
               _mm_load_ps(m)
    };
}

static inline void gmx_simdcall
store4(float *m, Simd4Float a)
{
    assert(std::size_t(m) % 16 == 0);
    _mm_store_ps(m, a.simdInternal_);
}

static inline Simd4Float gmx_simdcall
load4U(const float *m)
{
    return {
               _mm_loadu_ps(m)
    };
}

static inline void gmx_simdcall
store4U(float *m, Simd4Float a)
{
    _mm_storeu_ps(m, a.simdInternal_);
}

static inline Simd4Float gmx_simdcall
simd4SetZeroF()
{
    return {
               _mm_setzero_ps()
    };
}

static inline Simd4Float gmx_simdcall
operator&(Simd4Float a, Simd4Float b)
{
    return {
               _mm_and_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
andNot(Simd4Float a, Simd4Float b)
{
    return {
               _mm_andnot_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
operator|(Simd4Float a, Simd4Float b)
{
    return {
               _mm_or_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
operator^(Simd4Float a, Simd4Float b)
{
    return {
               _mm_xor_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
operator+(Simd4Float a, Simd4Float b)
{
    return {
               _mm_add_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
operator-(Simd4Float a, Simd4Float b)
{
    return {
               _mm_sub_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
operator-(Simd4Float x)
{
    return {
               _mm_xor_ps(x.simdInternal_, _mm_set1_ps(GMX_FLOAT_NEGZERO))
    };
}

static inline Simd4Float gmx_simdcall
operator*(Simd4Float a, Simd4Float b)
{
    return {
               _mm_mul_ps(a.simdInternal_, b.simdInternal_)
    };
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
static inline Simd4Float gmx_simdcall
fma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               _mm_add_ps(_mm_mul_ps(a.simdInternal_, b.simdInternal_), c.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
fms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               _mm_sub_ps(_mm_mul_ps(a.simdInternal_, b.simdInternal_), c.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
fnma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               _mm_sub_ps(c.simdInternal_, _mm_mul_ps(a.simdInternal_, b.simdInternal_))
    };
}

static inline Simd4Float gmx_simdcall
fnms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               _mm_sub_ps(_mm_setzero_ps(), _mm_add_ps(_mm_mul_ps(a.simdInternal_, b.simdInternal_), c.simdInternal_))
    };
}
#endif

static inline Simd4Float gmx_simdcall
rsqrt(Simd4Float x)
{
    return {
               _mm_rsqrt_ps(x.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
abs(Simd4Float x)
{
    return {
               _mm_andnot_ps( _mm_set1_ps(GMX_FLOAT_NEGZERO), x.simdInternal_ )
    };
}

static inline Simd4Float gmx_simdcall
max(Simd4Float a, Simd4Float b)
{
    return {
               _mm_max_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
min(Simd4Float a, Simd4Float b)
{
    return {
               _mm_min_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
round(Simd4Float x)
{
    return {
               _mm_round_ps(x.simdInternal_, _MM_FROUND_NINT)
    };
}

static inline Simd4Float gmx_simdcall
trunc(Simd4Float x)
{
    return {
               _mm_round_ps(x.simdInternal_, _MM_FROUND_TRUNC)
    };
}

static inline float gmx_simdcall
dotProduct(Simd4Float a, Simd4Float b)
{
    __m128 c, d;
    c = _mm_mul_ps(a.simdInternal_, b.simdInternal_);
    d = _mm_add_ps(c, _mm_permute_ps(c, _MM_SHUFFLE(2, 1, 2, 1)));
    d = _mm_add_ps(d, _mm_permute_ps(c, _MM_SHUFFLE(3, 2, 3, 2)));
    return *reinterpret_cast<float *>(&d);
}

static inline void gmx_simdcall
transpose(Simd4Float * v0, Simd4Float * v1,
          Simd4Float * v2, Simd4Float * v3)
{
    _MM_TRANSPOSE4_PS(v0->simdInternal_, v1->simdInternal_, v2->simdInternal_, v3->simdInternal_);
}

static inline Simd4FBool gmx_simdcall
operator==(Simd4Float a, Simd4Float b)
{
    return {
               _mm_cmp_ps(a.simdInternal_, b.simdInternal_, _CMP_EQ_OQ)
    };
}

static inline Simd4FBool gmx_simdcall
operator!=(Simd4Float a, Simd4Float b)
{
    return {
               _mm_cmp_ps(a.simdInternal_, b.simdInternal_, _CMP_NEQ_OQ)
    };
}

static inline Simd4FBool gmx_simdcall
operator<(Simd4Float a, Simd4Float b)
{
    return {
               _mm_cmp_ps(a.simdInternal_, b.simdInternal_, _CMP_LT_OQ)
    };
}

static inline Simd4FBool gmx_simdcall
operator<=(Simd4Float a, Simd4Float b)
{
    return {
               _mm_cmp_ps(a.simdInternal_, b.simdInternal_, _CMP_LE_OQ)
    };
}

static inline Simd4FBool gmx_simdcall
operator&&(Simd4FBool a, Simd4FBool b)
{
    return {
               _mm_and_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4FBool gmx_simdcall
operator||(Simd4FBool a, Simd4FBool b)
{
    return {
               _mm_or_ps(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(Simd4FBool a) { return _mm_movemask_ps(a.simdInternal_) != 0; }

static inline Simd4Float gmx_simdcall
selectByMask(Simd4Float a, Simd4FBool mask)
{
    return {
               _mm_and_ps(a.simdInternal_, mask.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
selectByNotMask(Simd4Float a, Simd4FBool mask)
{
    return {
               _mm_andnot_ps(mask.simdInternal_, a.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
blend(Simd4Float a, Simd4Float b, Simd4FBool sel)
{
    return {
               _mm_blendv_ps(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

static inline float gmx_simdcall
reduce(Simd4Float a)
{
    __m128 b;
    b = _mm_add_ps(a.simdInternal_, _mm_permute_ps(a.simdInternal_, _MM_SHUFFLE(1, 0, 3, 2)));
    b = _mm_add_ss(b, _mm_permute_ps(b, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float *>(&b);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_256_SIMD4_FLOAT_H
