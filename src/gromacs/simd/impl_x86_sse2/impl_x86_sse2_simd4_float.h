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
#ifndef GMX_SIMD_IMPL_X86_SSE2_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_X86_SSE2_SIMD4_FLOAT_H

#include "config.h"

#include <cassert>
#include <cstddef>

#include <emmintrin.h>

namespace gmx
{

struct Simd4Float
{
    __m128 r;
};

struct Simd4FBool
{
    __m128 b;
};

static inline Simd4Float gmx_simdcall
simd4LoadF(const float *m)
{
    assert(size_t(m) % 16 == 0);
    return {
               _mm_load_ps(m)
    };
}

static inline Simd4Float gmx_simdcall
simd4Load1F(const float *m)
{
    return {
               _mm_load1_ps(m)
    };
}

static inline Simd4Float gmx_simdcall
simd4Set1F(float r)
{
    return {
               _mm_set1_ps(r)
    };
}

static inline Simd4Float gmx_simdcall
simd4SetZeroF()
{
    return {
               _mm_setzero_ps()
    };
}

static inline void gmx_simdcall
simd4StoreF(float *m, Simd4Float a)
{
    assert(size_t(m) % 16 == 0);
    _mm_store_ps(m, a.r);
}

static inline Simd4Float gmx_simdcall
simd4LoadUF(const float *m)
{
    return {
               _mm_loadu_ps(m)
    };
}

static inline void gmx_simdcall
simd4StoreUF(float *m, Simd4Float a)
{
    _mm_storeu_ps(m, a.r);
}

static inline Simd4Float gmx_simdcall
simd4AndF(Simd4Float a, Simd4Float b)
{
    return {
               _mm_and_ps(a.r, b.r)
    };
}

static inline Simd4Float gmx_simdcall
simd4AndNotF(Simd4Float a, Simd4Float b)
{
    return {
               _mm_andnot_ps(a.r, b.r)
    };
}

static inline Simd4Float gmx_simdcall
simd4OrF(Simd4Float a, Simd4Float b)
{
    return {
               _mm_or_ps(a.r, b.r)
    };
}

static inline Simd4Float gmx_simdcall
simd4XorF(Simd4Float a, Simd4Float b)
{
    return {
               _mm_xor_ps(a.r, b.r)
    };
}

static inline Simd4Float gmx_simdcall
simd4AddF(Simd4Float a, Simd4Float b)
{
    return {
               _mm_add_ps(a.r, b.r)
    };
}

static inline Simd4Float gmx_simdcall
simd4SubF(Simd4Float a, Simd4Float b)
{
    return {
               _mm_sub_ps(a.r, b.r)
    };
}

static inline Simd4Float gmx_simdcall
simd4MulF(Simd4Float a, Simd4Float b)
{
    return {
               _mm_mul_ps(a.r, b.r)
    };
}

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline Simd4Float gmx_simdcall
simd4FmaddF(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               _mm_add_ps(_mm_mul_ps(a.r, b.r), c.r)
    };
}

static inline Simd4Float gmx_simdcall
simd4FmsubF(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               _mm_sub_ps(_mm_mul_ps(a.r, b.r), c.r)
    };
}

static inline Simd4Float gmx_simdcall
simd4FnmaddF(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               _mm_sub_ps(c.r, _mm_mul_ps(a.r, b.r))
    };
}

static inline Simd4Float gmx_simdcall
simd4FnmsubF(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               _mm_sub_ps(_mm_setzero_ps(), _mm_add_ps(_mm_mul_ps(a.r, b.r), c.r))
    };
}
#endif

static inline Simd4Float gmx_simdcall
simd4RsqrtF(Simd4Float x)
{
    return {
               _mm_rsqrt_ps(x.r)
    };
}

static inline Simd4Float gmx_simdcall
simd4AbsF(Simd4Float x)
{
    return {
               _mm_andnot_ps( _mm_set1_ps(GMX_FLOAT_NEGZERO), x.r )
    };
}

static inline Simd4Float gmx_simdcall
simd4NegF(Simd4Float x)
{
    return {
               _mm_xor_ps(x.r, _mm_set1_ps(GMX_FLOAT_NEGZERO))
    };
}

static inline Simd4Float gmx_simdcall
simd4MaxF(Simd4Float a, Simd4Float b)
{
    return {
               _mm_max_ps(a.r, b.r)
    };
}

static inline Simd4Float gmx_simdcall
simd4MinF(Simd4Float a, Simd4Float b)
{
    return {
               _mm_min_ps(a.r, b.r)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline Simd4Float gmx_simdcall
simd4RoundF(Simd4Float x)
{
    return {
               _mm_cvtepi32_ps( _mm_cvtps_epi32(x.r) )
    };
}

static inline Simd4Float gmx_simdcall
simd4TruncF(Simd4Float x)
{
    return {
               _mm_cvtepi32_ps( _mm_cvttps_epi32(x.r) )
    };
}

static inline float gmx_simdcall
simd4DotProductF(Simd4Float a, Simd4Float b)
{
    __m128 c, d;
    c = _mm_mul_ps(a.r, b.r);
    d = _mm_add_ps(c, _mm_shuffle_ps(c, c, _MM_SHUFFLE(2, 1, 2, 1)));
    d = _mm_add_ps(d, _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 2, 3, 2)));
    return *reinterpret_cast<float *>(&d);
}
#endif

static inline void gmx_simdcall
simd4Transpose(Simd4Float * v0, Simd4Float * v1,
               Simd4Float * v2, Simd4Float * v3)
{
    _MM_TRANSPOSE4_PS(v0->r, v1->r, v2->r, v3->r);
}

static inline Simd4FBool gmx_simdcall
simd4CmpEqF(Simd4Float a, Simd4Float b)
{
    return {
               _mm_cmpeq_ps(a.r, b.r)
    };
}

static inline Simd4FBool gmx_simdcall
simd4CmpLtF(Simd4Float a, Simd4Float b)
{
    return {
               _mm_cmplt_ps(a.r, b.r)
    };
}

static inline Simd4FBool gmx_simdcall
simd4CmpLeF(Simd4Float a, Simd4Float b)
{
    return {
               _mm_cmple_ps(a.r, b.r)
    };
}

static inline Simd4FBool gmx_simdcall
simd4AndFB(Simd4FBool a, Simd4FBool b)
{
    return {
               _mm_and_ps(a.b, b.b)
    };
}

static inline Simd4FBool gmx_simdcall
simd4OrFB(Simd4FBool a, Simd4FBool b)
{
    return {
               _mm_or_ps(a.b, b.b)
    };
}

static inline bool gmx_simdcall
simd4AnyTrueFB(Simd4FBool a) { return _mm_movemask_ps(a.b) != 0; }

static inline Simd4Float gmx_simdcall
simd4MaskF(Simd4Float a, Simd4FBool mask)
{
    return {
               _mm_and_ps(a.r, mask.b)
    };
}

static inline Simd4Float gmx_simdcall
simd4MaskNotF(Simd4Float a, Simd4FBool mask)
{
    return {
               _mm_andnot_ps(mask.b, a.r)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline Simd4Float gmx_simdcall
simd4BlendF(Simd4Float a, Simd4Float b, Simd4FBool sel)
{
    return {
               _mm_or_ps(_mm_andnot_ps(sel.b, a.r), _mm_and_ps(sel.b, b.r))
    };
}
#endif

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline float gmx_simdcall
simd4ReduceF(Simd4Float a)
{
    __m128 b;
    b = _mm_add_ps(a.r, _mm_shuffle_ps(a.r, a.r, _MM_SHUFFLE(1, 0, 3, 2)));
    b = _mm_add_ss(b, _mm_shuffle_ps(b, b, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float *>(&b);
}
#endif

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE2_SIMD4_FLOAT_H
