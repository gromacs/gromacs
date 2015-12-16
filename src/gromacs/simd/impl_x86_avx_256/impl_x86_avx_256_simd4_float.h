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

#include <immintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_avx_256_common.h"
#include "impl_x86_avx_256_simd_float.h"

/****************************************************
 *      SINGLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
#define Simd4Float          __m128
#define simd4LoadF           _mm_load_ps
#define simd4Load1F          _mm_broadcast_ss
#define simd4Set1F           _mm_set1_ps
#define simd4StoreF          _mm_store_ps
#define simd4LoadUF          _mm_loadu_ps
#define simd4StoreUF         _mm_storeu_ps
#define simd4SetZeroF        _mm_setzero_ps
#define simd4AddF            _mm_add_ps
#define simd4SubF            _mm_sub_ps
#define simd4MulF            _mm_mul_ps
#define simd4FmaddF(a, b, c)   _mm_add_ps(_mm_mul_ps(a, b), c)
#define simd4FmsubF(a, b, c)   _mm_sub_ps(_mm_mul_ps(a, b), c)
#define simd4FnmaddF(a, b, c)  _mm_sub_ps(c, _mm_mul_ps(a, b))
#define simd4FnmsubF(a, b, c)  _mm_sub_ps(_mm_setzero_ps(), simd4FmaddF(a, b, c))
#define simd4AndF            _mm_and_ps
#define simd4AndNotF         _mm_andnot_ps
#define simd4OrF             _mm_or_ps
#define simd4XorF            _mm_xor_ps
#define simd4RsqrtF          _mm_rsqrt_ps
#define simd4AbsF(x)        _mm_andnot_ps(_mm_set1_ps(GMX_FLOAT_NEGZERO), x)
#define simd4NegF(x)        _mm_xor_ps(x, _mm_set1_ps(GMX_FLOAT_NEGZERO))
#define simd4MaxF            _mm_max_ps
#define simd4MinF            _mm_min_ps
#define simd4RoundF(x)       _mm_round_ps(x, _MM_FROUND_NINT)
#define simd4TruncF(x)       _mm_round_ps(x, _MM_FROUND_TRUNC)
#define simd4DotProductF    simd4DotProductF_avx_256
#define Simd4FBool           __m128
#define simd4CmpEqF          _mm_cmpeq_ps
#define simd4CmpLtF          _mm_cmplt_ps
#define simd4CmpLeF          _mm_cmple_ps
#define simd4AndFB           _mm_and_ps
#define simd4OrFB            _mm_or_ps
#define simd4AnyTrueFB       _mm_movemask_ps
#define simd4MaskF      _mm_and_ps
#define simd4MaskNotF(a, sel)  _mm_andnot_ps(sel, a)
#define simd4BlendF         _mm_blendv_ps
#define simd4ReduceF         simd4ReduceF_avx_256


/* SIMD4 reduce helper */
static inline float gmx_simdcall
simd4ReduceF_avx_256(__m128 a)
{
    float f;
    a = _mm_hadd_ps(a, a);
    a = _mm_hadd_ps(a, a);
    _mm_store_ss(&f, a);
    return f;
}

/* SIMD4 Dot product helper function */
static inline float gmx_simdcall
simd4DotProductF_avx_256(__m128 a, __m128 b)
{
    float  f;
    __m128 c;
    a = _mm_mul_ps(a, b);
    c = _mm_add_ps(a, _mm_permute_ps(a, _MM_SHUFFLE(0, 3, 2, 1)));
    c = _mm_add_ps(c, _mm_permute_ps(a, _MM_SHUFFLE(1, 0, 3, 2)));
    _mm_store_ss(&f, c);
    return f;
}

#endif /* GMX_SIMD_IMPL_X86_AVX_256_SIMD4_FLOAT_H */
