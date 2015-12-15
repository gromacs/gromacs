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

#ifndef GMX_SIMD_IMPL_X86_AVX_128_FMA_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_128_FMA_SIMD4_FLOAT_H

#include "config.h"

#include <immintrin.h>
#include <x86intrin.h>

#include "gromacs/simd/impl_x86_sse4_1/impl_x86_sse4_1_simd4_float.h"

namespace gmx
{

static inline Simd4Float gmx_simdcall
fma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               _mm_macc_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
fms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               _mm_msub_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
fnma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               _mm_nmacc_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Float gmx_simdcall
fnms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return {
               _mm_nmsub_ps(a.simdInternal_, b.simdInternal_, c.simdInternal_)
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

#endif // GMX_SIMD_IMPL_X86_AVX_128_FMA_SIMD4_FLOAT_H
