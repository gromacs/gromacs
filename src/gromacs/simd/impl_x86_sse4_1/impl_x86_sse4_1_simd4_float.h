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

#ifndef GMX_SIMD_IMPL_X86_SSE4_1_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_X86_SSE4_1_SIMD4_FLOAT_H

#include "config.h"

#include <smmintrin.h>

#include "gromacs/simd/impl_x86_sse2/impl_x86_sse2_simd4_float.h"

namespace gmx
{

static inline Simd4Float gmx_simdcall round(Simd4Float x)
{
    return { _mm_round_ps(x.simdInternal_, _MM_FROUND_NINT) };
}

static inline Simd4Float gmx_simdcall trunc(Simd4Float x)
{
    return { _mm_round_ps(x.simdInternal_, _MM_FROUND_TRUNC) };
}

static inline float gmx_simdcall dotProduct(Simd4Float a, Simd4Float b)
{
    __m128 res = _mm_dp_ps(a.simdInternal_, b.simdInternal_, 0x71);
    return *reinterpret_cast<float*>(&res);
}

static inline Simd4Float gmx_simdcall blend(Simd4Float a, Simd4Float b, Simd4FBool sel)
{
    return { _mm_blendv_ps(a.simdInternal_, b.simdInternal_, sel.simdInternal_) };
}

} // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE4_1_SIMD4_FLOAT_H
