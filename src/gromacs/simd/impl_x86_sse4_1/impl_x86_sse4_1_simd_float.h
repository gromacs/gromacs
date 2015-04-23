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

#ifndef GMX_SIMD_IMPL_X86_SSE4_1_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_SSE4_1_SIMD_FLOAT_H

#include "config.h"

#include <smmintrin.h>

#include "gromacs/simd/impl_x86_sse2/impl_x86_sse2_simd_float.h"

namespace gmx
{

template<int index> gmx_simdcall
static inline std::int32_t
simdExtractFI(SimdFInt32 a)
{
    return _mm_extract_epi32(a.i, index);
}

static inline SimdFloat
simdRsqrtMaskF(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.r = _mm_blendv_ps(_mm_set1_ps(1.0f), x.r, m.b);
#endif
    return {
               _mm_and_ps(_mm_rsqrt_ps(x.r), m.b)
    };
}

static inline SimdFloat
simdRcpMaskF(SimdFloat x, SimdFBool m)
{
#ifndef NDEBUG
    x.r = _mm_blendv_ps(_mm_set1_ps(1.0f), x.r, m.b);
#endif
    return {
               _mm_and_ps(_mm_rcp_ps(x.r), m.b)
    };
}

static inline SimdFloat gmx_simdcall
simdRoundF(SimdFloat x)
{
    return {
               _mm_round_ps(x.r, _MM_FROUND_NINT)
    };
}

static inline SimdFloat gmx_simdcall
simdTruncF(SimdFloat x)
{
    return {
               _mm_round_ps(x.r, _MM_FROUND_TRUNC)
    };
}

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE4_1
static inline SimdFloat gmx_simdcall
simdFractionF(SimdFloat x)
{
    return {
               _mm_sub_ps(x.r, _mm_round_ps(x.r, _MM_FROUND_TRUNC))
    };
}
#endif

static inline SimdFloat gmx_simdcall
simdBlendF(SimdFloat a, SimdFloat b, SimdFBool sel)
{
    return {
               _mm_blendv_ps(a.r, b.r, sel.b)
    };
}

static inline SimdFInt32 gmx_simdcall
simdMulFI(SimdFInt32 a, SimdFInt32 b)
{
    return {
               _mm_mullo_epi32(a.i, b.i)
    };
}

static inline SimdFInt32 gmx_simdcall
simdBlendFI(SimdFInt32 a, SimdFInt32 b, SimdFIBool sel)
{
    return {
               _mm_blendv_epi8(a.i, b.i, sel.b)
    };
}


}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE4_1_SIMD_FLOAT_H
