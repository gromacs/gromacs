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

#ifndef GMX_SIMD_IMPL_X86_SSE4_1_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_X86_SSE4_1_SIMD_DOUBLE_H

#include "config.h"

#include <smmintrin.h>

#include "impl_x86_sse4_1_common.h"

/* Almost all SSE4.1 instructions already exist in SSE2, but a few of them
 * can be implemented more efficiently in SSE4.1.
 */
#undef  simdRoundD
#define simdRoundD(x)       _mm_round_pd(x, _MM_FROUND_NINT)

#undef  simdTruncD
#define simdTruncD(x)       _mm_round_pd(x, _MM_FROUND_TRUNC)

#undef  simdExtractDI
#define simdExtractDI       _mm_extract_epi32

#undef  simdMulDI
#define simdMulDI           _mm_mullo_epi32

#undef  simdBlendD
#define simdBlendD         _mm_blendv_pd

#undef  simdReduceD
#define simdReduceD(a)      simdReduceD_sse4_1(a)

#undef  simdBlendDI
#define simdBlendDI        _mm_blendv_epi8

static inline double gmx_simdcall
simdReduceD_sse4_1(__m128d a)
{
    double  f;

    a = _mm_hadd_pd(a, a);
    _mm_store_sd(&f, a);
    return f;
}

#endif /* GMX_SIMD_IMPL_X86_SSE4_1_SIMD_DOUBLE_H */
