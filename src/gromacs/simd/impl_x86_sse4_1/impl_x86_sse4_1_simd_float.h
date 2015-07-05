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

#include "impl_x86_sse4_1_common.h"

/* Almost all SSE4.1 instructions already exist in SSE2, but a few of them
 * can be implemented more efficiently in SSE4.1.
 */
#undef  gmx_simd_round_f
#define gmx_simd_round_f(x)       _mm_round_ps(x, _MM_FROUND_NINT)

#undef  gmx_simd_trunc_f
#define gmx_simd_trunc_f(x)       _mm_round_ps(x, _MM_FROUND_TRUNC)

#undef  gmx_simd_extract_fi
#define gmx_simd_extract_fi       _mm_extract_epi32

#undef  gmx_simd_mul_fi
#define gmx_simd_mul_fi           _mm_mullo_epi32

#undef  gmx_simd_blendv_f
#define gmx_simd_blendv_f         _mm_blendv_ps

#undef  gmx_simd_reduce_f
#define gmx_simd_reduce_f(a)      gmx_simd_reduce_f_sse4_1(a)

#undef  gmx_simd_blendv_fi
#define gmx_simd_blendv_fi        _mm_blendv_epi8

/* SIMD reduction function */
static gmx_inline float gmx_simdcall
gmx_simd_reduce_f_sse4_1(__m128 a)
{
    float  f;

    a = _mm_hadd_ps(a, a);
    a = _mm_hadd_ps(a, a);
    _mm_store_ss(&f, a);
    return f;
}

#endif /* GMX_SIMD_IMPL_X86_SSE4_1_SIMD_FLOAT_H */
