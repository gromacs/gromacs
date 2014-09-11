/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_X86_SSE4_1_H
#define GMX_SIMD_IMPL_X86_SSE4_1_H

#include "config.h"

#include <math.h>

#include <smmintrin.h>

/* x86 SSE4.1 SIMD instruction wrappers
 *
 * Please see documentation in gromacs/simd/simd.h for the available
 * defines.
 */

/* Inherit most of SSE4.1 from SSE2 */
#include "gromacs/simd/impl_x86_sse2/impl_x86_sse2.h"
/* Increment over SSE2 capabilities */
#define GMX_SIMD_X86_SSE4_1_OR_HIGHER


/* Override capability definitions from SSE2 */
#define  GMX_SIMD4_HAVE_FLOAT_DOTPRODUCT3

/* Almost all SSE4.1 instructions already exist in SSE2, but a few of them
 * can be implemented more efficiently in SSE4.1.
 */
#undef  gmx_simd_round_f
#define gmx_simd_round_f(x)       _mm_round_ps(x, _MM_FROUND_NINT)
#undef  gmx_simd_trunc_f
#define gmx_simd_trunc_f(x)       _mm_round_ps(x, _MM_FROUND_TRUNC)
#undef  gmx_simd_round_d
#define gmx_simd_round_d(x)       _mm_round_pd(x, _MM_FROUND_NINT)
#undef  gmx_simd_trunc_d
#define gmx_simd_trunc_d(x)       _mm_round_pd(x, _MM_FROUND_TRUNC)

#undef  gmx_simd_extract_fi
#define gmx_simd_extract_fi       _mm_extract_epi32
#undef  gmx_simd_mul_fi
#define gmx_simd_mul_fi           _mm_mullo_epi32

#undef  gmx_simd_extract_di
#define gmx_simd_extract_di       _mm_extract_epi32
#undef  gmx_simd_mul_di
#define gmx_simd_mul_di           _mm_mullo_epi32

#undef  gmx_simd_blendv_f
#define gmx_simd_blendv_f         _mm_blendv_ps
#undef  gmx_simd_blendv_d
#define gmx_simd_blendv_d         _mm_blendv_pd

#undef  gmx_simd_reduce_f
#define gmx_simd_reduce_f(a)      gmx_simd_reduce_f_sse4_1(a)
#undef  gmx_simd_reduce_d
#define gmx_simd_reduce_d(a)      gmx_simd_reduce_d_sse4_1(a)

#undef  gmx_simd_blendv_fi
#define gmx_simd_blendv_fi        _mm_blendv_epi8
#undef  gmx_simd_blendv_di
#define gmx_simd_blendv_di        _mm_blendv_epi8

#undef  gmx_simd4_dotproduct3_f
#define gmx_simd4_dotproduct3_f   gmx_simd4_dotproduct3_f_sse4_1

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

/* SIMD4 Dotproduct helper function */
static gmx_inline float gmx_simdcall
gmx_simd4_dotproduct3_f_sse4_1(__m128 a, __m128 b)
{
    float f;
    _MM_EXTRACT_FLOAT(f, _mm_dp_ps(a, b, 0x71), 0);
    return f;
}

static gmx_inline double gmx_simdcall
gmx_simd_reduce_d_sse4_1(__m128d a)
{
    double  f;

    a = _mm_hadd_pd(a, a);
    _mm_store_sd(&f, a);
    return f;
}

#endif /* GMX_SIMD_IMPL_X86_SSE4_1_H */
