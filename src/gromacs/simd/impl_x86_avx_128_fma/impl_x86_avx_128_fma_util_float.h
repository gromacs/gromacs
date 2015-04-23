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

#ifndef GMX_SIMD_IMPL_X86_AVX_128_FMA_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_128_FMA_UTIL_FLOAT_H

#include "config.h"

#include <immintrin.h>
#include <x86intrin.h>

#include "impl_x86_avx_128_fma_common.h"

#undef  gmx_simd_gather_loadu_transpose_f
#define gmx_simd_gather_loadu_transpose_f       gmx_simd_gather_loadu_transpose_f_avx_128_fma
#undef  gmx_simd_transpose_scatter_storeu_f
#define gmx_simd_transpose_scatter_storeu_f     gmx_simd_transpose_scatter_storeu_f_avx_128_fma
#undef  gmx_simd_expand_scalars_to_triplets_f
#define gmx_simd_expand_scalars_to_triplets_f   gmx_simd_expand_scalars_to_triplets_f_avx_128_fma

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus
template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_f_avx_128_fma(const float *        base,
                                              const gmx_int32_t    offset[],
                                              gmx_simd_float_t    &v0,
                                              gmx_simd_float_t    &v1,
                                              gmx_simd_float_t    &v2)
{
    __m128  t3;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);
    
    if ( (align & 0x3) == 0)
    {
        v0  = _mm_loadu_ps( base + align * offset[0] );
        v1  = _mm_loadu_ps( base + align * offset[1] );
        v2  = _mm_loadu_ps( base + align * offset[2] );
        t3  = _mm_loadu_ps( base + align * offset[3] );
    }
    else
    {
        v0  = gmx_mm_maskload_ps(base + align * offset[0], mask);
        v1  = gmx_mm_maskload_ps(base + align * offset[1], mask);
        v2  = gmx_mm_maskload_ps(base + align * offset[2], mask);
        t3  = gmx_mm_maskload_ps(base + align * offset[3], mask);
    }
    _MM_TRANSPOSE4_PS(v0, v1, v2, t3);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_f_avx_128_fma(float *              base,
                                                const gmx_int32_t    offset[],
                                                gmx_simd_float_t     v0,
                                                gmx_simd_float_t     v1,
                                                gmx_simd_float_t     v2)
{
    __m128  t3   = _mm_setzero_ps();
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);
    _MM_TRANSPOSE4_PS(v0, v1, v2, t3);
    gmx_mm_maskstore_ps(base + align * offset[0], mask, v0);
    gmx_mm_maskstore_ps(base + align * offset[1], mask, v1);
    gmx_mm_maskstore_ps(base + align * offset[2], mask, v2);
    gmx_mm_maskstore_ps(base + align * offset[3], mask, t3);
    
}

/* In the old group kernels we found it more efficient to transpose the data to store rather
 * than using maskload and maskstore. It might be worth to test again, but for now we assume
 * this is still the case, and rely on those version inherited from the SSE2 header.
 *
 * It is also worth testing if changing _mm_shuffle_ps() to _mm_permute_ps() could improve
 * throughput just-so-slightly.
 */

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_f_avx_128_fma(gmx_simd_float_t   scalar,
                                                  gmx_simd_float_t  &triplets0,
                                                  gmx_simd_float_t  &triplets1,
                                                  gmx_simd_float_t  &triplets2)
{
    triplets0 = _mm_permute_ps(scalar, _MM_SHUFFLE(1, 0, 0, 0));
    triplets1 = _mm_permute_ps(scalar, _MM_SHUFFLE(2, 2, 1, 1));
    triplets2 = _mm_permute_ps(scalar, _MM_SHUFFLE(3, 3, 3, 2));
}

/* Permute might be a small help for gmx_simd_reduce_incr_4_return_sum_f(). Test in kernels later. */
#endif /* __cplusplus */


#endif /* GMX_SIMD_IMPL_X86_AVX_128_FMA_UTIL_FLOAT_H */
