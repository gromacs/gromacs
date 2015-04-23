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

#ifndef GMX_SIMD_IMPL_X86_SSE4_1_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_X86_SSE4_1_UTIL_FLOAT_H

#include "config.h"

#include <smmintrin.h>

#include "impl_x86_sse4_1_common.h"

#undef  gmx_simd_gather_load_bysimdint_transpose_f
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_sse4_1
#undef  gmx_simd_reduce_incr_4_return_sum_f
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_sse4_1

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus
template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_f_sse4_1(const float *       base,
                                                  gmx_simd_fint32_t   offset,
                                                  gmx_simd_float_t   &v0,
                                                  gmx_simd_float_t   &v1,
                                                  gmx_simd_float_t   &v2,
                                                  gmx_simd_float_t   &v3)
{
    /* Use optimized bit-shift multiply for the most common alignments */
    if (align == 4)
    {
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 8)
    {
        offset = _mm_slli_epi32(offset, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset = _mm_add_epi32(offset, _mm_slli_epi32(offset, 1));
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 16)
    {
        offset = _mm_slli_epi32(offset, 4);
    }

    if (align == 4 || align == 8 || align == 12 || align == 16)
    {
        v0  = _mm_load_ps(base + _mm_extract_epi32(offset, 0) );
        v1  = _mm_load_ps(base + _mm_extract_epi32(offset, 1) );
        v2  = _mm_load_ps(base + _mm_extract_epi32(offset, 2) );
        v3  = _mm_load_ps(base + _mm_extract_epi32(offset, 3) );
    }
    else
    {
        v0  = _mm_load_ps(base + align * _mm_extract_epi32(offset, 0) );
        v1  = _mm_load_ps(base + align * _mm_extract_epi32(offset, 1) );
        v2  = _mm_load_ps(base + align * _mm_extract_epi32(offset, 2) );
        v3  = _mm_load_ps(base + align * _mm_extract_epi32(offset, 3) );
    }
    _MM_TRANSPOSE4_PS(v0, v1, v2, v3);
}


template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_f_sse4_1(const float *       base,
                                                  gmx_simd_fint32_t   offset,
                                                  gmx_simd_float_t   &v0,
                                                  gmx_simd_float_t   &v1)
{
    __m128 t1, t2;

    /* Use optimized bit-shift multiply for the most common alignments */
    if (align == 2)
    {
        offset = _mm_slli_epi32(offset, 1);
    }
    else if (align == 4)
    {
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 6)
    {
        /* multiply by 3, then by 2 */
        offset = _mm_add_epi32(offset, _mm_slli_epi32(offset, 1));
        offset = _mm_slli_epi32(offset, 1);
    }
    else if (align == 8)
    {
        offset = _mm_slli_epi32(offset, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset = _mm_add_epi32(offset, _mm_slli_epi32(offset, 1));
        offset = _mm_slli_epi32(offset, 2);
    }
    else if (align == 16)
    {
        offset = _mm_slli_epi32(offset, 4);
    }

    if (align == 2 || align == 4 || align == 6 ||
        align == 8 || align == 12 || align == 16)
    {
        v0  = _mm_castpd_ps(_mm_load_sd((const double *)(base + _mm_extract_epi32(offset, 0))));
        v1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + _mm_extract_epi32(offset, 1))));
        t1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + _mm_extract_epi32(offset, 2))));
        t2  = _mm_castpd_ps(_mm_load_sd((const double *)(base + _mm_extract_epi32(offset, 3))));
    }
    else
    {
        v0  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * _mm_extract_epi32(offset, 0))));
        v1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * _mm_extract_epi32(offset, 1))));
        t1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * _mm_extract_epi32(offset, 2))));
        t2  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * _mm_extract_epi32(offset, 3))));
    }
    t1  = _mm_unpacklo_ps(v0, t1);
    t2  = _mm_unpacklo_ps(v1, t2);
    v0  = _mm_unpacklo_ps(t1, t2);
    v1  = _mm_unpackhi_ps(t1, t2);
}

static gmx_inline float gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_f_sse4_1(float *           m,
                                           gmx_simd_float_t  v0,
                                           gmx_simd_float_t  v1,
                                           gmx_simd_float_t  v2,
                                           gmx_simd_float_t  v3)
{
    /* The present implementation is quite standardized: transpose and sum.
     * However, there are a couple of alternatives worth trying and benchmarking later:
     *
     * 1) Use 3*hadd. Much fewer instructions, but high latency (5) and bad throughput (2).
     * 2) Mix hadd with shuffle/blend (blend is better t-put than shuffle).
     *
     * In a standalone test all these had similar performance, so we need to test what
     * works best in the actual kernels where it is mixed with other instructions.
     */
    _MM_TRANSPOSE4_PS(v0, v1, v2, v3);
    v0 = _mm_add_ps(v0, v1);
    v2 = _mm_add_ps(v2, v3);
    v0 = _mm_add_ps(v0, v2);
    v2 = _mm_add_ps(v0, _mm_load_ps(m));
    _mm_store_ps(m, v2);

    return gmx_simd_reduce_f_sse4_1(v0);
}
#endif

#endif /* GMX_SIMD_IMPL_X86_SSE4_1_UTIL_FLOAT_H */
