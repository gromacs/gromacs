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

#ifndef GMX_SIMD_IMPL_X86_SSE4_1_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_X86_SSE4_1_UTIL_DOUBLE_H

#include "config.h"

#include <smmintrin.h>

#include "impl_x86_sse4_1_common.h"

#undef  gmx_simd_gather_load_bysimdint_transpose_d
#define gmx_simd_gather_load_bysimdint_transpose_d   gmx_simd_gather_load_bysimdint_transpose_d_sse4_1
#undef  gmx_simd_reduce_incr_4_return_sum_d
#define gmx_simd_reduce_incr_4_return_sum_d          gmx_simd_reduce_incr_4_return_sum_d_sse4_1

/****************************************************
 * Double precision higher-level utility functions  *
 ****************************************************/



#ifdef __cplusplus
template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_d_sse4_1(const double *       base,
                                                  gmx_simd_dint32_t    offset,
                                                  gmx_simd_double_t   &v0,
                                                  gmx_simd_double_t   &v1,
                                                  gmx_simd_double_t   &v2,
                                                  gmx_simd_double_t   &v3)
{
    __m128d t1, t2, t3, t4;
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
        t1  = _mm_load_pd(base + _mm_extract_epi32(offset, 0));
        t2  = _mm_load_pd(base + _mm_extract_epi32(offset, 1));
        t3  = _mm_load_pd(base + _mm_extract_epi32(offset, 0) + 2);
        t4  = _mm_load_pd(base + _mm_extract_epi32(offset, 1) + 2);
    }
    else
    {
        t1  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 0));
        t2  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 1));
        t3  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 0) + 2);
        t4  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 1) + 2);
    }
    v0  = _mm_unpacklo_pd(t1, t2);
    v1  = _mm_unpackhi_pd(t1, t2);
    v2  = _mm_unpacklo_pd(t3, t4);
    v3  = _mm_unpackhi_pd(t3, t4);
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_d_sse4_1(const double *       base,
                                                  gmx_simd_dint32_t    offset,
                                                  gmx_simd_double_t   &v0,
                                                  gmx_simd_double_t   &v1)
{
    __m128d t1, t2;
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
        t1  = _mm_load_pd(base + _mm_extract_epi32(offset, 0));
        t2  = _mm_load_pd(base + _mm_extract_epi32(offset, 1));
    }
    else
    {
        t1  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 0));
        t2  = _mm_load_pd(base + align * _mm_extract_epi32(offset, 1));
    }
    v0  = _mm_unpacklo_pd(t1, t2);
    v1  = _mm_unpackhi_pd(t1, t2);
}

static gmx_inline double gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_d_sse4_1(double *           m,
                                           gmx_simd_double_t  v0,
                                           gmx_simd_double_t  v1,
                                           gmx_simd_double_t  v2,
                                           gmx_simd_double_t  v3)
{
    __m128d t1, t2, t3, t4;

    t1 = _mm_unpacklo_pd(v0, v1);
    t2 = _mm_unpackhi_pd(v0, v1);
    t3 = _mm_unpacklo_pd(v2, v3);
    t4 = _mm_unpackhi_pd(v2, v3);

    t1 = _mm_add_pd(t1, t2);
    t3 = _mm_add_pd(t3, t4);

    t2 = _mm_add_pd(t1, _mm_load_pd(m));
    t4 = _mm_add_pd(t3, _mm_load_pd(m + 2));
    _mm_store_pd(m, t2);
    _mm_store_pd(m + 2, t4);

    t1 = _mm_add_pd(t1, t3);
    return gmx_simd_reduce_d_sse4_1(t1);
}
#endif

#endif /* GMX_SIMD_IMPL_X86_SSE4_1_UTIL_DOUBLE_H */
