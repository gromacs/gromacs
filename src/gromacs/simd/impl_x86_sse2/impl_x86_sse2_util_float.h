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

#ifndef GMX_SIMD_IMPL_X86_SSE2_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_X86_SSE2_UTIL_FLOAT_H

#include "config.h"

#include <emmintrin.h>

#include "impl_x86_sse2_common.h"

/* Higher-level utility functions (single) */
#define gmx_simd_gather_load_transpose_f             gmx_simd_gather_load_transpose_f_sse2
#define gmx_simd_best_pair_alignment_f               gmx_simd_best_pair_alignment_f_sse2
#define gmx_simd_gather_loadu_transpose_f            gmx_simd_gather_loadu_transpose_f_sse2
#define gmx_simd_transpose_scatter_storeu_f          gmx_simd_transpose_scatter_storeu_f_sse2
#define gmx_simd_transpose_scatter_incru_f           gmx_simd_transpose_scatter_incru_f_sse2
#define gmx_simd_transpose_scatter_decru_f           gmx_simd_transpose_scatter_decru_f_sse2
#define gmx_simd_expand_scalars_to_triplets_f        gmx_simd_expand_scalars_to_triplets_f_sse2
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_sse2
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_sse2
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_sse2

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus
template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_sse2(const float *        base,
                                      const gmx_int32_t    offset[],
                                      gmx_simd_float_t    &v0,
                                      gmx_simd_float_t    &v1,
                                      gmx_simd_float_t    &v2,
                                      gmx_simd_float_t    &v3)
{
    /* This will work with any value for align, as long as the resulting pointer is aligned */
    v0  = _mm_load_ps( base + align * offset[0] );
    v1  = _mm_load_ps( base + align * offset[1] );
    v2  = _mm_load_ps( base + align * offset[2] );
    v3  = _mm_load_ps( base + align * offset[3] );
    _MM_TRANSPOSE4_PS(v0, v1, v2, v3);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_sse2(const float *        base,
                                      const gmx_int32_t    offset[],
                                      gmx_simd_float_t    &v0,
                                      gmx_simd_float_t    &v1)
{
    __m128 t1, t2;
    /* This will work with any value for align, as long as the resulting pointer is aligned */
    v0  = _mm_castpd_ps(_mm_load_sd( (const double *)( base + align * offset[0] ) ));
    v1  = _mm_castpd_ps(_mm_load_sd( (const double *)( base + align * offset[1] ) ));
    t1  = _mm_castpd_ps(_mm_load_sd( (const double *)( base + align * offset[2] ) ));
    t2  = _mm_castpd_ps(_mm_load_sd( (const double *)( base + align * offset[3] ) ));
    t1  = _mm_unpacklo_ps(v0, t1);
    t2  = _mm_unpacklo_ps(v1, t2);
    v0  = _mm_unpacklo_ps(t1, t2);
    v1  = _mm_unpackhi_ps(t1, t2);
}

static const int gmx_simd_best_pair_alignment_f_sse2 = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_f_sse2(const float *        base,
                                       const gmx_int32_t    offset[],
                                       gmx_simd_float_t    &v0,
                                       gmx_simd_float_t    &v1,
                                       gmx_simd_float_t    &v2)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;
    
    /* The template conditional should be optimized away at compile time */
    if ( (align & 0x3) == 0)
    {
        /* With alignment 4 or better we can read a byte beyond triplets and use _mm_loadu_ps() */
        t1  = _mm_loadu_ps( base + align * offset[0] );
        t2  = _mm_loadu_ps( base + align * offset[1] );
        t3  = _mm_loadu_ps( base + align * offset[2] );
        t4  = _mm_loadu_ps( base + align * offset[3] );
        t5  = _mm_unpacklo_ps(t1, t2);
        t6  = _mm_unpacklo_ps(t3, t4);
        t7  = _mm_unpackhi_ps(t1, t2);
        t8  = _mm_unpackhi_ps(t3, t4);
        v0  = _mm_movelh_ps(t5, t6);
        v1  = _mm_movehl_ps(t6, t5);
        v2  = _mm_movelh_ps(t7, t8);
    }
    else
    {
        t1  = _mm_castpd_ps(_mm_load_sd( (const double *)( base + align * offset[0] ) ));
        t2  = _mm_castpd_ps(_mm_load_sd( (const double *)( base + align * offset[1] ) ));
        t3  = _mm_castpd_ps(_mm_load_sd( (const double *)( base + align * offset[2] ) ));
        t4  = _mm_castpd_ps(_mm_load_sd( (const double *)( base + align * offset[3] ) ));
        t5  = _mm_load_ss( base + align * offset[0] + 2 );
        t6  = _mm_load_ss( base + align * offset[1] + 2 );
        t7  = _mm_load_ss( base + align * offset[2] + 2 );
        t8  = _mm_load_ss( base + align * offset[3] + 2 );
        t1  = _mm_unpacklo_ps(t1, t2);
        t3  = _mm_unpacklo_ps(t3, t4);
        v0  = _mm_movelh_ps(t1, t3);
        v1  = _mm_movehl_ps(t3, t1);
        t5  = _mm_unpacklo_ps(t5, t6);
        t7  = _mm_unpacklo_ps(t7, t8);
        v2  = _mm_movelh_ps(t5, t7);
    }
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_f_sse2(float *              base,
                                         const gmx_int32_t    offset[],
                                         gmx_simd_float_t     v0,
                                         gmx_simd_float_t     v1,
                                         gmx_simd_float_t     v2)
{
    __m128 t1, t2;
    
    t1   = _mm_unpacklo_ps(v0, v1);
    t2   = _mm_unpackhi_ps(v0, v1);
    _mm_storel_pi( (__m64 *)( base + align * offset[0] ), t1);
    _mm_store_ss(base + align * offset[0] + 2, v2);
    _mm_storeh_pi( (__m64 *)( base + align * offset[1] ), t1);
    _mm_store_ss(base + align * offset[1] + 2, _mm_shuffle_ps(v2, v2, _MM_SHUFFLE(1, 1, 1, 1)));
    _mm_storel_pi( (__m64 *)( base + align * offset[2] ), t2);
    _mm_store_ss(base + align * offset[2] + 2, _mm_shuffle_ps(v2, v2, _MM_SHUFFLE(2, 2, 2, 2)));
    _mm_storeh_pi( (__m64 *)( base + align * offset[3] ), t2);
    _mm_store_ss(base + align * offset[3] + 2, _mm_shuffle_ps(v2, v2, _MM_SHUFFLE(3, 3, 3, 3)));
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_f_sse2(float *              base,
                                        const gmx_int32_t    offset[],
                                        gmx_simd_float_t     v0,
                                        gmx_simd_float_t     v1,
                                        gmx_simd_float_t     v2)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    t5          = _mm_unpacklo_ps(v1, v2);
    t6          = _mm_unpackhi_ps(v1, v2);
    t7          = _mm_shuffle_ps(v0, t5, _MM_SHUFFLE(1, 0, 0, 0));
    t8          = _mm_shuffle_ps(v0, t5, _MM_SHUFFLE(3, 2, 0, 1));
    t9          = _mm_shuffle_ps(v0, t6, _MM_SHUFFLE(1, 0, 0, 2));
    t10         = _mm_shuffle_ps(v0, t6, _MM_SHUFFLE(3, 2, 0, 3));
    t1          = _mm_load_ss(base + align * offset[0]);
    t1          = _mm_loadh_pi(t1, (__m64 *)(base + align * offset[0] + 1));
    t1          = _mm_add_ps(t1, t7);
    _mm_store_ss(base + align * offset[0], t1);
    _mm_storeh_pi( (__m64 *)(base + align * offset[0] + 1), t1);
    t2          = _mm_load_ss(base + align * offset[1]);
    t2          = _mm_loadh_pi(t2, (__m64 *)(base + align * offset[1] + 1));
    t2          = _mm_add_ps(t2, t8);
    _mm_store_ss(base + align * offset[1], t2);
    _mm_storeh_pi((__m64 *)(base + align * offset[1] + 1), t2);
    t3          = _mm_load_ss(base + align * offset[2]);
    t3          = _mm_loadh_pi(t3, (__m64 *)(base + align * offset[2] + 1));
    t3          = _mm_add_ps(t3, t9);
    _mm_store_ss(base + align * offset[2], t3);
    _mm_storeh_pi((__m64 *)(base + align * offset[2] + 1), t3);
    t4          = _mm_load_ss(base + align * offset[3]);
    t4          = _mm_loadh_pi(t4, (__m64 *)(base + align * offset[3] + 1));
    t4          = _mm_add_ps(t4, t10);
    _mm_store_ss(base + align * offset[3], t4);
    _mm_storeh_pi((__m64 *)(base + align * offset[3] + 1), t4);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_f_sse2(float *              base,
                                        const gmx_int32_t    offset[],
                                        gmx_simd_float_t     v0,
                                        gmx_simd_float_t     v1,
                                        gmx_simd_float_t     v2)
{
    /* This implementation is identical to the increment version, apart from using subtraction instead */
    __m128 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    t5          = _mm_unpacklo_ps(v1, v2);
    t6          = _mm_unpackhi_ps(v1, v2);
    t7          = _mm_shuffle_ps(v0, t5, _MM_SHUFFLE(1, 0, 0, 0));
    t8          = _mm_shuffle_ps(v0, t5, _MM_SHUFFLE(3, 2, 0, 1));
    t9          = _mm_shuffle_ps(v0, t6, _MM_SHUFFLE(1, 0, 0, 2));
    t10         = _mm_shuffle_ps(v0, t6, _MM_SHUFFLE(3, 2, 0, 3));
    t1          = _mm_load_ss(base + align * offset[0]);
    t1          = _mm_loadh_pi(t1, (__m64 *)(base + align * offset[0] + 1));
    t1          = _mm_sub_ps(t1, t7);
    _mm_store_ss(base + align * offset[0], t1);
    _mm_storeh_pi((__m64 *)(base + align * offset[0] + 1), t1);
    t2          = _mm_load_ss(base + align * offset[1]);
    t2          = _mm_loadh_pi(t2, (__m64 *)(base + align * offset[1] + 1));
    t2          = _mm_sub_ps(t2, t8);
    _mm_store_ss(base + align * offset[1], t2);
    _mm_storeh_pi((__m64 *)(base + align * offset[1] + 1), t2);
    t3          = _mm_load_ss(base + align * offset[2]);
    t3          = _mm_loadh_pi(t3, (__m64 *)(base + align * offset[2] + 1));
    t3          = _mm_sub_ps(t3, t9);
    _mm_store_ss(base + align * offset[2], t3);
    _mm_storeh_pi((__m64 *)(base + align * offset[2] + 1), t3);
    t4          = _mm_load_ss(base + align * offset[3]);
    t4          = _mm_loadh_pi(t4, (__m64 *)(base + align * offset[3] + 1));
    t4          = _mm_sub_ps(t4, t10);
    _mm_store_ss(base + align * offset[3], t4);
    _mm_storeh_pi((__m64 *)(base + align * offset[3] + 1), t4);
}

static gmx_inline void
gmx_simd_expand_scalars_to_triplets_f_sse2(gmx_simd_float_t    scalar,
                                           gmx_simd_float_t   &triplets0,
                                           gmx_simd_float_t   &triplets1,
                                           gmx_simd_float_t   &triplets2)
{
    triplets0 = _mm_shuffle_ps(scalar, scalar, _MM_SHUFFLE(1, 0, 0, 0));
    triplets1 = _mm_shuffle_ps(scalar, scalar, _MM_SHUFFLE(2, 2, 1, 1));
    triplets2 = _mm_shuffle_ps(scalar, scalar, _MM_SHUFFLE(3, 3, 3, 2));
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_sse2(const float *       base,
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
        v0  = _mm_load_ps(base + gmx_simd_extract_fi(offset, 0) );
        v1  = _mm_load_ps(base + gmx_simd_extract_fi(offset, 1) );
        v2  = _mm_load_ps(base + gmx_simd_extract_fi(offset, 2) );
        v3  = _mm_load_ps(base + gmx_simd_extract_fi(offset, 3) );
    }
    else
    {
        v0  = _mm_load_ps(base + align * gmx_simd_extract_fi(offset, 0) );
        v1  = _mm_load_ps(base + align * gmx_simd_extract_fi(offset, 1) );
        v2  = _mm_load_ps(base + align * gmx_simd_extract_fi(offset, 2) );
        v3  = _mm_load_ps(base + align * gmx_simd_extract_fi(offset, 3) );
    }
    _MM_TRANSPOSE4_PS(v0, v1, v2, v3);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_sse2(const float *       base,
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
        v0  = _mm_castpd_ps(_mm_load_sd((const double *)(base + gmx_simd_extract_fi(offset, 0))));
        v1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + gmx_simd_extract_fi(offset, 1))));
        t1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + gmx_simd_extract_fi(offset, 2))));
        t2  = _mm_castpd_ps(_mm_load_sd((const double *)(base + gmx_simd_extract_fi(offset, 3))));
    }
    else
    {
        v0  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * gmx_simd_extract_fi(offset, 0))));
        v1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * gmx_simd_extract_fi(offset, 1))));
        t1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * gmx_simd_extract_fi(offset, 2))));
        t2  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * gmx_simd_extract_fi(offset, 3))));
    }
    t1  = _mm_unpacklo_ps(v0, t1);
    t2  = _mm_unpacklo_ps(v1, t2);
    v0  = _mm_unpacklo_ps(t1, t2);
    v1  = _mm_unpackhi_ps(t1, t2);
}

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_f_sse2(const float *       base,
                                                 gmx_simd_fint32_t   offset,
                                                 gmx_simd_float_t   &v0,
                                                 gmx_simd_float_t   &v1)
{
    __m128 t1, t2;
    
    /* Use optimized bit-shift multiply for the most common alignments. */
    /* Do nothing for align == 1 */
    if (align == 2)
    {
        offset = _mm_slli_epi32(offset, 1);
    }
    else if (align == 4)
    {
        offset = _mm_slli_epi32(offset, 2);
    }
    
    if (align == 1 || align == 2 || align == 4)
    {
        v0  = _mm_castpd_ps(_mm_load_sd((const double *)(base + gmx_simd_extract_fi(offset, 0))));
        v1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + gmx_simd_extract_fi(offset, 1))));
        t1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + gmx_simd_extract_fi(offset, 2))));
        t2  = _mm_castpd_ps(_mm_load_sd((const double *)(base + gmx_simd_extract_fi(offset, 3))));
    }
    else
    {
        v0  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * gmx_simd_extract_fi(offset, 0))));
        v1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * gmx_simd_extract_fi(offset, 1))));
        t1  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * gmx_simd_extract_fi(offset, 2))));
        t2  = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * gmx_simd_extract_fi(offset, 3))));
    }
    t1  = _mm_unpacklo_ps(v0, t1);
    t2  = _mm_unpacklo_ps(v1, t2);
    v0  = _mm_unpacklo_ps(t1, t2);
    v1  = _mm_unpackhi_ps(t1, t2);
}


static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_f_sse2(float *           m,
                                         gmx_simd_float_t  v0,
                                         gmx_simd_float_t  v1,
                                         gmx_simd_float_t  v2,
                                         gmx_simd_float_t  v3)
{
    _MM_TRANSPOSE4_PS(v0, v1, v2, v3);
    v0 = _mm_add_ps(v0, v1);
    v2 = _mm_add_ps(v2, v3);
    v0 = _mm_add_ps(v0, v2);
    v2 = _mm_add_ps(v0, _mm_load_ps(m));
    _mm_store_ps(m, v2);
    
    return gmx_simd_reduce_f_sse2(v0);
}
#endif /* __cplusplus */

#endif /* GMX_SIMD_IMPL_X86_SSE2_UTIL_FLOAT_H */
