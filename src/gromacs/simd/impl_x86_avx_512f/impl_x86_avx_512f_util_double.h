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

#ifndef GMX_SIMD_IMPL_X86_AVX_512F_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_512F_UTIL_DOUBLE_H

#include "config.h"

#include <assert.h>
#include <immintrin.h>
#include <math.h>

#include "gromacs/utility/real.h"

#include "impl_x86_avx_512f_common.h"

/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_d             gmx_simd_gather_load_transpose_d_x86_avx_512f
#define gmx_simd_best_pair_alignment_d               gmx_simd_best_pair_alignment_d_x86_avx_512f
#define gmx_simd_gather_loadu_transpose_d            gmx_simd_gather_loadu_transpose_d_x86_avx_512f
#define gmx_simd_transpose_scatter_storeu_d          gmx_simd_transpose_scatter_storeu_d_x86_avx_512f
#define gmx_simd_transpose_scatter_incru_d           gmx_simd_transpose_scatter_incru_d_x86_avx_512f
#define gmx_simd_transpose_scatter_decru_d           gmx_simd_transpose_scatter_decru_d_x86_avx_512f
#define gmx_simd_expand_scalars_to_triplets_d        gmx_simd_expand_scalars_to_triplets_d_x86_avx_512f
#define gmx_simd_gather_load_bysimdint_transpose_d   gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f
#define gmx_simd_gather_loadu_bysimdint_transpose_d  gmx_simd_gather_loadu_bysimdint_transpose_d_x86_avx_512f
#define gmx_simd_reduce_incr_4_return_sum_d          gmx_simd_reduce_incr_4_return_sum_d_x86_avx_512f
/* Half-simd-width utilities */
#define gmx_simd_load_dual_hsimd_d                   gmx_simd_load_dual_hsimd_d_x86_avx_512f
#define gmx_simd_loaddup_hsimd_d                     gmx_simd_loaddup_hsimd_d_x86_avx_512f
#define gmx_simd_load1_dual_hsimd_d                  gmx_simd_load1_dual_hsimd_d_x86_avx_512f
#define gmx_simd_store_dual_hsimd_d                  gmx_simd_store_dual_hsimd_d_x86_avx_512f
#define gmx_simd_decr_hsimd_d                        gmx_simd_decr_hsimd_d_x86_avx_512f
#define gmx_simd_gather_load_transpose_hsimd_d       gmx_simd_gather_load_transpose_hsimd_d_x86_avx_512f
#define gmx_simd_reduce_incr_4_return_sum_hsimd_d    gmx_simd_reduce_incr_4_return_sum_hsimd_d_x86_avx_512f

/****************************************************
 * Double precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

/* Declare functions implemented further down that are used here */
template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f(const double *       base,
                                                        gmx_simd_dint32_t    simdoffset,
                                                        gmx_simd_double_t   &v0,
                                                        gmx_simd_double_t   &v1,
                                                        gmx_simd_double_t   &v2,
                                                        gmx_simd_double_t   &v3);

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f(const double *       base,
                                                        gmx_simd_dint32_t    simdoffset,
                                                        gmx_simd_double_t   &v0,
                                                        gmx_simd_double_t   &v1);

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d_x86_avx_512f(const double *        base,
                                              const gmx_int32_t     offset[],
                                              gmx_simd_double_t    &v0,
                                              gmx_simd_double_t    &v1,
                                              gmx_simd_double_t    &v2,
                                              gmx_simd_double_t    &v3)
{
    assert((size_t)offset % 32 == 0);
    gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f<align>(base,
                                                                   gmx_simd_load_di(offset),
                                                                   v0, v1, v2, v3);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d_x86_avx_512f(const double *        base,
                                              const gmx_int32_t     offset[],
                                              gmx_simd_double_t    &v0,
                                              gmx_simd_double_t    &v1)
{
    assert((size_t)offset % 32 == 0);
    gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f<align>(base,
                                                                   gmx_simd_load_di(offset),
                                                                   v0, v1);
}

static const int gmx_simd_best_pair_alignment_d_x86_avx_512f = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_d_x86_avx_512f(const double *        base,
                                               const gmx_int32_t     offset[],
                                               gmx_simd_double_t    &v0,
                                               gmx_simd_double_t    &v1,
                                               gmx_simd_double_t    &v2)
{
    __m256i simdoffset;

    assert((size_t)offset % 32 == 0);

    simdoffset = gmx_simd_load_di(offset);

    if (align == 4)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm256_mullo_epi32(simdoffset, _mm256_set1_epi32(align));
    }

    v0  = _mm512_i32logather_pd(_mm512_castsi256_si512(simdoffset), base,   sizeof(double));
    v1  = _mm512_i32logather_pd(_mm512_castsi256_si512(simdoffset), base+1, sizeof(double));
    v2  = _mm512_i32logather_pd(_mm512_castsi256_si512(simdoffset), base+2, sizeof(double));
}



template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_d_x86_avx_512f(double *              base,
                                                 const gmx_int32_t     offset[],
                                                 gmx_simd_double_t     v0,
                                                 gmx_simd_double_t     v1,
                                                 gmx_simd_double_t     v2)
{
    __m256i simdoffset;

    assert((size_t)offset % 32 == 0);

    simdoffset = gmx_simd_load_di(offset);

    if (align == 4)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm256_mullo_epi32(simdoffset, _mm256_set1_epi32(align));
    }

    _mm512_i32loscatter_pd(base,   _mm512_castsi256_si512(simdoffset), v0, sizeof(double));
    _mm512_i32loscatter_pd(base+1, _mm512_castsi256_si512(simdoffset), v1, sizeof(double));
    _mm512_i32loscatter_pd(base+2, _mm512_castsi256_si512(simdoffset), v2, sizeof(double));
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_d_x86_avx_512f(double *              base,
                                                const gmx_int32_t     offset[],
                                                gmx_simd_double_t     v0,
                                                gmx_simd_double_t     v1,
                                                gmx_simd_double_t     v2)
{
    __m512d t0, t1, t2;
    gmx_simd_gather_loadu_transpose_d_x86_avx_512f<align>(base, offset, t0, t1, t2);
    t0 = _mm512_add_pd(t0, v0);
    t1 = _mm512_add_pd(t1, v1);
    t2 = _mm512_add_pd(t2, v2);
    gmx_simd_transpose_scatter_storeu_d_x86_avx_512f<align>(base, offset, t0, t1, t2);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_d_x86_avx_512f(double *              base,
                                                const gmx_int32_t     offset[],
                                                gmx_simd_double_t     v0,
                                                gmx_simd_double_t     v1,
                                                gmx_simd_double_t     v2)
{
    __m512d t0, t1, t2;
    gmx_simd_gather_loadu_transpose_d_x86_avx_512f<align>(base, offset, t0, t1, t2);
    t0 = _mm512_sub_pd(t0, v0);
    t1 = _mm512_sub_pd(t1, v1);
    t2 = _mm512_sub_pd(t2, v2);
    gmx_simd_transpose_scatter_storeu_d_x86_avx_512f<align>(base, offset, t0, t1, t2);
}


static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_d_x86_avx_512f(gmx_simd_double_t   scalar,
                                                   gmx_simd_double_t  &triplets0,
                                                   gmx_simd_double_t  &triplets1,
                                                   gmx_simd_double_t  &triplets2)
{
    triplets0 = _mm512_castsi512_pd(_mm512_permutevar_epi32(_mm512_set_epi32(5, 4, 5, 4, 3, 2, 3, 2, 3, 2, 1, 0, 1, 0, 1, 0),
                                                            _mm512_castpd_si512(scalar)));
    triplets1 = _mm512_castsi512_pd(_mm512_permutevar_epi32(_mm512_set_epi32(11, 10, 9, 8, 9, 8, 9, 8, 7, 6, 7, 6, 7, 6, 5, 4),
                                                            _mm512_castpd_si512(scalar)));
    triplets2 = _mm512_castsi512_pd(_mm512_permutevar_epi32(_mm512_set_epi32(15, 14, 15, 14, 15, 14, 13, 12, 13, 12, 13, 12, 11, 10, 11, 10),
                                                            _mm512_castpd_si512(scalar)));
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f(const double *       base,
                                                        gmx_simd_dint32_t    simdoffset,
                                                        gmx_simd_double_t   &v0,
                                                        gmx_simd_double_t   &v1,
                                                        gmx_simd_double_t   &v2,
                                                        gmx_simd_double_t   &v3)
{
    assert((size_t)base % 32 == 0);
    assert(align % 4 == 0);

    if (align == 4)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm256_mullo_epi32(simdoffset, _mm256_set1_epi32(align));
    }

    v0  = _mm512_i32logather_pd(_mm512_castsi256_si512(simdoffset), base,   sizeof(double));
    v1  = _mm512_i32logather_pd(_mm512_castsi256_si512(simdoffset), base+1, sizeof(double));
    v2  = _mm512_i32logather_pd(_mm512_castsi256_si512(simdoffset), base+2, sizeof(double));
    v3  = _mm512_i32logather_pd(_mm512_castsi256_si512(simdoffset), base+3, sizeof(double));
}


template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_d_x86_avx_512f(const double *       base,
                                                         gmx_simd_dint32_t    simdoffset,
                                                         gmx_simd_double_t   &v0,
                                                         gmx_simd_double_t   &v1)
{
    if (align == 2)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 1);
    }
    else if (align == 4)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm256_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm256_mullo_epi32(simdoffset, _mm256_set1_epi32(align));
    }

    v0  = _mm512_i32logather_pd(_mm512_castsi256_si512(simdoffset), base,   sizeof(double));
    v1  = _mm512_i32logather_pd(_mm512_castsi256_si512(simdoffset), base+1, sizeof(double));
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_x86_avx_512f(const double *       base,
                                                        gmx_simd_dint32_t    simdoffset,
                                                        gmx_simd_double_t   &v0,
                                                        gmx_simd_double_t   &v1)
{
    assert((size_t)base % 32 == 0);
    assert(align % 2 == 0);
    gmx_simd_gather_loadu_bysimdint_transpose_d_x86_avx_512f<align>(base, simdoffset, v0, v1);
}


static gmx_inline double gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_d_x86_avx_512f(double * m,
                                                 gmx_simd_double_t v0, gmx_simd_double_t v1,
                                                 gmx_simd_double_t v2, gmx_simd_double_t v3)
{
    double  d;
    __m512d t0, t2;

    assert((size_t)m % 32 == 0);

    /* Not optimal. Will optimize when we know more about AVX-512 latencies */
    t0 = _mm512_add_pd(v0, _mm512_permute_pd(v0, 0x55));
    t2 = _mm512_add_pd(v2, _mm512_permute_pd(v2, 0x55));
    t0 = _mm512_mask_add_pd(t0, _mm512_int2mask(0xAA), v1, _mm512_permute_pd(v1, 0x55));
    t2 = _mm512_mask_add_pd(t2, _mm512_int2mask(0xAA), v3, _mm512_permute_pd(v3, 0x55));
    t0 = _mm512_add_pd(t0, _mm512_shuffle_f64x2(t0, t0, 0x4E));
    t0 = _mm512_mask_add_pd(t0, _mm512_int2mask(0xF0), t2, _mm512_shuffle_f64x2(t2, t2, 0x4E));
    t0 = _mm512_add_pd(t0, _mm512_shuffle_f64x2(t0, t0, 0xB1));
    t0 = _mm512_mask_shuffle_f64x2(t0, _mm512_int2mask(0x0C), t0, t0, 0xEE);

    t2 = _mm512_maskz_loadu_pd(_mm512_int2mask(0xF), m);
    t2 = _mm512_add_pd(t2, t0);
    _mm512_mask_storeu_pd(m, _mm512_int2mask(0xF), t2);

    t0 = _mm512_add_pd(t0, _mm512_permutex_pd(t0, 0x4E));
    t0 = _mm512_add_pd(t0, _mm512_permutex_pd(t0, 0xB1));

    _mm512_mask_storeu_pd(&d, _mm512_int2mask(0x0001), t0);
    return d;
}

/******************************************************
 * Double precision half-simd-width utility functions *
 ******************************************************/
static gmx_inline gmx_simd_double_t
gmx_simd_load_dual_hsimd_d_x86_avx_512f(const double * m0,
                                        const double * m1)
{
    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    return _mm512_mask_loadu_pd(_mm512_maskz_loadu_pd(_mm512_int2mask(0x0F), m0), _mm512_int2mask(0xF0), m1-4);
}

static gmx_inline gmx_simd_double_t
gmx_simd_loaddup_hsimd_d_x86_avx_512f(const double * m)
{
    __m512d tmp;

    assert((size_t)m % 32 == 0);

    tmp = _mm512_maskz_loadu_pd(_mm512_int2mask(0x0F), m);
    return _mm512_shuffle_f64x2(tmp, tmp, 0x44);
}

static gmx_inline gmx_simd_double_t
gmx_simd_load1_dual_hsimd_d_x86_avx_512f(const double * m)
{
    __m512d tmp;
    tmp = _mm512_maskz_expandloadu_pd(_mm512_int2mask(0x11), m);
    return _mm512_permutex_pd(tmp, 0);
}


static gmx_inline void
gmx_simd_store_dual_hsimd_d_x86_avx_512f(double *           m0,
                                         double *           m1,
                                         gmx_simd_double_t  a)
{
    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    _mm512_mask_storeu_pd(m0,   _mm512_int2mask(0x0F), a);
    _mm512_mask_storeu_pd(m1-4, _mm512_int2mask(0xF0), a);
}

static gmx_inline void
gmx_simd_decr_hsimd_d_x86_avx_512f(double *          m,
                                   gmx_simd_double_t a)
{
    __m512d t;

    assert((size_t)m % 32 == 0);

    a = _mm512_add_pd(a, _mm512_shuffle_f64x2(a, a, 0xEE));
    t = _mm512_maskz_loadu_pd(_mm512_int2mask(0x0F), m);
    t = _mm512_sub_pd(t, a);
    _mm512_mask_storeu_pd(m, _mm512_int2mask(0x0F), t);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_hsimd_d_x86_avx_512f(const double *       base0,
                                                    const double *       base1,
                                                    gmx_int32_t          offset[],
                                                    gmx_simd_double_t   &v0,
                                                    gmx_simd_double_t   &v1)
{
    __m512i  idx0, idx1, idx;
    __m512d  tmp1, tmp2;

    assert((size_t)offset % 16 == 0);
    assert((size_t)base0 % 16 == 0);
    assert((size_t)base1 % 16 == 0);
    assert((size_t)align % 2  == 0);

    idx0 = _mm512_maskz_loadu_epi32(_mm512_int2mask(0x000F), offset);

    idx0 = _mm512_mullo_epi32(idx0, _mm512_set1_epi32(align));
    idx1 = _mm512_add_epi32(idx0, _mm512_set1_epi32(1));

    idx = _mm512_mask_shuffle_i32x4(idx0, _mm512_int2mask(0x00F0), idx1, idx1, 0x00);

    tmp1 = _mm512_i32logather_pd(idx, base0, sizeof(double));
    tmp2 = _mm512_i32logather_pd(idx, base1, sizeof(double));

    v0 = _mm512_shuffle_f64x2(tmp1, tmp2, 0x44 );
    v1 = _mm512_shuffle_f64x2(tmp1, tmp2, 0xEE );
}

static gmx_inline double
gmx_simd_reduce_incr_4_return_sum_hsimd_d_x86_avx_512f(double *           m,
                                                       gmx_simd_double_t  v0,
                                                       gmx_simd_double_t  v1)
{
    double   d;
    __m512d  t0, t1;

    assert((size_t)m % 32 == 0);

    t0 = _mm512_add_pd(v0, _mm512_permutex_pd(v0, 0x4E));
    t0 = _mm512_mask_add_pd(t0, _mm512_int2mask(0xCC), v1, _mm512_permutex_pd(v1, 0x4E));
    t0 = _mm512_add_pd(t0, _mm512_permutex_pd(t0, 0xB1));
    t0 = _mm512_mask_shuffle_f64x2(t0, _mm512_int2mask(0xAA), t0, t0, 0xEE);

    t1 = _mm512_maskz_loadu_pd(_mm512_int2mask(0xF), m);
    t1 = _mm512_add_pd(t1, t0);
    _mm512_mask_storeu_pd(m, _mm512_int2mask(0xF), t1);

    t0 = _mm512_add_pd(t0, _mm512_permutex_pd(t0, 0x4E));
    t0 = _mm512_add_pd(t0, _mm512_permutex_pd(t0, 0xB1));

    _mm512_mask_storeu_pd(&d, _mm512_int2mask(0x01), t0);
    return d;
}
#endif


#endif /* GMX_SIMD_IMPL_X86_AVX_512F_UTIL_DOUBLE_H */
