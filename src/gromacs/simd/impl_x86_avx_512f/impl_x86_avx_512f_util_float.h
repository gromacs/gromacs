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

#ifndef GMX_SIMD_IMPL_X86_AVX_512F_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_512F_UTIL_FLOAT_H

#include "config.h"

#include <assert.h>
#include <math.h>

#include <immintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_avx_512f_common.h"

/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_f             gmx_simd_gather_load_transpose_f_x86_avx_512f
#define gmx_simd_best_pair_alignment_f               gmx_simd_best_pair_alignment_f_x86_avx_512f
#define gmx_simd_gather_loadu_transpose_f            gmx_simd_gather_loadu_transpose_f_x86_avx_512f
#define gmx_simd_transpose_scatter_storeu_f          gmx_simd_transpose_scatter_storeu_f_x86_avx_512f
#define gmx_simd_transpose_scatter_incru_f           gmx_simd_transpose_scatter_incru_f_x86_avx_512f
#define gmx_simd_transpose_scatter_decru_f           gmx_simd_transpose_scatter_decru_f_x86_avx_512f
#define gmx_simd_expand_scalars_to_triplets_f        gmx_simd_expand_scalars_to_triplets_f_x86_avx_512f
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_x86_avx_512f
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_x86_avx_512f
/* Half-simd-width utilities */
#define gmx_simd_load_dual_hsimd_f                   gmx_simd_load_dual_hsimd_f_x86_avx_512f
#define gmx_simd_loaddup_hsimd_f                     gmx_simd_loaddup_hsimd_f_x86_avx_512f
#define gmx_simd_load1_dual_hsimd_f                  gmx_simd_load1_dual_hsimd_f_x86_avx_512f
#define gmx_simd_store_dual_hsimd_f                  gmx_simd_store_dual_hsimd_f_x86_avx_512f
#define gmx_simd_decr_hsimd_f                        gmx_simd_decr_hsimd_f_x86_avx_512f
#define gmx_simd_gather_load_transpose_hsimd_f       gmx_simd_gather_load_transpose_hsimd_f_x86_avx_512f
#define gmx_simd_reduce_incr_4_return_sum_hsimd_f    gmx_simd_reduce_incr_4_return_sum_hsimd_f_x86_avx_512f

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

/* Declare functions implemented further down that are used here */
template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f(const float *       base,
                                                        gmx_simd_fint32_t   simdoffset,
                                                        gmx_simd_float_t   &v0,
                                                        gmx_simd_float_t   &v1,
                                                        gmx_simd_float_t   &v2,
                                                        gmx_simd_float_t   &v3);

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f(const float *       base,
                                                        gmx_simd_fint32_t   simdoffset,
                                                        gmx_simd_float_t   &v0,
                                                        gmx_simd_float_t   &v1);


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_x86_avx_512f(const float *        base,
                                              const gmx_int32_t    offset[],
                                              gmx_simd_float_t    &v0,
                                              gmx_simd_float_t    &v1,
                                              gmx_simd_float_t    &v2,
                                              gmx_simd_float_t    &v3)
{
    assert((size_t)offset % 64 == 0);
    gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f<align>(base,
                                                                   gmx_simd_load_fi(offset),
                                                                   v0, v1, v2, v3);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_x86_avx_512f(const float *        base,
                                              const gmx_int32_t    offset[],
                                              gmx_simd_float_t    &v0,
                                              gmx_simd_float_t    &v1)
{
    assert((size_t)offset % 64 == 0);
    gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f<align>(base,
                                                                   gmx_simd_load_fi(offset),
                                                                   v0, v1);
}


static const int gmx_simd_best_pair_alignment_f_x86_avx_512f = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_f_x86_avx_512f(const float *        base,
                                               const gmx_int32_t    offset[],
                                               gmx_simd_float_t    &v0,
                                               gmx_simd_float_t    &v1,
                                               gmx_simd_float_t    &v2)
{
    __m512i simdoffset;

    assert((size_t)offset % 64 == 0);

    simdoffset = gmx_simd_load_fi(offset);

    if (align == 4)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm512_mullo_epi32(simdoffset, _mm512_set1_epi32(align));
    }

    v0  = _mm512_i32gather_ps(simdoffset, base,   sizeof(float));
    v1  = _mm512_i32gather_ps(simdoffset, base+1, sizeof(float));
    v2  = _mm512_i32gather_ps(simdoffset, base+2, sizeof(float));
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_f_x86_avx_512f(float *              base,
                                                 const gmx_int32_t    offset[],
                                                 gmx_simd_float_t     v0,
                                                 gmx_simd_float_t     v1,
                                                 gmx_simd_float_t     v2)
{
    __m512i simdoffset;

    assert((size_t)offset % 64 == 0);

    simdoffset = gmx_simd_load_fi(offset);

    if (align == 4)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm512_mullo_epi32(simdoffset, _mm512_set1_epi32(align));
    }

    _mm512_i32scatter_ps(base,   simdoffset, v0, sizeof(float));
    _mm512_i32scatter_ps(base+1, simdoffset, v1, sizeof(float));
    _mm512_i32scatter_ps(base+2, simdoffset, v2, sizeof(float));
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_f_x86_avx_512f(float *              base,
                                                const gmx_int32_t    offset[],
                                                gmx_simd_float_t     v0,
                                                gmx_simd_float_t     v1,
                                                gmx_simd_float_t     v2)
{
    __m512 t0, t1, t2;
    gmx_simd_gather_loadu_transpose_f_x86_avx_512f<align>(base, offset, t0, t1, t2);
    t0 = _mm512_add_ps(t0, v0);
    t1 = _mm512_add_ps(t1, v1);
    t2 = _mm512_add_ps(t2, v2);
    gmx_simd_transpose_scatter_storeu_f_x86_avx_512f<align>(base, offset, t0, t1, t2);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_f_x86_avx_512f(float *              base,
                                                const gmx_int32_t    offset[],
                                                gmx_simd_float_t     v0,
                                                gmx_simd_float_t     v1,
                                                gmx_simd_float_t     v2)
{
    __m512 t0, t1, t2;
    gmx_simd_gather_loadu_transpose_f_x86_avx_512f<align>(base, offset, t0, t1, t2);
    t0 = _mm512_sub_ps(t0, v0);
    t1 = _mm512_sub_ps(t1, v1);
    t2 = _mm512_sub_ps(t2, v2);
    gmx_simd_transpose_scatter_storeu_f_x86_avx_512f<align>(base, offset, t0, t1, t2);
}

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_f_x86_avx_512f(gmx_simd_float_t   scalar,
                                                   gmx_simd_float_t  &triplets0,
                                                   gmx_simd_float_t  &triplets1,
                                                   gmx_simd_float_t  &triplets2)
{
    triplets0 = _mm512_castsi512_ps(_mm512_permutevar_epi32(_mm512_set_epi32(5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, 0),
                                                            _mm512_castps_si512(scalar)));
    triplets1 = _mm512_castsi512_ps(_mm512_permutevar_epi32(_mm512_set_epi32(10, 10, 9, 9, 9, 8, 8, 8, 7, 7, 7, 6, 6, 6, 5, 5),
                                                            _mm512_castps_si512(scalar)));
    triplets2 = _mm512_castsi512_ps(_mm512_permutevar_epi32(_mm512_set_epi32(15, 15, 15, 14, 14, 14, 13, 13, 13, 12, 12, 12, 11, 11, 11, 10),
                                                            _mm512_castps_si512(scalar)));
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f(const float *       base,
                                                        gmx_simd_fint32_t   simdoffset,
                                                        gmx_simd_float_t   &v0,
                                                        gmx_simd_float_t   &v1,
                                                        gmx_simd_float_t   &v2,
                                                        gmx_simd_float_t   &v3)
{
    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);

    if (align == 4)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 2);
    }
    else if (align == 8)
    {
        simdoffset = _mm512_slli_epi32(simdoffset, 3);
    }
    else
    {
        simdoffset = _mm512_mullo_epi32(simdoffset, _mm512_set1_epi32(align));
    }

    v0  = _mm512_i32gather_ps(simdoffset, base,   sizeof(float));
    v1  = _mm512_i32gather_ps(simdoffset, base+1, sizeof(float));
    v2  = _mm512_i32gather_ps(simdoffset, base+2, sizeof(float));
    v3  = _mm512_i32gather_ps(simdoffset, base+3, sizeof(float));
}

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_f_x86_avx_512f(const float *       base,
                                                         gmx_simd_fint32_t   simdoffset,
                                                         gmx_simd_float_t   &v0,
                                                         gmx_simd_float_t   &v1)
{
    if (align == 2)
    {
        v0  = _mm512_i32gather_ps(simdoffset, base,   align * sizeof(float));
        v1  = _mm512_i32gather_ps(simdoffset, base+1, align * sizeof(float));
    }
    else
    {
        if (align == 4)
        {
            simdoffset = _mm512_slli_epi32(simdoffset, 2);
        }
        else if (align == 8)
        {
            simdoffset = _mm512_slli_epi32(simdoffset, 3);
        }
        else
        {
            simdoffset = _mm512_mullo_epi32(simdoffset, _mm512_set1_epi32(align));
        }
        v0  = _mm512_i32gather_ps(simdoffset, base,   sizeof(float));
        v1  = _mm512_i32gather_ps(simdoffset, base+1, sizeof(float));
    }
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_x86_avx_512f(const float *       base,
                                                        gmx_simd_fint32_t   simdoffset,
                                                        gmx_simd_float_t   &v0,
                                                        gmx_simd_float_t   &v1)
{
    assert((size_t)base % 8 == 0);
    assert(align % 2 == 0);
    gmx_simd_gather_loadu_bysimdint_transpose_f_x86_avx_512f<align>(base, simdoffset, v0, v1);
}

static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_f_x86_avx_512f(float *           m,
                                                 gmx_simd_float_t  v0,
                                                 gmx_simd_float_t  v1,
                                                 gmx_simd_float_t  v2,
                                                 gmx_simd_float_t  v3)
{
    float  f;
    __m512 t0, t1, t2, t3;

    assert((size_t)m % 16 == 0);

    t0 = _mm512_add_ps(v0, _mm512_permute_ps(v0, 0x4E));
    t0 = _mm512_mask_add_ps(t0, _mm512_int2mask(0xCCCC), v2, _mm512_permute_ps(v2, 0x4E));
    t1 = _mm512_add_ps(v1, _mm512_permute_ps(v1, 0x4E));
    t1 = _mm512_mask_add_ps(t1, _mm512_int2mask(0xCCCC), v3, _mm512_permute_ps(v3, 0x4E));
    t2 = _mm512_add_ps(t0, _mm512_permute_ps(t0, 0xB1));
    t2 = _mm512_mask_add_ps(t2, _mm512_int2mask(0xAAAA), t1, _mm512_permute_ps(t1, 0xB1));

    t2 = _mm512_add_ps(t2, _mm512_shuffle_f32x4(t2, t2, 0x4E));
    t2 = _mm512_add_ps(t2, _mm512_shuffle_f32x4(t2, t2, 0xB1));

    t0 = _mm512_maskz_loadu_ps(_mm512_int2mask(0xF), m);
    t0 = _mm512_add_ps(t0, t2);
    _mm512_mask_storeu_ps(m, _mm512_int2mask(0xF), t0);

    t2 = _mm512_add_ps(t2, _mm512_permute_ps(t2, 0x4E));
    t2 = _mm512_add_ps(t2, _mm512_permute_ps(t2, 0xB1));

    _mm512_mask_storeu_ps(&f, _mm512_int2mask(0x0001), t2);
    return f;
}

/******************************************************
 * Single precision half-simd-width utility functions *
 ******************************************************/
static gmx_inline gmx_simd_float_t
gmx_simd_load_dual_hsimd_f_x86_avx_512f(const float * m0,
                                        const float * m1)
{
    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    return _mm512_mask_loadu_ps(_mm512_maskz_loadu_ps(_mm512_int2mask(0x00FF), m0), _mm512_int2mask(0xFF00), m1-8);
}

static gmx_inline gmx_simd_float_t
gmx_simd_loaddup_hsimd_f_x86_avx_512f(const float * m)
{
    assert((size_t)m % 32 == 0);

    __m512 tmp;

    tmp = _mm512_maskz_loadu_ps(_mm512_int2mask(0x00FF), m);
    return _mm512_shuffle_f32x4(tmp, tmp, 0x44);
}

static gmx_inline gmx_simd_float_t
gmx_simd_load1_dual_hsimd_f_x86_avx_512f(const float * m)
{
    __m512 tmp;
    tmp = _mm512_maskz_expandloadu_ps(_mm512_int2mask(0x0101), m);
    tmp = _mm512_permute_ps(tmp, 0x00);
    return _mm512_shuffle_f32x4(tmp, tmp, 0xA0);
}


static gmx_inline void
gmx_simd_store_dual_hsimd_f_x86_avx_512f(float *           m0,
                                         float *           m1,
                                         gmx_simd_float_t  a)
{
    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    _mm512_mask_storeu_ps(m0,   _mm512_int2mask(0x00FF), a);
    _mm512_mask_storeu_ps(m1-8, _mm512_int2mask(0xFF00), a);
}

static gmx_inline void
gmx_simd_decr_hsimd_f_x86_avx_512f(float *          m,
                                   gmx_simd_float_t a)
{
    __m512 t;

    assert((size_t)m % 32 == 0);

    a = _mm512_add_ps(a, _mm512_shuffle_f32x4(a, a, 0xEE));
    t = _mm512_maskz_loadu_ps(_mm512_int2mask(0x00FF), m);
    t = _mm512_sub_ps(t, a);
    _mm512_mask_storeu_ps(m, _mm512_int2mask(0x00FF), t);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_hsimd_f_x86_avx_512f(const float *       base0,
                                                    const float *       base1,
                                                    gmx_int32_t         offset[],
                                                    gmx_simd_float_t   &v0,
                                                    gmx_simd_float_t   &v1)
{
    __m512i idx0, idx1, idx;
    __m512  tmp1, tmp2;

    assert((size_t)offset % 32 == 0);
    assert((size_t)base0 % 8 == 0);
    assert((size_t)base1 % 8 == 0);
    assert((size_t)align % 2 == 0);

    idx0 = _mm512_maskz_loadu_epi32(_mm512_int2mask(0x00FF), offset);

    idx0 = _mm512_mullo_epi32(idx0, _mm512_set1_epi32(align));
    idx1 = _mm512_add_epi32(idx0, _mm512_set1_epi32(1));

    idx = _mm512_mask_shuffle_i32x4(idx0, _mm512_int2mask(0xFF00), idx1, idx1, 0x44);

    tmp1 = _mm512_i32gather_ps(idx, base0, sizeof(float));
    tmp2 = _mm512_i32gather_ps(idx, base1, sizeof(float));

    v0 = _mm512_shuffle_f32x4(tmp1, tmp2, 0x44 );
    v1 = _mm512_shuffle_f32x4(tmp1, tmp2, 0xEE );
}

static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_hsimd_f_x86_avx_512f(float *            m,
                                                       gmx_simd_float_t   v0,
                                                       gmx_simd_float_t   v1)
{
    float  f;
    __m512 t0, t1;

    assert((size_t)m % 32 == 0);

    /* This is not optimal, but no point optimizing until we know AVX-512 latencies */
    t0 = _mm512_add_ps(v0, _mm512_permute_ps(v0, 0x4E));
    t1 = _mm512_add_ps(v1, _mm512_permute_ps(v1, 0x4E));
    t0 = _mm512_add_ps(t0, _mm512_permute_ps(t0, 0xB1));
    t0 = _mm512_mask_add_ps(t0, _mm512_int2mask(0xCCCC), t1, _mm512_permute_ps(t1, 0xB1));
    t0 = _mm512_add_ps(t0, _mm512_shuffle_f32x4(t0, t0, 0xB1));
    t0 = _mm512_mask_shuffle_f32x4(t0, _mm512_int2mask(0xAAAA), t0, t0, 0xEE);

    t1 = _mm512_maskz_loadu_ps(_mm512_int2mask(0xF), m);
    t1 = _mm512_add_ps(t1, t0);
    _mm512_mask_storeu_ps(m, _mm512_int2mask(0xF), t1);

    t0 = _mm512_add_ps(t0, _mm512_permute_ps(t0, 0x4E));
    t0 = _mm512_add_ps(t0, _mm512_permute_ps(t0, 0xB1));

    _mm512_mask_storeu_ps(&f, _mm512_int2mask(0x0001), t0);
    return f;
}
#endif

#endif /* GMX_SIMD_IMPL_X86_AVX_512F_UTIL_FLOAT_H */
