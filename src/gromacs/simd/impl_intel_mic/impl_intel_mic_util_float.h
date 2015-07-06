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

#ifndef GMX_SIMD_IMPL_INTEL_MIC_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_INTEL_MIC_UTIL_FLOAT_H

#include "config.h"

#include <math.h>

#include <immintrin.h>

#include "impl_intel_mic_common.h"

/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_f             gmx_simd_gather_load_transpose_f_mic
#define gmx_simd_best_pair_alignment_f               gmx_simd_best_pair_alignment_f_mic
#define gmx_simd_gather_loadu_transpose_f            gmx_simd_gather_loadu_transpose_f_mic
#define gmx_simd_transpose_scatter_storeu_f          gmx_simd_transpose_scatter_storeu_f_mic
#define gmx_simd_transpose_scatter_incru_f           gmx_simd_transpose_scatter_incru_f_mic
#define gmx_simd_transpose_scatter_decru_f           gmx_simd_transpose_scatter_decru_f_mic
#define gmx_simd_expand_scalars_to_triplets_f        gmx_simd_expand_scalars_to_triplets_f_mic
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_mic
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_mic
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_mic
/* Half-simd-width utilities */
#define gmx_simd_load_dual_hsimd_f                   gmx_simd_load_dual_hsimd_f_mic
#define gmx_simd_loaddup_hsimd_f                     gmx_simd_loaddup_hsimd_f_mic
#define gmx_simd_load1_dual_hsimd_f                  gmx_simd_load1_dual_hsimd_f_mic
#define gmx_simd_store_dual_hsimd_f                  gmx_simd_store_dual_hsimd_f_mic
#define gmx_simd_decr_hsimd_f                        gmx_simd_decr_hsimd_f_mic
#define gmx_simd_gather_load_transpose_hsimd_f       gmx_simd_gather_load_transpose_hsimd_f_mic
#define gmx_simd_reduce_incr_4_return_sum_hsimd_f    gmx_simd_reduce_incr_4_return_sum_hsimd_f_mic

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

/* Declare functions implemented further down that are used here */
template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_f_mic(const float *       base,
                                               gmx_simd_fint32_t   simdoffset,
                                               gmx_simd_float_t   &v0,
                                               gmx_simd_float_t   &v1,
                                               gmx_simd_float_t   &v2,
                                               gmx_simd_float_t   &v3);

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_f_mic(const float *       base,
                                               gmx_simd_fint32_t   simdoffset,
                                               gmx_simd_float_t   &v0,
                                               gmx_simd_float_t   &v1);


template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_transpose_f_mic(const float *        base,
                                     const gmx_int32_t    offset[],
                                     gmx_simd_float_t    &v0,
                                     gmx_simd_float_t    &v1,
                                     gmx_simd_float_t    &v2,
                                     gmx_simd_float_t    &v3)
{
    assert((size_t)offset % 64 == 0);
    gmx_simd_gather_load_bysimdint_transpose_f_mic<align>(base,
                                                          gmx_simd_load_fi(offset),
                                                          v0, v1, v2, v3);
}


template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_transpose_f_mic(const float *        base,
                                     const gmx_int32_t    offset[],
                                     gmx_simd_float_t    &v0,
                                     gmx_simd_float_t    &v1)
{
    assert((size_t)offset % 64 == 0);
    gmx_simd_gather_load_bysimdint_transpose_f_mic<align>(base,
                                                          gmx_simd_load_fi(offset),
                                                          v0, v1);
}


static const int gmx_simd_best_pair_alignment_f_mic = 2;

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_loadu_transpose_f_mic(const float *        base,
                                      const gmx_int32_t    offset[],
                                      gmx_simd_float_t    &v0,
                                      gmx_simd_float_t    &v1,
                                      gmx_simd_float_t    &v2)
{
    __m512i simdoffset;

    assert((size_t)offset % 64 == 0);

    simdoffset = gmx_simd_load_fi(offset);

    /* All instructions might be latency ~4 on MIC, so we use shifts where we
     * only need a single instruction (since the shift parameter is an immediate),
     * but multiplication otherwise.
     */
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
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_storeu_f_mic(float *              base,
                                        const gmx_int32_t    offset[],
                                        gmx_simd_float_t     v0,
                                        gmx_simd_float_t     v1,
                                        gmx_simd_float_t     v2)
{
    __m512i simdoffset;

    assert((size_t)offset % 64 == 0);

    simdoffset = gmx_simd_load_fi(offset);

    /* All instructions might be latency ~4 on MIC, so we use shifts where we
     * only need a single instruction (since the shift parameter is an immediate),
     * but multiplication otherwise.
     */
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
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_incru_f_mic(float *              base,
                                       const gmx_int32_t    offset[],
                                       gmx_simd_float_t     v0,
                                       gmx_simd_float_t     v1,
                                       gmx_simd_float_t     v2)
{
    __m512 t0, t1, t2;
    gmx_simd_gather_loadu_transpose_f_mic<align>(base, offset, t0, t1, t2);
    t0 = _mm512_add_ps(t0, v0);
    t1 = _mm512_add_ps(t1, v1);
    t2 = _mm512_add_ps(t2, v2);
    gmx_simd_transpose_scatter_storeu_f_mic<align>(base, offset, t0, t1, t2);
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_decru_f_mic(float *              base,
                                       const gmx_int32_t    offset[],
                                       gmx_simd_float_t     v0,
                                       gmx_simd_float_t     v1,
                                       gmx_simd_float_t     v2)
{
    __m512 t0, t1, t2;
    gmx_simd_gather_loadu_transpose_f_mic<align>(base, offset, t0, t1, t2);
    t0 = _mm512_sub_ps(t0, v0);
    t1 = _mm512_sub_ps(t1, v1);
    t2 = _mm512_sub_ps(t2, v2);
    gmx_simd_transpose_scatter_storeu_f_mic<align>(base, offset, t0, t1, t2);
}

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_f_mic(gmx_simd_float_t   scalar,
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
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_f_mic(const float *       base,
                                               gmx_simd_fint32_t   simdoffset,
                                               gmx_simd_float_t   &v0,
                                               gmx_simd_float_t   &v1,
                                               gmx_simd_float_t   &v2,
                                               gmx_simd_float_t   &v3)
{
    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);

    /* All instructions might be latency ~4 on MIC, so we use shifts where we
     * only need a single instruction (since the shift parameter is an immediate),
     * but multiplication otherwise.
     */
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
static gmx_inline void gmx_simdcall
gmx_simd_gather_loadu_bysimdint_transpose_f_mic(const float *       base,
                                                gmx_simd_fint32_t   simdoffset,
                                                gmx_simd_float_t   &v0,
                                                gmx_simd_float_t   &v1)
{
    /* All instructions might be latency ~4 on MIC, so we use shifts where we
     * only need a single instruction (since the shift parameter is an immediate),
     * but multiplication otherwise.
     * For align == 2 we can merge the constant into the scale parameter,
     * which can take constants up to 8 in total.
     */
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
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_f_mic(const float *       base,
                                               gmx_simd_fint32_t   simdoffset,
                                               gmx_simd_float_t   &v0,
                                               gmx_simd_float_t   &v1)
{
    assert((size_t)base % 8 == 0);
    assert(align % 2 == 0);
    gmx_simd_gather_loadu_bysimdint_transpose_f_mic<align>(base, simdoffset, v0, v1);
}

static gmx_inline float gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_f_mic(float *           m,
                                        gmx_simd_float_t  v0,
                                        gmx_simd_float_t  v1,
                                        gmx_simd_float_t  v2,
                                        gmx_simd_float_t  v3)
{
    float  f;
    __m512 t0, t1, t2, t3;

    assert((size_t)m % 16 == 0);

    t0 = _mm512_add_ps(v0, _mm512_swizzle_ps(v0, _MM_SWIZ_REG_BADC));
    t0 = _mm512_mask_add_ps(t0, _mm512_int2mask(0xCCCC), v2, _mm512_swizzle_ps(v2, _MM_SWIZ_REG_BADC));
    t1 = _mm512_add_ps(v1, _mm512_swizzle_ps(v1, _MM_SWIZ_REG_BADC));
    t1 = _mm512_mask_add_ps(t1, _mm512_int2mask(0xCCCC), v3, _mm512_swizzle_ps(v3, _MM_SWIZ_REG_BADC));
    t2 = _mm512_add_ps(t0, _mm512_swizzle_ps(t0, _MM_SWIZ_REG_CDAB));
    t2 = _mm512_mask_add_ps(t2, _mm512_int2mask(0xAAAA), t1, _mm512_swizzle_ps(t1, _MM_SWIZ_REG_CDAB));

    t2 = _mm512_add_ps(t2, _mm512_permute4f128_ps(t2, _MM_PERM_BADC));
    t2 = _mm512_add_ps(t2, _mm512_permute4f128_ps(t2, _MM_PERM_CDAB));

    t0 = gmx_simd4_load_f(m);
    t0 = gmx_simd4_add_f(t0, t2);
    gmx_simd4_store_f(m, t0);

    t2 = _mm512_add_ps(t2, _mm512_swizzle_ps(t2, _MM_SWIZ_REG_BADC));
    t2 = _mm512_add_ps(t2, _mm512_swizzle_ps(t2, _MM_SWIZ_REG_CDAB));

    _mm512_mask_packstorelo_ps(&f, _mm512_mask2int(0x1), t2);
    return f;
}

/******************************************************
 * Single precision half-simd-width utility functions *
 ******************************************************/
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_load_dual_hsimd_f_mic(const float * m0,
                               const float * m1)
{
    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    return _mm512_castpd_ps(_mm512_mask_extload_pd(_mm512_extload_pd((const double*)m0, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE),
                                                   _mm512_int2mask(0xF0), (const double*)m1, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE));
}

static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_loaddup_hsimd_f_mic(const float * m)
{
    assert((size_t)m % 32 == 0);

    return _mm512_castpd_ps(_mm512_extload_pd((const double*)m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE));
}

static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_load1_dual_hsimd_f_mic(const float * m)
{
    return _mm512_mask_extload_ps(_mm512_extload_ps(m, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE), _mm512_int2mask(0xFF00),
                                  m+1, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
}


static gmx_inline void gmx_simdcall
gmx_simd_store_dual_hsimd_f_mic(float *           m0,
                                float *           m1,
                                gmx_simd_float_t  a)
{
    __m512 t0;

    assert((size_t)m0 % 32 == 0);
    assert((size_t)m1 % 32 == 0);

    _mm512_mask_packstorelo_ps(m0, _mm512_int2mask(0x00FF), a);
    _mm512_mask_packstorelo_ps(m1, _mm512_int2mask(0xFF00), a);
}

static gmx_inline void gmx_simdcall
gmx_simd_decr_hsimd_f_mic(float *          m,
                          gmx_simd_float_t a)
{
    __m512 t;

    assert((size_t)m % 32 == 0);

    t = _mm512_castpd_ps(_mm512_extload_pd((const double*)m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE));
    a = _mm512_add_ps(a, _mm512_permute4f128_ps(a, _MM_PERM_BADC));
    t = _mm512_sub_ps(t, a);
    _mm512_mask_packstorelo_ps(m, _mm512_int2mask(0x00FF), t);
}


template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_transpose_hsimd_f_mic(const float *       base0,
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

    idx0 = _mm512_loadunpacklo_epi32(_mm512_undefined_epi32(), offset);

    idx0 = _mm512_mullo_epi32(idx0, _mm512_set1_epi32(align));
    idx1 = _mm512_add_epi32(idx0, _mm512_set1_epi32(1));

    idx = _mm512_mask_permute4f128_epi32(idx0, _mm512_int2mask(0xFF00), idx1, PERM_LOW2HIGH);

    tmp1 = _mm512_i32gather_ps(idx, base0, sizeof(float));
    tmp2 = _mm512_i32gather_ps(idx, base1, sizeof(float));

    v0 = _mm512_mask_permute4f128_ps(tmp1, _mm512_int2mask(0xFF00), tmp2, PERM_LOW2HIGH);
    v1 = _mm512_mask_permute4f128_ps(tmp2, _mm512_int2mask(0x00FF), tmp1, PERM_HIGH2LOW);
}

static gmx_inline float gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_hsimd_f_mic(float *            m,
                                              gmx_simd_float_t   v0,
                                              gmx_simd_float_t   v1)
{
    float  f;
    __m512 t0, t1;

    assert((size_t)m % 32 == 0);

    t0 = _mm512_add_ps(v0, _mm512_swizzle_ps(v0, _MM_SWIZ_REG_BADC));
    t0 = _mm512_mask_add_ps(t0, _mm512_int2mask(0xCCCC), v1, _mm512_swizzle_ps(v1, _MM_SWIZ_REG_BADC));
    t0 = _mm512_add_ps(t0, _mm512_swizzle_ps(t0, _MM_SWIZ_REG_CDAB));
    t0 = _mm512_add_ps(t0, _mm512_castpd_ps(_mm512_swizzle_pd(_mm512_castps_pd(t0), _MM_SWIZ_REG_BADC)));
    t0 = _mm512_mask_permute4f128_ps(t0, _mm512_int2mask(0b1010101010101010), t0, _MM_PERM_BADC);
    t1 = gmx_simd4_load_f(m);
    t1 = _mm512_add_ps(t1, t0);
    gmx_simd4_store_f(m, t1);

    t0 = _mm512_add_ps(t0, _mm512_swizzle_ps(t0, _MM_SWIZ_REG_BADC));
    t0 = _mm512_add_ps(t0, _mm512_swizzle_ps(t0, _MM_SWIZ_REG_CDAB));

    _mm512_mask_packstorelo_ps(&f, _mm512_mask2int(0x1), t0);
    return f;
}
#endif


#endif /* GMX_SIMD_IMPL_INTEL_MIC_UTIL_FLOAT_H */
