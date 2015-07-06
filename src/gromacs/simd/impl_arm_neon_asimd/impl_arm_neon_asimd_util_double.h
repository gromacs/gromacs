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

#ifndef GMX_SIMD_IMPL_ARM_NEON_ASIMD_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_ARM_NEON_ASIMD_UTIL_DOUBLE_H

#include "config.h"

#include <assert.h>
#include <math.h>

#include <arm_neon.h>

#include "impl_arm_neon_asimd_common.h"

/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_d             gmx_simd_gather_load_transpose_d_arm_neon_asimd
#define gmx_simd_best_pair_alignment_d               gmx_simd_best_pair_alignment_d_arm_neon_asimd
#define gmx_simd_gather_loadu_transpose_d            gmx_simd_gather_loadu_transpose_d_arm_neon_asimd
#define gmx_simd_transpose_scatter_storeu_d          gmx_simd_transpose_scatter_storeu_d_arm_neon_asimd
#define gmx_simd_transpose_scatter_incru_d           gmx_simd_transpose_scatter_incru_d_arm_neon_asimd
#define gmx_simd_transpose_scatter_decru_d           gmx_simd_transpose_scatter_decru_d_arm_neon_asimd
#define gmx_simd_expand_scalars_to_triplets_d        gmx_simd_expand_scalars_to_triplets_d_arm_neon_asimd
#define gmx_simd_gather_load_bysimdint_transpose_d   gmx_simd_gather_load_bysimdint_transpose_d_arm_neon_asimd
#define gmx_simd_gather_loadu_bysimdint_transpose_d  gmx_simd_gather_loadu_bysimdint_transpose_d_arm_neon_asimd
#define gmx_simd_reduce_incr_4_return_sum_d          gmx_simd_reduce_incr_4_return_sum_d_arm_neon_asimd

/****************************************************
 * Double precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus
template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_transpose_d_arm_neon_asimd(const double *        base,
                                                const gmx_int32_t     offset[],
                                                gmx_simd_double_t    &v0,
                                                gmx_simd_double_t    &v1,
                                                gmx_simd_double_t    &v2,
                                                gmx_simd_double_t    &v3)
{
    float64x2_t t1, t2, t3, t4;

    assert((size_t)offset % 8 == 0);
    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);

    t1  = gmx_simd_load_d(base + align * offset[0]);
    t2  = gmx_simd_load_d(base + align * offset[1]);
    t3  = gmx_simd_load_d(base + align * offset[0] + 2);
    t4  = gmx_simd_load_d(base + align * offset[1] + 2);
    v0  = vuzp1q_f64(t1, t2);
    v1  = vuzp2q_f64(t1, t2);
    v2  = vuzp1q_f64(t3, t4);
    v3  = vuzp2q_f64(t3, t4);
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_transpose_d_arm_neon_asimd(const double *        base,
                                                const gmx_int32_t     offset[],
                                                gmx_simd_double_t    &v0,
                                                gmx_simd_double_t    &v1)
{
    float64x2_t t1, t2;

    assert((size_t)offset % 8 == 0);
    assert((size_t)base % 16 == 0);
    assert(align % 2 == 0);

    t1  = gmx_simd_load_d(base + align * offset[0]);
    t2  = gmx_simd_load_d(base + align * offset[1]);
    v0  = vuzp1q_f64(t1, t2);
    v1  = vuzp2q_f64(t1, t2);
}

static const int gmx_simd_best_pair_alignment_d_arm_neon_asimd = 2;

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_loadu_transpose_d_arm_neon_asimd(const double *        base,
                                                 const gmx_int32_t     offset[],
                                                 gmx_simd_double_t    &v0,
                                                 gmx_simd_double_t    &v1,
                                                 gmx_simd_double_t    &v2)
{
    float64x2_t t1, t2;
    float64x1_t t3, t4;

    assert((size_t)offset % 8 == 0);

    t1  = gmx_simd_loadu_d(base + align * offset[0]);
    t2  = gmx_simd_loadu_d(base + align * offset[1]);
    t3  = vld1_f64(base + align * offset[0] + 2);
    t4  = vld1_f64(base + align * offset[1] + 2);
    v0  = vuzp1q_f64(t1, t2);
    v1  = vuzp2q_f64(t1, t2);
    v2  = vcombine_f64(t3, t4);
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_storeu_d_arm_neon_asimd(double *              base,
                                                   const gmx_int32_t     offset[],
                                                   gmx_simd_double_t     v0,
                                                   gmx_simd_double_t     v1,
                                                   gmx_simd_double_t     v2)
{
    float64x2_t t0, t1;

    assert((size_t)offset % 8 == 0);

    t0  = vuzp1q_f64(v0, v1);
    t1  = vuzp2q_f64(v0, v1);
    vst1q_f64(base + align * offset[0], t0);
    vst1q_f64(base + align * offset[1], t1);
    vst1_f64(base + align * offset[0] + 2, vget_low_f64(v2));
    vst1_f64(base + align * offset[1] + 2, vget_high_f64(v2));
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_incru_d_arm_neon_asimd(double *              base,
                                                  const gmx_int32_t     offset[],
                                                  gmx_simd_double_t     v0,
                                                  gmx_simd_double_t     v1,
                                                  gmx_simd_double_t     v2)
{
    float64x2_t t0, t1, t2, t3;

    assert((size_t)offset % 8 == 0);

    t0  = vuzp1q_f64(v0, v1);
    t1  = vuzp2q_f64(v0, v1);

    t2  = gmx_simd_loadu_d(base + align * offset[0]);
    t3  = gmx_simd_loadu_d(base + align * offset[1]);
    t2  = gmx_simd_add_d(t2, t0);
    t3  = gmx_simd_add_d(t3, t1);
    gmx_simd_storeu_d(base + align * offset[0], t2);
    gmx_simd_storeu_d(base + align * offset[1], t3);

    base[ align * offset[0] + 2 ] += vget_low_f64(v2);
    base[ align * offset[1] + 2 ] += vget_high_f64(v2);
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_transpose_scatter_decru_d_arm_neon_asimd(double *              base,
                                                  const gmx_int32_t     offset[],
                                                  gmx_simd_double_t     v0,
                                                  gmx_simd_double_t     v1,
                                                  gmx_simd_double_t     v2)
{
    float64x2_t t0, t1, t2, t3;

    assert((size_t)offset % 8 == 0);

    t0  = vuzp1q_f64(v0, v1);
    t1  = vuzp2q_f64(v0, v1);

    t2  = gmx_simd_loadu_d(base + align * offset[0]);
    t3  = gmx_simd_loadu_d(base + align * offset[1]);
    t2  = gmx_simd_sub_d(t2, t0);
    t3  = gmx_simd_sub_d(t3, t1);
    gmx_simd_storeu_d(base + align * offset[0], t2);
    gmx_simd_storeu_d(base + align * offset[1], t3);

    base[ align * offset[0] + 2 ] -= vget_low_f64(v2);
    base[ align * offset[1] + 2 ] -= vget_high_f64(v2);
}

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_d_arm_neon_asimd(gmx_simd_double_t    scalar,
                                                     gmx_simd_double_t   &triplets0,
                                                     gmx_simd_double_t   &triplets1,
                                                     gmx_simd_double_t   &triplets2)
{
    triplets0 = vuzp1q_f64(scalar, scalar);
    triplets1 = scalar;
    triplets2 = vuzp2q_f64(scalar, scalar);
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_d_arm_neon_asimd(const double *       base,
                                                          gmx_simd_dint32_t    offset,
                                                          gmx_simd_double_t   &v0,
                                                          gmx_simd_double_t   &v1,
                                                          gmx_simd_double_t   &v2,
                                                          gmx_simd_double_t   &v3)
{
    int ioffset[4] __attribute__((aligned(16)));

    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);

    gmx_simd_store_di(ioffset, offset);
    gmx_simd_gather_load_transpose_d_arm_neon_asimd<align>(base, ioffset, v0, v1, v2, v3);
}


template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_load_bysimdint_transpose_d_arm_neon_asimd(const double *        base,
                                                          gmx_simd_dint32_t     offset,
                                                          gmx_simd_double_t    &v0,
                                                          gmx_simd_double_t    &v1)
{
    int ioffset[4] __attribute__((aligned(16)));
    gmx_simd_store_di(ioffset, offset);
    gmx_simd_gather_load_transpose_d_arm_neon_asimd<align>(base, ioffset, v0, v1);
}

template <int align>
static gmx_inline void gmx_simdcall
gmx_simd_gather_loadu_bysimdint_transpose_d_arm_neon_asimd(const double *       base,
                                                           gmx_simd_dint32_t    offset,
                                                           gmx_simd_double_t   &v0,
                                                           gmx_simd_double_t   &v1)
{
    float64x2_t t1, t2;
    int         ioffset[4] __attribute__((aligned(16)));

    gmx_simd_store_di(ioffset, offset);
    t1  = gmx_simd_loadu_d(base + align * ioffset[0]);
    t2  = gmx_simd_loadu_d(base + align * ioffset[1]);
    v0  = vuzp1q_f64(t1, t2);
    v1  = vuzp2q_f64(t1, t2);
}


static gmx_inline double gmx_simdcall
gmx_simd_reduce_incr_4_return_sum_d_arm_neon_asimd(double *           m,
                                                   gmx_simd_double_t  v0,
                                                   gmx_simd_double_t  v1,
                                                   gmx_simd_double_t  v2,
                                                   gmx_simd_double_t  v3)
{
    float64x2_t t1, t2, t3, t4;

    assert((size_t)m % 16 == 0);

    t1 = vuzp1q_f64(v0, v1);
    t2 = vuzp2q_f64(v0, v1);
    t3 = vuzp1q_f64(v2, v3);
    t4 = vuzp2q_f64(v2, v3);

    t1 = gmx_simd_add_d(t1, t2);
    t3 = gmx_simd_add_d(t3, t4);

    t2 = gmx_simd_add_d(t1, gmx_simd_load_d(m));
    t4 = gmx_simd_add_d(t3, gmx_simd_load_d(m + 2));
    gmx_simd_store_d(m, t2);
    gmx_simd_store_d(m + 2, t4);

    t1 = gmx_simd_add_d(t1, t3);
    return gmx_simd_reduce_d_arm_neon_asimd(t1);
}
#endif /* __cplusplus */

#endif /* GMX_SIMD_IMPL_ARM_NEON_ASIMD_UTIL_DOUBLE_H */
