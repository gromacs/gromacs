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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_QPX_UTIL_DOUBLE_H
#define GMX_SIMD_IMPLEMENTATION_IBM_QPX_UTIL_DOUBLE_H

#if !defined(__ibmxl__) && !defined(__xlC__)
/* xlc-13 (at least) seems to be buggy with the asserts at high optimization */
#    include <assert.h>
#endif
#include <math.h>
#ifdef __clang__
#    include <qpxmath.h>
#endif

#include "impl_ibm_qpx_common.h"

/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_d             gmx_simd_gather_load_transpose_d_ibm_qpx
#define gmx_simd_best_pair_alignment_d               gmx_simd_best_pair_alignment_d_ibm_qpx
#define gmx_simd_gather_loadu_transpose_d            gmx_simd_gather_loadu_transpose_d_ibm_qpx
#define gmx_simd_transpose_scatter_storeu_d          gmx_simd_transpose_scatter_storeu_d_ibm_qpx
#define gmx_simd_transpose_scatter_incru_d           gmx_simd_transpose_scatter_incru_d_ibm_qpx
#define gmx_simd_transpose_scatter_decru_d           gmx_simd_transpose_scatter_decru_d_ibm_qpx
#define gmx_simd_expand_scalars_to_triplets_d        gmx_simd_expand_scalars_to_triplets_ibm_qpx
#define gmx_simd_gather_load_bysimdint_transpose_d   gmx_simd_gather_load_bysimdint_transpose_d_ibm_qpx
#define gmx_simd_gather_loadu_bysimdint_transpose_d  gmx_simd_gather_loadu_bysimdint_transpose_d_ibm_qpx
#define gmx_simd_reduce_incr_4_return_sum_d          gmx_simd_reduce_incr_4_return_sum_d_ibm_qpx

/****************************************************
 * Double precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

#define GMX_QPX_TRANSPOSE4(v0, v1, v2, v3)                        \
    {                                                                 \
        vector4double gmx_qpx_t0 = vec_perm(v0, v2, vec_gpci(00415)); \
        vector4double gmx_qpx_t1 = vec_perm(v0, v2, vec_gpci(02637)); \
        vector4double gmx_qpx_t2 = vec_perm(v1, v3, vec_gpci(00415)); \
        vector4double gmx_qpx_t3 = vec_perm(v1, v3, vec_gpci(02637)); \
        v0 = vec_perm(gmx_qpx_t0, gmx_qpx_t2, vec_gpci(00415));       \
        v1 = vec_perm(gmx_qpx_t0, gmx_qpx_t2, vec_gpci(02637));       \
        v2 = vec_perm(gmx_qpx_t1, gmx_qpx_t3, vec_gpci(00415));       \
        v3 = vec_perm(gmx_qpx_t1, gmx_qpx_t3, vec_gpci(02637));       \
    }


template <int align>
static __attribute__((always_inline)) void
gmx_simd_gather_load_transpose_d_ibm_qpx(const double *        base,
                                         const gmx_int32_t     offset[],
                                         gmx_simd_double_t    &v0,
                                         gmx_simd_double_t    &v1,
                                         gmx_simd_double_t    &v2,
                                         gmx_simd_double_t    &v3)
{
#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 32 == 0);
    assert(align % 4 == 0);
#endif

    v0  = gmx_simd_load_d( base + align * offset[0] );
    v1  = gmx_simd_load_d( base + align * offset[1] );
    v2  = gmx_simd_load_d( base + align * offset[2] );
    v3  = gmx_simd_load_d( base + align * offset[3] );
    GMX_QPX_TRANSPOSE4(v0, v1, v2, v3);
}

template <int align>
static __attribute__((always_inline)) void
gmx_simd_gather_load_transpose_d_ibm_qpx(const double *       base,
                                         const gmx_int32_t    offset[],
                                         gmx_simd_double_t   &v0,
                                         gmx_simd_double_t   &v1)
{
#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 16 == 0);
    assert(align % 2 == 0);
#endif

    gmx_simd_double_t t0, t1, t2, t3;
    t0 = vec_ld2(0, (double *)(base + align * offset[0]) );
    t1 = vec_ld2(0, (double *)(base + align * offset[1]) );
    t2 = vec_ld2(0, (double *)(base + align * offset[2]) );
    t3 = vec_ld2(0, (double *)(base + align * offset[3]) );
    t0 = vec_perm(t0, t2, vec_gpci(00415));
    t1 = vec_perm(t1, t3, vec_gpci(00415));
    v0 = vec_perm(t0, t1, vec_gpci(00415));
    v1 = vec_perm(t0, t1, vec_gpci(02637));
}

static const int gmx_simd_best_pair_alignment_d_ibm_qpx = 2;

template <int align>
static __attribute__((always_inline)) void
gmx_simd_gather_loadu_transpose_d_ibm_qpx(const double *       base,
                                          const gmx_int32_t    offset[],
                                          gmx_simd_double_t   &v0,
                                          gmx_simd_double_t   &v1,
                                          gmx_simd_double_t   &v2)
{
    gmx_simd_double_t t1, t2, t3, t4, t5, t6, t7, t8;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
#endif

    if ( (align & 0x3) == 0)
    {
        gmx_simd_double_t t3;
        gmx_simd_gather_load_transpose_d_ibm_qpx<align>(base, offset, v0, v1, v2, t3);
    }
    else
    {
        t1  = vec_perm(vec_splats(base[align * offset[0]]), vec_splats(base[align * offset[0] + 1]), vec_gpci(00415));
        t2  = vec_perm(vec_splats(base[align * offset[1]]), vec_splats(base[align * offset[1] + 1]), vec_gpci(00415));
        t3  = vec_perm(vec_splats(base[align * offset[2]]), vec_splats(base[align * offset[2] + 1]), vec_gpci(00415));
        t4  = vec_perm(vec_splats(base[align * offset[3]]), vec_splats(base[align * offset[3] + 1]), vec_gpci(00415));

        t5  = vec_splats( *(base + align * offset[0] + 2) );
        t6  = vec_splats( *(base + align * offset[1] + 2) );
        t7  = vec_splats( *(base + align * offset[2] + 2) );
        t8  = vec_splats( *(base + align * offset[3] + 2) );

        t1  = vec_perm(t1, t2, vec_gpci(00415));
        t3  = vec_perm(t3, t4, vec_gpci(00415));
        v0  = vec_perm(t1, t3, vec_gpci(00145));
        v1  = vec_perm(t3, t1, vec_gpci(06723));
        t5  = vec_perm(t5, t6, vec_gpci(00415));
        t7  = vec_perm(t7, t8, vec_gpci(00415));
        v2  = vec_perm(t5, t7, vec_gpci(00145));
    }
}

template <int align>
static __attribute__((always_inline)) void
gmx_simd_transpose_scatter_storeu_d_ibm_qpx(double *            base,
                                            const gmx_int32_t   offset[],
                                            gmx_simd_double_t   v0,
                                            gmx_simd_double_t   v1,
                                            gmx_simd_double_t   v2)
{
#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
#endif

    base[align * offset[0]    ] = vec_extract(v0, 0);
    base[align * offset[0] + 1] = vec_extract(v1, 0);
    base[align * offset[0] + 2] = vec_extract(v2, 0);
    base[align * offset[1]    ] = vec_extract(v0, 1);
    base[align * offset[1] + 1] = vec_extract(v1, 1);
    base[align * offset[1] + 2] = vec_extract(v2, 1);
    base[align * offset[2]    ] = vec_extract(v0, 2);
    base[align * offset[2] + 1] = vec_extract(v1, 2);
    base[align * offset[2] + 2] = vec_extract(v2, 2);
    base[align * offset[3]    ] = vec_extract(v0, 3);
    base[align * offset[3] + 1] = vec_extract(v1, 3);
    base[align * offset[3] + 2] = vec_extract(v2, 3);
}

template <int align>
static __attribute__((always_inline)) void
gmx_simd_transpose_scatter_incru_d_ibm_qpx(double *            base,
                                           const gmx_int32_t   offset[],
                                           gmx_simd_double_t   v0,
                                           gmx_simd_double_t   v1,
                                           gmx_simd_double_t   v2)
{
#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
#endif

    base[align * offset[0]    ] += vec_extract(v0, 0);
    base[align * offset[0] + 1] += vec_extract(v1, 0);
    base[align * offset[0] + 2] += vec_extract(v2, 0);
    base[align * offset[1]    ] += vec_extract(v0, 1);
    base[align * offset[1] + 1] += vec_extract(v1, 1);
    base[align * offset[1] + 2] += vec_extract(v2, 1);
    base[align * offset[2]    ] += vec_extract(v0, 2);
    base[align * offset[2] + 1] += vec_extract(v1, 2);
    base[align * offset[2] + 2] += vec_extract(v2, 2);
    base[align * offset[3]    ] += vec_extract(v0, 3);
    base[align * offset[3] + 1] += vec_extract(v1, 3);
    base[align * offset[3] + 2] += vec_extract(v2, 3);
}

template <int align>
static __attribute__((always_inline)) void
gmx_simd_transpose_scatter_decru_d_ibm_qpx(double *            base,
                                           const gmx_int32_t   offset[],
                                           gmx_simd_double_t   v0,
                                           gmx_simd_double_t   v1,
                                           gmx_simd_double_t   v2)
{
#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
#endif

    base[align * offset[0]    ] -= vec_extract(v0, 0);
    base[align * offset[0] + 1] -= vec_extract(v1, 0);
    base[align * offset[0] + 2] -= vec_extract(v2, 0);
    base[align * offset[1]    ] -= vec_extract(v0, 1);
    base[align * offset[1] + 1] -= vec_extract(v1, 1);
    base[align * offset[1] + 2] -= vec_extract(v2, 1);
    base[align * offset[2]    ] -= vec_extract(v0, 2);
    base[align * offset[2] + 1] -= vec_extract(v1, 2);
    base[align * offset[2] + 2] -= vec_extract(v2, 2);
    base[align * offset[3]    ] -= vec_extract(v0, 3);
    base[align * offset[3] + 1] -= vec_extract(v1, 3);
    base[align * offset[3] + 2] -= vec_extract(v2, 3);
}

template <int align>
static __attribute__((always_inline)) void
gmx_simd_gather_load_bysimdint_transpose_d_ibm_qpx(const double *      base,
                                                   gmx_simd_dint32_t   offset,
                                                   gmx_simd_double_t  &v0,
                                                   gmx_simd_double_t  &v1,
                                                   gmx_simd_double_t  &v2,
                                                   gmx_simd_double_t  &v3)
{
    int __attribute__ ((aligned (16))) ioffset[4];

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 32 == 0);
    assert(align % 4 == 0);
#endif

    gmx_simd_store_di(ioffset, offset);
    gmx_simd_gather_load_transpose_d_ibm_qpx<align>(base, ioffset, v0, v1, v2, v3);
}


template <int align>
static __attribute__((always_inline)) void
gmx_simd_gather_load_bysimdint_transpose_d_ibm_qpx(const double  *      base,
                                                   gmx_simd_dint32_t    offset,
                                                   gmx_simd_double_t   &v0,
                                                   gmx_simd_double_t   &v1)
{
    int __attribute__ ((aligned (16))) ioffset[4];

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 16 == 0);
    assert(align % 2 == 0);
#endif

    gmx_simd_store_di(ioffset, offset);
    gmx_simd_gather_load_transpose_d_ibm_qpx<align>(base, ioffset, v0, v1);
}

static __attribute__((always_inline)) double
gmx_simd_reduce_incr_4_return_sum_d_ibm_qpx(double *           m,
                                            gmx_simd_double_t  v0,
                                            gmx_simd_double_t  v1,
                                            gmx_simd_double_t  v2,
                                            gmx_simd_double_t  v3)
{
#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)m % 32 == 0);
#endif

    GMX_QPX_TRANSPOSE4(v0, v1, v2, v3);
    v0 = vec_add(v0, v1);
    v2 = vec_add(v2, v3);
    v0 = vec_add(v0, v2);
    v2 = vec_add(v0, gmx_simd_load_d(m));
    gmx_simd_store_d(m, v2);

    return gmx_simd_reduce_ibm_qpx(v0);
}
#endif /* __cplusplus */


#endif /* GMX_SIMD_IMPLEMENTATION_IBM_QPX_UTIL_DOUBLE_H */
