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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VSX_UTIL_DOUBLE_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VSX_UTIL_DOUBLE_H

#include "config.h"

#if !defined(__ibmxl__) && !defined(__xlC__)
/* xlc-13 (at least) seems to be buggy with the asserts at high optimization */
#    include <assert.h>
#endif

#include <math.h>

#include <altivec.h>

#include "impl_ibm_vsx_common.h"

/* Make sure we do not screw up c++ - undefine vector/bool, and rely on __vector,
 * which is present both on gcc and xlc.
 */
#undef vector

/* g++ is also unhappy with the clash of vector bool and the C++ reserved 'bool',
 * which is solved by undefining bool and reyling on __bool. However, that does
 * not work with xlc, which requires us to use bool. Solve the conflict by
 * defining a new gmx_vsx_bool.
 */
#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__)
#    define gmx_vsx_bool __bool
#    undef  bool
#else
#    define gmx_vsx_bool bool
#endif

/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_d             gmx_simd_gather_load_transpose_d_ibm_vsx
#define gmx_simd_best_pair_alignment_d               gmx_simd_best_pair_alignment_d_ibm_vsx
#define gmx_simd_gather_loadu_transpose_d            gmx_simd_gather_loadu_transpose_d_ibm_vsx
#define gmx_simd_transpose_scatter_storeu_d          gmx_simd_transpose_scatter_storeu_d_ibm_vsx
#define gmx_simd_transpose_scatter_incru_d           gmx_simd_transpose_scatter_incru_d_ibm_vsx
#define gmx_simd_transpose_scatter_decru_d           gmx_simd_transpose_scatter_decru_d_ibm_vsx
#define gmx_simd_expand_scalars_to_triplets_d        gmx_simd_expand_scalars_to_triplets_d_ibm_vsx
#define gmx_simd_gather_load_bysimdint_transpose_d   gmx_simd_gather_load_bysimdint_transpose_d_ibm_vsx
#define gmx_simd_gather_loadu_bysimdint_transpose_d  gmx_simd_gather_loadu_bysimdint_transpose_d_ibm_vsx
#define gmx_simd_reduce_incr_4_return_sum_d          gmx_simd_reduce_incr_4_return_sum_d_ibm_vsx

/****************************************************
 * Double precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus
template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d_ibm_vsx(const double *        base,
                                         const gmx_int32_t     offset[],
                                         gmx_simd_double_t    &v0,
                                         gmx_simd_double_t    &v1,
                                         gmx_simd_double_t    &v2,
                                         gmx_simd_double_t    &v3)
{
    gmx_simd_double_t t1, t2, t3, t4;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 8 == 0);
    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);
#endif

    t1  = *(__vector double *)(base + align * offset[0]);
    t2  = *(__vector double *)(base + align * offset[1]);
    t3  = *(__vector double *)(base + align * offset[0] + 2);
    t4  = *(__vector double *)(base + align * offset[1] + 2);
    v0  = vec_mergeh(t1, t2);
    v1  = vec_mergel(t1, t2);
    v2  = vec_mergeh(t3, t4);
    v3  = vec_mergel(t3, t4);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d_ibm_vsx(const double *        base,
                                         const gmx_int32_t     offset[],
                                         gmx_simd_double_t    &v0,
                                         gmx_simd_double_t    &v1)
{
    gmx_simd_double_t t1, t2;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 8 == 0);
    assert((size_t)base % 16 == 0);
    assert(align % 2 == 0);
#endif

    t1  = *(__vector double *)(base + align * offset[0]);
    t2  = *(__vector double *)(base + align * offset[1]);
    v0  = vec_mergeh(t1, t2);
    v1  = vec_mergel(t1, t2);
}

static const int gmx_simd_best_pair_alignment_d_ibm_vsx = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_d_ibm_vsx(const double *        base,
                                          const gmx_int32_t     offset[],
                                          gmx_simd_double_t    &v0,
                                          gmx_simd_double_t    &v1,
                                          gmx_simd_double_t    &v2)
{
    gmx_simd_double_t t1, t2, t3, t4;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 8 == 0);
#endif

    t1  = gmx_simd_loadu_d(base + align * offset[0]);
    t2  = gmx_simd_loadu_d(base + align * offset[1]);
    t3  = vec_splats(*(base + align * offset[0] + 2));
    t4  = vec_splats(*(base + align * offset[1] + 2));
    v0  = vec_mergeh(t1, t2);
    v1  = vec_mergel(t1, t2);
    v2  = vec_mergeh(t3, t4);
}

/* gcc-4.9 does not detect that vec_extract() uses its argument */
template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_d_ibm_vsx(double *                      base,
                                            const gmx_int32_t             offset[],
                                            gmx_simd_double_t             v0,
                                            gmx_simd_double_t             v1,
                                            gmx_simd_double_t gmx_unused  v2)
{
    gmx_simd_double_t t1, t2;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 8 == 0);
#endif

    t1  = vec_mergeh(v0, v1);
    t2  = vec_mergel(v0, v1);

    gmx_simd_storeu_d(base + align * offset[0], t1);
    base[align * offset[0] + 2]  = vec_extract(v2, 0);
    gmx_simd_storeu_d(base + align * offset[1], t2);
    base[align * offset[1] + 2]  = vec_extract(v2, 1);
}

/* gcc-4.9 does not detect that vec_extract() uses its argument */
template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_d_ibm_vsx(double *                      base,
                                           const gmx_int32_t             offset[],
                                           gmx_simd_double_t             v0,
                                           gmx_simd_double_t             v1,
                                           gmx_simd_double_t gmx_unused  v2)
{
    gmx_simd_double_t t1, t2, t3, t4;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 8 == 0);
#endif

    t1  = vec_mergeh(v0, v1);
    t2  = vec_mergel(v0, v1);

    t3 = gmx_simd_loadu_d(base + align * offset[0]);
    t4 = gmx_simd_loadu_d(base + align * offset[1]);

    t3                           = vec_add(t3, t1);
    gmx_simd_storeu_d(base + align * offset[0], t3);
    base[align * offset[0] + 2] += vec_extract(v2, 0);
    t4                           = vec_add(t4, t2);
    gmx_simd_storeu_d(base + align * offset[1], t4);
    base[align * offset[1] + 2] += vec_extract(v2, 1);
}

/* gcc-4.9 does not detect that vec_extract() uses its argument */
template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_d_ibm_vsx(double *                      base,
                                           const gmx_int32_t             offset[],
                                           gmx_simd_double_t             v0,
                                           gmx_simd_double_t             v1,
                                           gmx_simd_double_t gmx_unused  v2)
{
    gmx_simd_double_t t1, t2, t3, t4;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 8 == 0);
#endif

    t1  = vec_mergeh(v0, v1);
    t2  = vec_mergel(v0, v1);

    t3 = gmx_simd_loadu_d(base + align * offset[0]);
    t4 = gmx_simd_loadu_d(base + align * offset[1]);

    t3                           = vec_sub(t3, t1);
    gmx_simd_storeu_d(base + align * offset[0], t3);
    base[align * offset[0] + 2] -= vec_extract(v2, 0);
    t4                           = vec_sub(t4, t2);
    gmx_simd_storeu_d(base + align * offset[1], t4);
    base[align * offset[1] + 2] -= vec_extract(v2, 1);
}

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_d_ibm_vsx(gmx_simd_double_t    scalar,
                                              gmx_simd_double_t   &triplets0,
                                              gmx_simd_double_t   &triplets1,
                                              gmx_simd_double_t   &triplets2)
{
    triplets0 = vec_mergeh(scalar, scalar);
    triplets1 = scalar;
    triplets2 = vec_mergel(scalar, scalar);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_ibm_vsx(const double *       base,
                                                   gmx_simd_dint32_t    offset,
                                                   gmx_simd_double_t   &v0,
                                                   gmx_simd_double_t   &v1,
                                                   gmx_simd_double_t   &v2,
                                                   gmx_simd_double_t   &v3)
{
    int               idx0, idx1;
    gmx_simd_double_t t1, t2, t3, t4;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);
#endif

    idx0 = gmx_simd_extract_di(offset, 0);
    idx1 = gmx_simd_extract_di(offset, 1);

    if (align == 1)
    {
        t1  = gmx_simd_load_d(base + idx0);
        t2  = gmx_simd_load_d(base + idx1);
        t3  = gmx_simd_load_d(base + idx0 + 2);
        t4  = gmx_simd_load_d(base + idx1 + 2);
    }
    else
    {
        t1  = gmx_simd_load_d(base + align * idx0);
        t2  = gmx_simd_load_d(base + align * idx1);
        t3  = gmx_simd_load_d(base + align * idx0 + 2);
        t4  = gmx_simd_load_d(base + align * idx1 + 2);
    }
    v0  = vec_mergeh(t1, t2);
    v1  = vec_mergel(t1, t2);
    v2  = vec_mergeh(t3, t4);
    v3  = vec_mergel(t3, t4);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_ibm_vsx(const double *       base,
                                                   gmx_simd_dint32_t    offset,
                                                   gmx_simd_double_t   &v0,
                                                   gmx_simd_double_t   &v1)
{
    int               idx0, idx1;
    gmx_simd_double_t t1, t2;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)base % 16 == 0);
    assert(align % 2 == 0);
#endif

    idx0 = gmx_simd_extract_di(offset, 0);
    idx1 = gmx_simd_extract_di(offset, 1);

    if (align == 1)
    {
        t1  = gmx_simd_load_d(base + idx0);
        t2  = gmx_simd_load_d(base + idx1);
    }
    else
    {
        t1  = gmx_simd_load_d(base + align * idx0);
        t2  = gmx_simd_load_d(base + align * idx1);
    }
    v0  = vec_mergeh(t1, t2);
    v1  = vec_mergel(t1, t2);
}

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_d_ibm_vsx(const double *       base,
                                                    gmx_simd_dint32_t    offset,
                                                    gmx_simd_double_t   &v0,
                                                    gmx_simd_double_t   &v1)
{
    /* Technically VSX uses the same instruction for aligned and unaligned
     * loads, but due to compiler bugs we need to use different function
     * names for now, and thus we cannot merge this with the corresponding
     * aligned-load function.
     */
    int               idx0, idx1;
    gmx_simd_double_t t1, t2;

    idx0 = gmx_simd_extract_di(offset, 0);
    idx1 = gmx_simd_extract_di(offset, 1);

    if (align == 1)
    {
        t1  = gmx_simd_loadu_d(base + idx0);
        t2  = gmx_simd_loadu_d(base + idx1);
    }
    else
    {
        t1  = gmx_simd_loadu_d(base + align * idx0);
        t2  = gmx_simd_loadu_d(base + align * idx1);
    }
    v0  = vec_mergeh(t1, t2);
    v1  = vec_mergel(t1, t2);
}

static gmx_inline double
gmx_simd_reduce_incr_4_return_sum_d_ibm_vsx(double *           m,
                                            gmx_simd_double_t  v0,
                                            gmx_simd_double_t  v1,
                                            gmx_simd_double_t  v2,
                                            gmx_simd_double_t  v3)
{
    gmx_simd_double_t t1, t2, t3, t4;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)m % 16 == 0);
#endif

    t1 = vec_mergeh(v0, v1);
    t2 = vec_mergel(v0, v1);
    t3 = vec_mergeh(v2, v3);
    t4 = vec_mergel(v2, v3);

    t1 = vec_add(t1, t2);
    t3 = vec_add(t3, t4);

    *(__vector double *)(m)   += t1;
    *(__vector double *)(m+2) += t3;

    t1 = vec_add(t1, t3);
    return gmx_simd_reduce_d_ibm_vsx(t1);
}
#endif /* __cplusplus */


#endif /* GMX_SIMD_IMPLEMENTATION_IBM_VSX_UTIL_DOUBLE_H */
