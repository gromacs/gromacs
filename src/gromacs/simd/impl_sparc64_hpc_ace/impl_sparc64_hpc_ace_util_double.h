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

#ifndef GMX_SIMD_IMPL_SPARC64_HPC_ACE_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_SPARC64_HPC_ACE_UTIL_DOUBLE_H

#include "config.h"

/* Fujitsu header borrows the name from SSE2, since some instructions have aliases.
 * Environment/compiler version GM-1.2.0-17 seems to be buggy; when -Xg is
 * defined to enable GNUC extensions, this sets _ISOC99_SOURCE, which in
 * turn causes all intrinsics to be declared inline _instead_ of static. This
 * leads to duplicate symbol errors at link time.
 * To work around this we unset this before including the HPC-ACE header, and
 * reset the value afterwards.
 */
#ifdef _ISOC99_SOURCE
#    undef _ISOC99_SOURCE
#    define SAVE_ISOC99_SOURCE
#endif

#include <emmintrin.h>

#ifdef SAVE_ISOC99_SOURCE
#    define _ISOC99_SOURCE
#    undef SAVE_ISOC99_SOURCE
#endif

#include <math.h>

#include "impl_sparc64_hpc_ace_common.h"

/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_d             gmx_simd_gather_load_transpose_d_sparc64_hpc_ace
#define gmx_simd_best_pair_alignment_d               gmx_simd_best_pair_alignment_d_sparc64_hpc_ace
#define gmx_simd_gather_loadu_transpose_d            gmx_simd_gather_loadu_transpose_d_sparc64_hpc_ace
#define gmx_simd_transpose_scatter_storeu_d          gmx_simd_transpose_scatter_storeu_d_sparc64_hpc_ace
#define gmx_simd_transpose_scatter_incru_d           gmx_simd_transpose_scatter_incru_d_sparc64_hpc_ace
#define gmx_simd_transpose_scatter_decru_d           gmx_simd_transpose_scatter_decru_d_sparc64_hpc_ace
#define gmx_simd_expand_scalars_to_triplets_d        gmx_simd_expand_scalars_to_triplets_d_sparc64_hpc_ace
#define gmx_simd_gather_load_bysimdint_transpose_d   gmx_simd_gather_load_bysimdint_transpose_d_sparc64_hpc_ace
#define gmx_simd_gather_loadu_bysimdint_transpose_d  gmx_simd_gather_loadu_bysimdint_transpose_d_sparc64_hpc_ace
#define gmx_simd_reduce_incr_4_return_sum_d          gmx_simd_reduce_incr_4_return_sum_d_sparc64_hpc_ace

/****************************************************
 * Double precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d_sparc64_hpc_ace(const double *        base,
                                                 const gmx_int32_t     offset[],
                                                 gmx_simd_double_t    &v0,
                                                 gmx_simd_double_t    &v1,
                                                 gmx_simd_double_t    &v2,
                                                 gmx_simd_double_t    &v3)
{
    _fjsp_v2r8 t1, t2, t3, t4;
    t1  = _fjsp_load_v2r8(base + align * offset[0]);
    t2  = _fjsp_load_v2r8(base + align * offset[1]);
    t3  = _fjsp_load_v2r8(base + align * offset[0] + 2);
    t4  = _fjsp_load_v2r8(base + align * offset[1] + 2);
    v0  = _fjsp_unpacklo_v2r8(t1, t2);
    v1  = _fjsp_unpackhi_v2r8(t1, t2);
    v2  = _fjsp_unpacklo_v2r8(t3, t4);
    v3  = _fjsp_unpackhi_v2r8(t3, t4);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d_sparc64_hpc_ace(const double *        base,
                                                 const gmx_int32_t     offset[],
                                                 gmx_simd_double_t    &v0,
                                                 gmx_simd_double_t    &v1)
{
    _fjsp_v2r8 t1, t2;
    t1  = _fjsp_load_v2r8(base + align * offset[0]);
    t2  = _fjsp_load_v2r8(base + align * offset[1]);
    v0  = _fjsp_unpacklo_v2r8(t1, t2);
    v1  = _fjsp_unpackhi_v2r8(t1, t2);
}

static const int gmx_simd_best_pair_alignment_d_sparc64_hpc_ace = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_d_sparc64_hpc_ace(const double *        base,
                                                  const gmx_int32_t     offset[],
                                                  gmx_simd_double_t    &v0,
                                                  gmx_simd_double_t    &v1,
                                                  gmx_simd_double_t    &v2)
{
    _fjsp_v2r8 t1, t2, t3, t4;
    t1           = _fjsp_load_v2r8(base + align * offset[0]);
    t2           = _fjsp_load_v2r8(base + align * offset[1]);
    t3           = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), base + align * offset[0] + 2);
    t4           = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), base + align * offset[1] + 2);
    GMX_FJSP_TRANSPOSE2_V2R8(t1, t2);
    v0           = t1;
    v1           = t2;
    v2           = _fjsp_unpacklo_v2r8(t3, t4);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_d_sparc64_hpc_ace(double *              base,
                                                    const gmx_int32_t     offset[],
                                                    gmx_simd_double_t     v0,
                                                    gmx_simd_double_t     v1,
                                                    gmx_simd_double_t     v2)
{
    _fjsp_storel_v2r8(base + align * offset[0], v0);
    _fjsp_storel_v2r8(base + align * offset[0] + 1, v1);
    _fjsp_storel_v2r8(base + align * offset[0] + 2, v2);
    _fjsp_storeh_v2r8(base + align * offset[1], v0);
    _fjsp_storeh_v2r8(base + align * offset[1] + 1, v1);
    _fjsp_storeh_v2r8(base + align * offset[1] + 2, v2);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_d_sparc64_hpc_ace(double *              base,
                                                   const gmx_int32_t     offset[],
                                                   gmx_simd_double_t     v0,
                                                   gmx_simd_double_t     v1,
                                                   gmx_simd_double_t     v2)
{
    _fjsp_v2r8 t1, t2, t3, t4, t5, t6, t7;

    t1          = _fjsp_load_v2r8(base + align * offset[0]);
    t2          = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), base + align * offset[0]+2);
    t3          = _fjsp_load_v2r8(base + align * offset[1]);
    t4          = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), base + align * offset[1]+2);

    t5          = _fjsp_unpacklo_v2r8(v0, v1);
    t6          = _fjsp_unpackhi_v2r8(v0, v1);
    t7          = _fjsp_unpackhi_v2r8(v2, v2);

    t1          = _fjsp_add_v2r8(t1, t5);
    t2          = _fjsp_add_v2r8(t2, v2);

    t3          = _fjsp_add_v2r8(t3, t6);
    t4          = _fjsp_add_v2r8(t4, t7);

    _fjsp_storel_v2r8(base + align * offset[0], t1);
    _fjsp_storeh_v2r8(base + align * offset[0] + 1, t1);
    _fjsp_storel_v2r8(base + align * offset[0] + 2, t2);
    _fjsp_storel_v2r8(base + align * offset[1], t3);
    _fjsp_storeh_v2r8(base + align * offset[1] + 1, t3);
    _fjsp_storel_v2r8(base + align * offset[1] + 2, t4);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_d_sparc64_hpc_ace(double *              base,
                                                   const gmx_int32_t     offset[],
                                                   gmx_simd_double_t     v0,
                                                   gmx_simd_double_t     v1,
                                                   gmx_simd_double_t     v2)
{
    _fjsp_v2r8 t1, t2, t3, t4, t5, t6, t7;

    t1          = _fjsp_load_v2r8(base + align * offset[0]);
    t2          = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), base + align * offset[0]+2);
    t3          = _fjsp_load_v2r8(base + align * offset[1]);
    t4          = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(), base + align * offset[1]+2);

    t5          = _fjsp_unpacklo_v2r8(v0, v1);
    t6          = _fjsp_unpackhi_v2r8(v0, v1);
    t7          = _fjsp_unpackhi_v2r8(v2, v2);

    t1          = _fjsp_sub_v2r8(t1, t5);
    t2          = _fjsp_sub_v2r8(t2, v2);

    t3          = _fjsp_sub_v2r8(t3, t6);
    t4          = _fjsp_sub_v2r8(t4, t7);

    _fjsp_storel_v2r8(base + align * offset[0], t1);
    _fjsp_storeh_v2r8(base + align * offset[0] + 1, t1);
    _fjsp_storel_v2r8(base + align * offset[0] + 2, t2);
    _fjsp_storel_v2r8(base + align * offset[1], t3);
    _fjsp_storeh_v2r8(base + align * offset[1] + 1, t3);
    _fjsp_storel_v2r8(base + align * offset[1] + 2, t4);
}

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_d_sparc64_hpc_ace(gmx_simd_double_t    scalar,
                                                      gmx_simd_double_t   &triplets0,
                                                      gmx_simd_double_t   &triplets1,
                                                      gmx_simd_double_t   &triplets2)
{
    triplets0 = _fjsp_unpacklo_v2r8(scalar, scalar);
    triplets1 = scalar;
    triplets2 = _fjsp_unpackhi_v2r8(scalar, scalar);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_sparc64_hpc_ace(const double *       base,
                                                           gmx_simd_dint32_t    offset,
                                                           gmx_simd_double_t   &v0,
                                                           gmx_simd_double_t   &v1,
                                                           gmx_simd_double_t   &v2,
                                                           gmx_simd_double_t   &v3)
{
    _fjsp_v2r8 t1, t2, t3, t4;
    long long int __attribute__ ((aligned (16))) itmp[2];

    _fjsp_store_v2r8( (double *) itmp, offset );

    t1  = _fjsp_load_v2r8(base + align * itmp[0]);
    t2  = _fjsp_load_v2r8(base + align * itmp[1]);
    t3  = _fjsp_load_v2r8(base + align * itmp[0] + 2);
    t4  = _fjsp_load_v2r8(base + align * itmp[1] + 2);
    v0  = _fjsp_unpacklo_v2r8(t1, t2);
    v1  = _fjsp_unpackhi_v2r8(t1, t2);
    v2  = _fjsp_unpacklo_v2r8(t3, t4);
    v3  = _fjsp_unpackhi_v2r8(t3, t4);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d_sparc64_hpc_ace(const double *       base,
                                                           gmx_simd_dint32_t    offset,
                                                           gmx_simd_double_t   &v0,
                                                           gmx_simd_double_t   &v1)
{
    _fjsp_v2r8 t1, t2;
    long long int __attribute__ ((aligned (16))) itmp[2];

    _fjsp_store_v2r8( (double *) itmp, offset );

    t1  = _fjsp_load_v2r8(base + align * itmp[0]);
    t2  = _fjsp_load_v2r8(base + align * itmp[1]);
    v0  = _fjsp_unpacklo_v2r8(t1, t2);
    v1  = _fjsp_unpackhi_v2r8(t1, t2);
}

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_d_sparc64_hpc_ace(const double *       base,
                                                            gmx_simd_dint32_t    offset,
                                                            gmx_simd_double_t   &v0,
                                                            gmx_simd_double_t   &v1)
{
    gmx_simd_gather_load_bysimdint_transpose_d_sparc64_hpc_ace<align>(base, offset, v0, v1);
}

static gmx_inline double
gmx_simd_reduce_incr_4_return_sum_d_sparc64_hpc_ace(double *           m,
                                                    gmx_simd_double_t  v0,
                                                    gmx_simd_double_t  v1,
                                                    gmx_simd_double_t  v2,
                                                    gmx_simd_double_t  v3)
{
    _fjsp_v2r8 t1, t2, t3, t4;

    t1 = _fjsp_unpacklo_v2r8(v0, v1);
    t2 = _fjsp_unpackhi_v2r8(v0, v1);
    t3 = _fjsp_unpacklo_v2r8(v2, v3);
    t4 = _fjsp_unpackhi_v2r8(v2, v3);

    t1 = _fjsp_add_v2r8(t1, t2);
    t3 = _fjsp_add_v2r8(t3, t4);

    t2 = _fjsp_add_v2r8(t1, _fjsp_load_v2r8(m));
    t4 = _fjsp_add_v2r8(t3, _fjsp_load_v2r8(m + 2));
    _fjsp_store_v2r8(m, t2);
    _fjsp_store_v2r8(m + 2, t4);

    t1 = _fjsp_add_v2r8(t1, t3);
    return gmx_simd_reduce_d_sparc64_hpc_ace(t1);
}
#endif /* __cplusplus */


#endif /* GMX_SIMD_IMPL_SPARC64_HPC_ACE_UTIL_DOUBLE_H */
