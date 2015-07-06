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

#ifndef GMX_SIMD_IMPL_SPARC64_HPC_ACE_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_SPARC64_HPC_ACE_UTIL_FLOAT_H

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
#include "impl_sparc64_hpc_ace_util_double.h"

/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_f             gmx_simd_gather_load_transpose_f_sparc64_hpc_ace
#define gmx_simd_best_pair_alignment_f               gmx_simd_best_pair_alignment_f_sparc64_hpc_ace
#define gmx_simd_gather_loadu_transpose_f            gmx_simd_gather_loadu_transpose_f_sparc64_hpc_ace
#define gmx_simd_transpose_scatter_storeu_f          gmx_simd_transpose_scatter_storeu_f_sparc64_hpc_ace
#define gmx_simd_transpose_scatter_incru_f           gmx_simd_transpose_scatter_incru_f_sparc64_hpc_ace
#define gmx_simd_transpose_scatter_decru_f           gmx_simd_transpose_scatter_decru_f_sparc64_hpc_ace
#define gmx_simd_expand_scalars_to_triplets_f        gmx_simd_expand_scalars_to_triplets_d_sparc64_hpc_ace
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_sparc64_hpc_ace
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_sparc64_hpc_ace
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_sparc64_hpc_ace

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

#define GMX_FJSP_TRANSPOSE2_V2R8(row0, row1) {                   \
        _fjsp_v2r8 gmx_fjsp_t1 = row0;                           \
        row0           = _fjsp_unpacklo_v2r8(row0, row1);        \
        row1           = _fjsp_unpackhi_v2r8(gmx_fjsp_t1, row1); \
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_sparc64_hpc_ace(const float *         base,
                                                 const gmx_int32_t     offset[],
                                                 gmx_simd_float_t     &v0,
                                                 gmx_simd_float_t     &v1,
                                                 gmx_simd_float_t     &v2,
                                                 gmx_simd_float_t     &v3)
{
    _fjsp_v2r8 t1, t2, t3, t4;
    t1  = gmx_simd_load_f(base + align * offset[0]);
    t2  = gmx_simd_load_f(base + align * offset[1]);
    t3  = gmx_simd_load_f(base + align * offset[0] + 2);
    t4  = gmx_simd_load_f(base + align * offset[1] + 2);
    v0  = _fjsp_unpacklo_v2r8(t1, t2);
    v1  = _fjsp_unpackhi_v2r8(t1, t2);
    v2  = _fjsp_unpacklo_v2r8(t3, t4);
    v3  = _fjsp_unpackhi_v2r8(t3, t4);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_sparc64_hpc_ace(const float *         base,
                                                 const gmx_int32_t     offset[],
                                                 gmx_simd_float_t     &v0,
                                                 gmx_simd_float_t     &v1)
{
    _fjsp_v2r8 t1, t2;
    t1  = gmx_simd_load_f(base + align * offset[0]);
    t2  = gmx_simd_load_f(base + align * offset[1]);
    v0  = _fjsp_unpacklo_v2r8(t1, t2);
    v1  = _fjsp_unpackhi_v2r8(t1, t2);
}

static const int gmx_simd_best_pair_alignment_f_sparc64_hpc_ace = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_f_sparc64_hpc_ace(const float *         base,
                                                  const gmx_int32_t     offset[],
                                                  gmx_simd_float_t     &v0,
                                                  gmx_simd_float_t     &v1,
                                                  gmx_simd_float_t     &v2)
{
    _fjsp_v2r8 t1, t2, t3, t4;
    /* Load elements 1+2 */
    t1           = gmx_simd_load_f(base + align * offset[0]);
    t2           = gmx_simd_load_f(base + align * offset[1]);
    /* We cannot load a single float (element 3), so load overlapping (elements 2+3) */
    t3           = gmx_simd_load_f(_fjsp_setzero_v2r8(), base + align * offset[0] + 1);
    t4           = gmx_simd_load_f(_fjsp_setzero_v2r8(), base + align * offset[1] + 1);
    GMX_FJSP_TRANSPOSE2_V2R8(t1, t2);
    v0           = t1;
    v1           = t2;
    v2           = _fjsp_unpackhi_v2r8(t3, t4);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_f_sparc64_hpc_ace(double *              base,
                                                    const gmx_int32_t     offset[],
                                                    gmx_simd_float_t      v0,
                                                    gmx_simd_float_t      v1,
                                                    gmx_simd_float_t      v2)
{
    float __attribute__ ((aligned (16))) f0[2], f1[2], f2[2];

    *(_fjsp_v2r4 *)f0 = _fjsp_dtos_v2r4(v0);
    *(_fjsp_v2r4 *)f1 = _fjsp_dtos_v2r4(v1);
    *(_fjsp_v2r4 *)f2 = _fjsp_dtos_v2r4(v2);

    base[align * offset[0]    ] = f0[0];
    base[align * offset[0] + 1] = f1[0];
    base[align * offset[0] + 2] = f2[0];
    base[align * offset[1]    ] = f0[1];
    base[align * offset[1] + 1] = f1[1];
    base[align * offset[1] + 2] = f2[1];
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_f_sparc64_hpc_ace(float *               base,
                                                   const gmx_int32_t     offset[],
                                                   gmx_simd_float_t      v0,
                                                   gmx_simd_float_t      v1,
                                                   gmx_simd_float_t      v2)
{
    float __attribute__ ((aligned (16))) f0[2], f1[2], f2[2];

    *(_fjsp_v2r4 *)f0 = _fjsp_dtos_v2r4(v0);
    *(_fjsp_v2r4 *)f1 = _fjsp_dtos_v2r4(v1);
    *(_fjsp_v2r4 *)f2 = _fjsp_dtos_v2r4(v2);

    base[align * offset[0]    ] += f0[0];
    base[align * offset[0] + 1] += f1[0];
    base[align * offset[0] + 2] += f2[0];
    base[align * offset[1]    ] += f0[1];
    base[align * offset[1] + 1] += f1[1];
    base[align * offset[1] + 2] += f2[1];
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_f_sparc64_hpc_ace(double *              base,
                                                   const gmx_int32_t     offset[],
                                                   gmx_simd_float_t      v0,
                                                   gmx_simd_float_t      v1,
                                                   gmx_simd_float_t      v2)
{
    float __attribute__ ((aligned (16))) f0[2], f1[2], f2[2];

    *(_fjsp_v2r4 *)f0 = _fjsp_dtos_v2r4(v0);
    *(_fjsp_v2r4 *)f1 = _fjsp_dtos_v2r4(v1);
    *(_fjsp_v2r4 *)f2 = _fjsp_dtos_v2r4(v2);

    base[align * offset[0]    ] -= f0[0];
    base[align * offset[0] + 1] -= f1[0];
    base[align * offset[0] + 2] -= f2[0];
    base[align * offset[1]    ] -= f0[1];
    base[align * offset[1] + 1] -= f1[1];
    base[align * offset[1] + 2] -= f2[1];
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_sparc64_hpc_ace(const float *        base,
                                                           gmx_simd_fint32_t    offset,
                                                           gmx_simd_float_t    &v0,
                                                           gmx_simd_float_t    &v1,
                                                           gmx_simd_float_t    &v2,
                                                           gmx_simd_float_t    &v3)
{
    _fjsp_v2r8 t1, t2, t3, t4;
    long long int __attribute__ ((aligned (16))) itmp[2];

    _fjsp_store_v2r8( (double *) itmp, offset );

    t1  = gmx_simd_load_f(base + align * itmp[0]);
    t2  = gmx_simd_load_f(base + align * itmp[1]);
    t3  = gmx_simd_load_f(base + align * itmp[0] + 2);
    t4  = gmx_simd_load_f(base + align * itmp[1] + 2);
    v0  = _fjsp_unpacklo_v2r8(t1, t2);
    v1  = _fjsp_unpackhi_v2r8(t1, t2);
    v2  = _fjsp_unpacklo_v2r8(t3, t4);
    v3  = _fjsp_unpackhi_v2r8(t3, t4);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_sparc64_hpc_ace(const float *        base,
                                                           gmx_simd_fint32_t    offset,
                                                           gmx_simd_float_t    &v0,
                                                           gmx_simd_float_t    &v1)
{
    _fjsp_v2r8 t1, t2;
    long long int __attribute__ ((aligned (16))) itmp[2];

    _fjsp_store_v2r8( (double *) ioffset, offset );

    t1  = gmx_simd_load_f(base + align * itmp[0]);
    t2  = gmx_simd_load_f(base + align * itmp[1]);
    v0  = _fjsp_unpacklo_v2r8(t1, t2);
    v1  = _fjsp_unpackhi_v2r8(t1, t2);
}

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_f_sparc64_hpc_ace(const float *        base,
                                                            gmx_simd_dint32_t    offset,
                                                            gmx_simd_float_t    &v0,
                                                            gmx_simd_float_t    &v1)
{
    gmx_simd_gather_load_bysimdint_transpose_f_sparc64_hpc_ace<align>(base, offset, v0, v1);
}

static gmx_inline double
gmx_simd_reduce_d_sparc64_hpc_ace(gmx_simd_double_t x)
{
    double d;
    x = _fjsp_add_v2r8(x, _fjsp_unpackhi_v2r8(x, x));
    _fjsp_storel_v2r8(&d, x);
    return d;
}

static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_f_sparc64_hpc_ace(float *            m,
                                                    gmx_simd_float_t   v0,
                                                    gmx_simd_float_t   v1,
                                                    gmx_simd_float_t   v2,
                                                    gmx_simd_float_t   v3)
{
    _fjsp_v2r8 t1, t2, t3, t4;

    t1 = _fjsp_unpacklo_v2r8(v0, v1);
    t2 = _fjsp_unpackhi_v2r8(v0, v1);
    t3 = _fjsp_unpacklo_v2r8(v2, v3);
    t4 = _fjsp_unpackhi_v2r8(v2, v3);

    t1 = _fjsp_add_v2r8(t1, t2);
    t3 = _fjsp_add_v2r8(t3, t4);

    t2 = _fjsp_add_v2r8(t1, gmx_simd_load_f(m));
    t4 = _fjsp_add_v2r8(t3, gmx_simd_load_f(m + 2));
    gmx_simd_store_f(m, t2);
    gmx_simd_store_f(m + 2, t4);

    t1 = _fjsp_add_v2r8(t1, t3);
    return gmx_simd_reduce_d_sparc64_hpc_ace(t1);
}
#endif /* __cplusplus */


#endif /* GMX_SIMD_IMPL_SPARC64_HPC_ACE_UTIL_FLOAT_H */
