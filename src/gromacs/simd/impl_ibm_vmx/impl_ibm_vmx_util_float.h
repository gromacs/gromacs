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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VMX_UTIL_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VMX_UTIL_FLOAT_H

#if !defined(__ibmxl__) && !defined(__xlC__)
/* xlc-13 (at least) seems to be buggy with the asserts at high optimization */
#    include <assert.h>
#endif

#include <math.h>

#include <altivec.h>

#include "impl_ibm_vmx_common.h"

#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__)
#    define gmx_vmx_bool __bool
#    undef  bool
#else
#    define gmx_vmx_bool bool
#endif

/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_f             gmx_simd_gather_load_transpose_f_ibm_vmx
#define gmx_simd_best_pair_alignment_f               gmx_simd_best_pair_alignment_f_ibm_vmx
#define gmx_simd_gather_loadu_transpose_f            gmx_simd_gather_loadu_transpose_f_ibm_vmx
#define gmx_simd_transpose_scatter_storeu_f          gmx_simd_transpose_scatter_storeu_f_ibm_vmx
#define gmx_simd_transpose_scatter_incru_f           gmx_simd_transpose_scatter_incru_f_ibm_vmx
#define gmx_simd_transpose_scatter_decru_f           gmx_simd_transpose_scatter_decru_f_ibm_vmx
#define gmx_simd_expand_scalars_to_triplets_f        gmx_simd_expand_scalars_to_triplets_f_ibm_vmx
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_ibm_vmx
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_ibm_vmx
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_ibm_vmx

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

/* Internal macro for transposes */
#define GMX_VMX_TRANSPOSE4(v0, v1, v2, v3)               \
    {                                                    \
        __vector float gmx_vmx_t0 = vec_mergeh(v0, v2);  \
        __vector float gmx_vmx_t1 = vec_mergel(v0, v2);  \
        __vector float gmx_vmx_t2 = vec_mergeh(v1, v3);  \
        __vector float gmx_vmx_t3 = vec_mergel(v1, v3);  \
        v0 = vec_mergeh(gmx_vmx_t0, gmx_vmx_t2);         \
        v1 = vec_mergel(gmx_vmx_t0, gmx_vmx_t2);         \
        v2 = vec_mergeh(gmx_vmx_t1, gmx_vmx_t3);         \
        v3 = vec_mergel(gmx_vmx_t1, gmx_vmx_t3);         \
    }


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_ibm_vmx(const float *        base,
                                         const gmx_int32_t    offset[],
                                         gmx_simd_float_t    &v0,
                                         gmx_simd_float_t    &v1,
                                         gmx_simd_float_t    &v2,
                                         gmx_simd_float_t    &v3)
{
#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);
#endif

    v0  = gmx_simd_load_f( base + align * offset[0] );
    v1  = gmx_simd_load_f( base + align * offset[1] );
    v2  = gmx_simd_load_f( base + align * offset[2] );
    v3  = gmx_simd_load_f( base + align * offset[3] );
    GMX_VMX_TRANSPOSE4(v0, v1, v2, v3);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_ibm_vmx(const float *        base,
                                         const gmx_int32_t    offset[],
                                         gmx_simd_float_t    &v0,
                                         gmx_simd_float_t    &v1)
{
    gmx_simd_float_t       t0, t1, t2, t3, t4, t5, t6, t7;
    __vector unsigned char p0, p1, p2, p3;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 8 == 0);
    assert(align % 2 == 0);
#endif

    if ((align & 0x3) == 0)
    {
        gmx_simd_gather_load_transpose_f_ibm_vmx<align>(base, offset, v0, v1, t2, t3);
    }
    else
    {
        /* This is REALLY slow, since we have no choice but to load individual
         * elements when we cannot guarantee that we can access beyond the end of
         * the memory. Fortunately, 99% of the usage should be the aligned-to-4
         * case above instead.
         */
        t0 = vec_lde(0, base + align * offset[0]);
        t1 = vec_lde(0, base + align * offset[1]);
        t2 = vec_lde(0, base + align * offset[2]);
        t3 = vec_lde(0, base + align * offset[3]);
        p0 = vec_lvsl(0, base + align * offset[0]);
        p1 = vec_lvsl(0, base + align * offset[1]);
        p2 = vec_lvsl(0, base + align * offset[2]);
        p3 = vec_lvsl(0, base + align * offset[3]);
        t0 = vec_perm(t0, t0, p0);
        t1 = vec_perm(t1, t1, p1);
        t2 = vec_perm(t2, t2, p2);
        t3 = vec_perm(t3, t3, p3);
        t0 = vec_mergeh(t0, t2);
        t1 = vec_mergeh(t1, t3);
        v0 = vec_mergeh(t0, t1);

        t4 = vec_lde(0, base + align * offset[0] + 1);
        t5 = vec_lde(0, base + align * offset[1] + 1);
        t6 = vec_lde(0, base + align * offset[2] + 1);
        t7 = vec_lde(0, base + align * offset[3] + 1);
        p0 = vec_lvsl(0, base + align * offset[0] + 1);
        p1 = vec_lvsl(0, base + align * offset[1] + 1);
        p2 = vec_lvsl(0, base + align * offset[2] + 1);
        p3 = vec_lvsl(0, base + align * offset[3] + 1);
        t4 = vec_perm(t4, t4, p0);
        t5 = vec_perm(t5, t5, p1);
        t6 = vec_perm(t6, t6, p2);
        t7 = vec_perm(t7, t7, p3);
        t4 = vec_mergeh(t4, t6);
        t5 = vec_mergeh(t5, t7);
        v1 = vec_mergeh(t4, t5);
    }
}

static const int gmx_simd_best_pair_alignment_f_ibm_vmx = 4;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_f_ibm_vmx(const float *        base,
                                          const gmx_int32_t    offset[],
                                          gmx_simd_float_t    &v0,
                                          gmx_simd_float_t    &v1,
                                          gmx_simd_float_t    &v2)
{
    gmx_simd_float_t       t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
    __vector unsigned char p0, p1, p2, p3;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
#endif

    if ((align & 0x3) == 0)
    {
        gmx_simd_gather_load_transpose_f_ibm_vmx<align>(base, offset, v0, v1, v2, t3);
    }
    else
    {
        /* This is REALLY slow, since we have no choice but to load individual
         * elements when we cannot guarantee that we can access beyond the end of
         * the memory. Unfortunately this is likely the most common case.
         */
        t0  = vec_lde(0, base + align * offset[0]);
        t1  = vec_lde(0, base + align * offset[1]);
        t2  = vec_lde(0, base + align * offset[2]);
        t3  = vec_lde(0, base + align * offset[3]);
        p0  = vec_lvsl(0, base + align * offset[0]);
        p1  = vec_lvsl(0, base + align * offset[1]);
        p2  = vec_lvsl(0, base + align * offset[2]);
        p3  = vec_lvsl(0, base + align * offset[3]);
        t0  = vec_perm(t0, t0, p0);
        t1  = vec_perm(t1, t1, p1);
        t2  = vec_perm(t2, t2, p2);
        t3  = vec_perm(t3, t3, p3);
        t0  = vec_mergeh(t0, t2);
        t1  = vec_mergeh(t1, t3);
        v0  = vec_mergeh(t0, t1);

        t4  = vec_lde(0, base + align * offset[0] + 1);
        t5  = vec_lde(0, base + align * offset[1] + 1);
        t6  = vec_lde(0, base + align * offset[2] + 1);
        t7  = vec_lde(0, base + align * offset[3] + 1);
        p0  = vec_lvsl(0, base + align * offset[0] + 1);
        p1  = vec_lvsl(0, base + align * offset[1] + 1);
        p2  = vec_lvsl(0, base + align * offset[2] + 1);
        p3  = vec_lvsl(0, base + align * offset[3] + 1);
        t4  = vec_perm(t4, t4, p0);
        t5  = vec_perm(t5, t5, p1);
        t6  = vec_perm(t6, t6, p2);
        t7  = vec_perm(t7, t7, p3);
        t4  = vec_mergeh(t4, t6);
        t5  = vec_mergeh(t5, t7);
        v1  = vec_mergeh(t4, t5);

        t8  = vec_lde(0, base + align * offset[0] + 2);
        t9  = vec_lde(0, base + align * offset[1] + 2);
        t10 = vec_lde(0, base + align * offset[2] + 2);
        t11 = vec_lde(0, base + align * offset[3] + 2);
        p0  = vec_lvsl(0, base + align * offset[0] + 2);
        p1  = vec_lvsl(0, base + align * offset[1] + 2);
        p2  = vec_lvsl(0, base + align * offset[2] + 2);
        p3  = vec_lvsl(0, base + align * offset[3] + 2);
        t8  = vec_perm(t8, t8, p0);
        t9  = vec_perm(t9, t9, p1);
        t10 = vec_perm(t10, t10, p2);
        t11 = vec_perm(t11, t11, p3);
        t8  = vec_mergeh(t8, t10);
        t9  = vec_mergeh(t9, t11);
        v2  = vec_mergeh(t8, t9);
    }
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_f_ibm_vmx(float *                        base,
                                            const gmx_int32_t              offset[],
                                            gmx_simd_float_t               v0,
                                            gmx_simd_float_t               v1,
                                            gmx_simd_float_t gmx_unused    v2)
{
    __vector float         t3;
    __vector unsigned char p0, p1, p2, p3;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
#endif

    t3 = (__vector float)vec_splat_u32(0);
    GMX_VMX_TRANSPOSE4(v0, v1, v2, t3);

    p0 = vec_lvsr(0, base + align * offset[0]);
    p1 = vec_lvsr(0, base + align * offset[1]);
    p2 = vec_lvsr(0, base + align * offset[2]);
    p3 = vec_lvsr(0, base + align * offset[3]);

    v0 = vec_perm(v0, v0, p0);
    v1 = vec_perm(v1, v1, p1);
    v2 = vec_perm(v2, v2, p2);
    t3 = vec_perm(t3, t3, p3);

    vec_ste(v0, 0, base + align * offset[0]);
    vec_ste(v0, 4, base + align * offset[0]);
    vec_ste(v0, 8, base + align * offset[0]);
    vec_ste(v1, 0, base + align * offset[1]);
    vec_ste(v1, 4, base + align * offset[1]);
    vec_ste(v1, 8, base + align * offset[1]);
    vec_ste(v2, 0, base + align * offset[2]);
    vec_ste(v2, 4, base + align * offset[2]);
    vec_ste(v2, 8, base + align * offset[2]);
    vec_ste(t3, 0, base + align * offset[3]);
    vec_ste(t3, 4, base + align * offset[3]);
    vec_ste(t3, 8, base + align * offset[3]);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_f_ibm_vmx(float *                       base,
                                           const gmx_int32_t             offset[],
                                           gmx_simd_float_t              v0,
                                           gmx_simd_float_t              v1,
                                           gmx_simd_float_t gmx_unused   v2)
{
    gmx_simd_float_t             t0, t1, t2;

    gmx_simd_gather_loadu_transpose_f_ibm_vmx<align>(base, offset, t0, t1, t2);
    t0 = vec_add(t0, v0);
    t1 = vec_add(t1, v1);
    t2 = vec_add(t2, v2);
    gmx_simd_transpose_scatter_storeu_f_ibm_vmx<align>(base, offset, t0, t1, t2);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_f_ibm_vmx(float *                       base,
                                           const gmx_int32_t             offset[],
                                           gmx_simd_float_t              v0,
                                           gmx_simd_float_t              v1,
                                           gmx_simd_float_t gmx_unused   v2)
{
    gmx_simd_float_t             t0, t1, t2;

    gmx_simd_gather_loadu_transpose_f_ibm_vmx<align>(base, offset, t0, t1, t2);
    t0 = vec_sub(t0, v0);
    t1 = vec_sub(t1, v1);
    t2 = vec_sub(t2, v2);
    gmx_simd_transpose_scatter_storeu_f_ibm_vmx<align>(base, offset, t0, t1, t2);
}

static gmx_inline void
gmx_simd_expand_scalars_to_triplets_f_ibm_vmx(gmx_simd_float_t    scalar,
                                              gmx_simd_float_t   &triplets0,
                                              gmx_simd_float_t   &triplets1,
                                              gmx_simd_float_t   &triplets2)
{
    __vector unsigned char perm0 = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7 };
    __vector unsigned char perm1 = { 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11, 8, 9, 10, 11 };
    __vector unsigned char perm2 = { 8, 9, 10, 11, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15 };

    triplets0 = vec_perm(scalar, scalar, perm0);
    triplets1 = vec_perm(scalar, scalar, perm1);
    triplets2 = vec_perm(scalar, scalar, perm2);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_ibm_vmx(const float *       base,
                                                   gmx_simd_fint32_t   offset,
                                                   gmx_simd_float_t   &v0,
                                                   gmx_simd_float_t   &v1,
                                                   gmx_simd_float_t   &v2,
                                                   gmx_simd_float_t   &v3)
{
#ifdef __GNUC__
    gmx_int32_t __attribute__ ((aligned (16))) ioffset[4];
    gmx_simd_store_fi(ioffset, offset);
#else
    struct {
        __vector signed int simd; gmx_int32_t i[4];
    } conv;
    gmx_int32_t *ioffset = conv.i;
    gmx_simd_store_fi((int *)&conv.simd, offset);
#endif
    gmx_simd_gather_load_transpose_f_ibm_vmx<align>(base, ioffset, v0, v1, v2, v3);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_ibm_vmx(const float *       base,
                                                   gmx_simd_fint32_t   offset,
                                                   gmx_simd_float_t   &v0,
                                                   gmx_simd_float_t   &v1)
{
#ifdef __GNUC__
    gmx_int32_t __attribute__ ((aligned (16))) ioffset[4];
    gmx_simd_store_fi(ioffset, offset);
#else
    struct {
        __vector signed int simd; gmx_int32_t i[4];
    } conv;
    gmx_int32_t *ioffset = conv.i;
    gmx_simd_store_fi((int *)&conv.simd, offset);
#endif
    gmx_simd_gather_load_transpose_f_ibm_vmx<align>(base, ioffset, v0, v1);
}

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_f_ibm_vmx(const float *       base,
                                                    gmx_simd_fint32_t   offset,
                                                    gmx_simd_float_t   &v0,
                                                    gmx_simd_float_t   &v1)
{
    gmx_simd_gather_load_bysimdint_transpose_f_ibm_vmx<align>(base, offset, v0, v1);
}


static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_f_ibm_vmx(float *           m,
                                            gmx_simd_float_t  v0,
                                            gmx_simd_float_t  v1,
                                            gmx_simd_float_t  v2,
                                            gmx_simd_float_t  v3)
{
#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)m % 16 == 0);
#endif

    GMX_VMX_TRANSPOSE4(v0, v1, v2, v3);
    v0 = vec_add(v0, v1);
    v2 = vec_add(v2, v3);
    v0 = vec_add(v0, v2);
    v2 = vec_add(v0, gmx_simd_load_f(m));
    gmx_simd_store_f(m, v2);

    return gmx_simd_reduce_f_ibm_vmx(v0);
}
#endif /* __cplusplus */


#endif /* GMX_SIMD_IMPLEMENTATION_IBM_VMX_UTIL_FLOAT_H */
