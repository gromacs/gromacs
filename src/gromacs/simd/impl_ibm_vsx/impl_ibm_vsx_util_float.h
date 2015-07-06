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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VSX_UTIL_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VSX_UTIL_FLOAT_H

#if !defined(__ibmxl__) && !defined(__xlC__)
/* xlc-13 (at least) seems to be buggy with the asserts at high optimization */
#    include <assert.h>
#endif

#include <math.h>

#include <altivec.h>

#include "impl_ibm_vsx_common.h"
#include "impl_ibm_vsx_simd4_float.h"

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
#define gmx_simd_gather_load_transpose_f             gmx_simd_gather_load_transpose_f_ibm_vsx
#define gmx_simd_best_pair_alignment_f               gmx_simd_best_pair_alignment_f_ibm_vsx
#define gmx_simd_gather_loadu_transpose_f            gmx_simd_gather_loadu_transpose_f_ibm_vsx
#define gmx_simd_transpose_scatter_storeu_f          gmx_simd_transpose_scatter_storeu_f_ibm_vsx
#define gmx_simd_transpose_scatter_incru_f           gmx_simd_transpose_scatter_incru_f_ibm_vsx
#define gmx_simd_transpose_scatter_decru_f           gmx_simd_transpose_scatter_decru_f_ibm_vsx
#define gmx_simd_expand_scalars_to_triplets_f        gmx_simd_expand_scalars_to_triplets_f_ibm_vsx
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_ibm_vsx
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_ibm_vsx
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_ibm_vsx

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_ibm_vsx(const float *        base,
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
    GMX_VSX_TRANSPOSE4(v0, v1, v2, v3);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_ibm_vsx(const float *        base,
                                         const gmx_int32_t    offset[],
                                         gmx_simd_float_t    &v0,
                                         gmx_simd_float_t    &v1)
{
    gmx_simd_float_t t0, t1, t2, t3;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
    assert((size_t)base % 8 == 0);
    assert(align % 2 == 0);
#endif

    t0 = (__vector float)vec_splats(*(double *)(base + align * offset[0]));
    t1 = (__vector float)vec_splats(*(double *)(base + align * offset[1]));
    t2 = (__vector float)vec_splats(*(double *)(base + align * offset[2]));
    t3 = (__vector float)vec_splats(*(double *)(base + align * offset[3]));
    t0 = vec_mergeh(t0, t2);
    t1 = vec_mergeh(t1, t3);
    v0 = vec_mergeh(t0, t1);
    v1 = vec_mergel(t0, t1);
}

static const int gmx_simd_best_pair_alignment_f_ibm_vsx = 2;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_f_ibm_vsx(const float *        base,
                                          const gmx_int32_t    offset[],
                                          gmx_simd_float_t    &v0,
                                          gmx_simd_float_t    &v1,
                                          gmx_simd_float_t    &v2)
{
    gmx_simd_float_t             t1, t2, t3, t4, t5, t6, t7, t8;
    const __vector unsigned char perm_lo2hi = { 0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23 };
    const __vector unsigned char perm_hi2lo = { 24, 25, 26, 27, 28, 29, 30, 31, 8, 9, 10, 11, 12, 13, 14, 15 };

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
#endif

    /* The template conditional should be optimized away at compile time */
    if ( (align & 0x3) == 0)
    {
        gmx_simd_float_t t3;
        gmx_simd_gather_load_transpose_f_ibm_vsx<align>(base, offset, v0, v1, v2, t3);
    }
    else
    {
#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__) && ((__GNUC__ < 4) || ((__GNUC__ == 4) && (__GNUC_MINOR < 6)))
        /* GCC prior to 4.6 is buggy for loading pairs-of-floats-as-doubles, use simple splat instead */
        t1 = vec_mergeh(vec_splats(base[align * offset[0]]), vec_splats(base[align * offset[0] + 1]));
        t2 = vec_mergeh(vec_splats(base[align * offset[1]]), vec_splats(base[align * offset[1] + 1]));
        t3 = vec_mergeh(vec_splats(base[align * offset[2]]), vec_splats(base[align * offset[2] + 1]));
        t4 = vec_mergeh(vec_splats(base[align * offset[3]]), vec_splats(base[align * offset[3] + 1]));
#else
        t1 = (__vector float)vec_splats(*(double *)(base + align * offset[0]));
        t2 = (__vector float)vec_splats(*(double *)(base + align * offset[1]));
        t3 = (__vector float)vec_splats(*(double *)(base + align * offset[2]));
        t4 = (__vector float)vec_splats(*(double *)(base + align * offset[3]));
#endif
        t5  = vec_splats( *(base + align * offset[0] + 2) );
        t6  = vec_splats( *(base + align * offset[1] + 2) );
        t7  = vec_splats( *(base + align * offset[2] + 2) );
        t8  = vec_splats( *(base + align * offset[3] + 2) );

        t1  = vec_mergeh(t1, t2);
        t3  = vec_mergeh(t3, t4);
        v0  = vec_perm(t1, t3, (__vector unsigned char)perm_lo2hi);
        v1  = vec_perm(t3, t1, (__vector unsigned char)perm_hi2lo);
        t5  = vec_mergeh(t5, t6);
        t7  = vec_mergeh(t7, t8);
        v2  = vec_perm(t5, t7, (__vector unsigned char)perm_lo2hi);
    }
}

/* gcc-4.9 does not detect that vec_extract() uses its argument */
template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_f_ibm_vsx(float *                        base,
                                            const gmx_int32_t              offset[],
                                            gmx_simd_float_t               v0,
                                            gmx_simd_float_t               v1,
                                            gmx_simd_float_t gmx_unused    v2)
{
    gmx_simd_float_t t1, t2;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
#endif

    t1   = vec_mergeh(v0, v1);
    t2   = vec_mergel(v0, v1);
    *(double *)( base + align * offset[0] ) = vec_extract((__vector double)t1, 0);
    base[align*offset[0] + 2]               = vec_extract(v2, 0);
    *(double *)( base + align * offset[1] ) = vec_extract((__vector double)t1, 1);
    base[align*offset[1] + 2]               = vec_extract(v2, 1);
    *(double *)( base + align * offset[2] ) = vec_extract((__vector double)t2, 0);
    base[align*offset[2] + 2]               = vec_extract(v2, 2);
    *(double *)( base + align * offset[3] ) = vec_extract((__vector double)t2, 1);
    base[align*offset[3] + 2]               = vec_extract(v2, 3);
}

/* gcc-4.9 does not detect that vec_extract() uses its argument */
template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_f_ibm_vsx(float *                       base,
                                           const gmx_int32_t             offset[],
                                           gmx_simd_float_t              v0,
                                           gmx_simd_float_t              v1,
                                           gmx_simd_float_t gmx_unused   v2)
{
    const __vector unsigned char perm_lo2hi = { 0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23 };
    gmx_simd_float_t             t0, t1, m0, m1, m2, m3;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
#endif

    t0          = vec_mergeh(v0, v1);
    t1          = vec_mergel(v0, v1);

    m0 = (__vector float)vec_splats(*(double *)(base + align * offset[0]));
    m1 = (__vector float)vec_splats(*(double *)(base + align * offset[1]));
    m2 = (__vector float)vec_splats(*(double *)(base + align * offset[2]));
    m3 = (__vector float)vec_splats(*(double *)(base + align * offset[3]));

    m0 = vec_perm(m0, m1, (__vector unsigned char)perm_lo2hi);
    m2 = vec_perm(m2, m3, (__vector unsigned char)perm_lo2hi);
    m0 = vec_add(m0, t0);
    m2 = vec_add(m2, t1);

    *(double *)(base + align * offset[0])  = vec_extract( (__vector double)m0, 0);
    base[align * offset[0] + 2]           += vec_extract( v2, 0);
    *(double *)(base + align * offset[1])  = vec_extract( (__vector double)m0, 1);
    base[align * offset[1] + 2]           += vec_extract( v2, 1);
    *(double *)(base + align * offset[2])  = vec_extract( (__vector double)m2, 0);
    base[align * offset[2] + 2]           += vec_extract( v2, 2);
    *(double *)(base + align * offset[3])  = vec_extract( (__vector double)m2, 1);
    base[align * offset[3] + 2]           += vec_extract( v2, 3);
}

/* gcc-4.9 does not detect that vec_extract() uses its argument */
template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_f_ibm_vsx(float *                       base,
                                           const gmx_int32_t             offset[],
                                           gmx_simd_float_t              v0,
                                           gmx_simd_float_t              v1,
                                           gmx_simd_float_t gmx_unused   v2)
{
    const __vector unsigned char perm_lo2hi = { 0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23 };
    gmx_simd_float_t             t0, t1, m0, m1, m2, m3;

#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)offset % 16 == 0);
#endif

    t0          = vec_mergeh(v0, v1);
    t1          = vec_mergel(v0, v1);

    m0 = (__vector float)vec_splats(*(double *)(base + align * offset[0]));
    m1 = (__vector float)vec_splats(*(double *)(base + align * offset[1]));
    m2 = (__vector float)vec_splats(*(double *)(base + align * offset[2]));
    m3 = (__vector float)vec_splats(*(double *)(base + align * offset[3]));

    m0 = vec_perm(m0, m1, (__vector unsigned char)perm_lo2hi);
    m2 = vec_perm(m2, m3, (__vector unsigned char)perm_lo2hi);
    m0 = vec_sub(m0, t0);
    m2 = vec_sub(m2, t1);

    *(double *)(base + align * offset[0])  = vec_extract( (__vector double)m0, 0);
    base[align * offset[0] + 2]           -= vec_extract( v2, 0);
    *(double *)(base + align * offset[1])  = vec_extract( (__vector double)m0, 1);
    base[align * offset[1] + 2]           -= vec_extract( v2, 1);
    *(double *)(base + align * offset[2])  = vec_extract( (__vector double)m2, 0);
    base[align * offset[2] + 2]           -= vec_extract( v2, 2);
    *(double *)(base + align * offset[3])  = vec_extract( (__vector double)m2, 1);
    base[align * offset[3] + 2]           -= vec_extract( v2, 3);
}

static gmx_inline void
gmx_simd_expand_scalars_to_triplets_f_ibm_vsx(gmx_simd_float_t    scalar,
                                              gmx_simd_float_t   &triplets0,
                                              gmx_simd_float_t   &triplets1,
                                              gmx_simd_float_t   &triplets2)
{
    /* These permutes will be translated to immediate permutes (xxpermdi)
     * since they operate on doublewords, which will be faster than loading
     * the constants required for fully flexible permutes.
     * (although the real reason was that the latter was buggy on xlc-13.1).
     */
    __vector unsigned char perm0 = { 0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23 };
    __vector unsigned char perm1 = { 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };
    __vector unsigned char perm2 = { 8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31 };
    __vector float         t0, t1;

    t0        = vec_mergeh(scalar, scalar);
    t1        = vec_mergel(scalar, scalar);
    triplets0 = vec_perm(t0, scalar, perm0);
    triplets1 = vec_perm(t0, t1, perm1);
    triplets2 = vec_perm(scalar, t1, perm2);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_ibm_vsx(const float *       base,
                                                   gmx_simd_fint32_t   offset,
                                                   gmx_simd_float_t   &v0,
                                                   gmx_simd_float_t   &v1,
                                                   gmx_simd_float_t   &v2,
                                                   gmx_simd_float_t   &v3)
{
#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)base % 16 == 0);
    assert(align % 4 == 0);
#endif

    /* Use optimized bit-shift multiply for the most common alignments */
    if (align == 4)
    {
        offset = gmx_simd_slli_fi(offset, 2);
    }
    else if (align == 8)
    {
        offset = gmx_simd_slli_fi(offset, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset = gmx_simd_add_fi(offset, gmx_simd_slli_fi(offset, 1));
        offset = gmx_simd_slli_fi(offset, 2);
    }
    else if (align == 16)
    {
        offset = gmx_simd_slli_fi(offset, 4);
    }

    if (align == 4 || align == 8 || align == 12 || align == 16)
    {
        v0  = gmx_simd_load_f(base + gmx_simd_extract_fi(offset, 0) );
        v1  = gmx_simd_load_f(base + gmx_simd_extract_fi(offset, 1) );
        v2  = gmx_simd_load_f(base + gmx_simd_extract_fi(offset, 2) );
        v3  = gmx_simd_load_f(base + gmx_simd_extract_fi(offset, 3) );
    }
    else
    {
        v0  = gmx_simd_load_f(base + align * gmx_simd_extract_fi(offset, 0) );
        v1  = gmx_simd_load_f(base + align * gmx_simd_extract_fi(offset, 1) );
        v2  = gmx_simd_load_f(base + align * gmx_simd_extract_fi(offset, 2) );
        v3  = gmx_simd_load_f(base + align * gmx_simd_extract_fi(offset, 3) );
    }
    GMX_VSX_TRANSPOSE4(v0, v1, v2, v3);
}


template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_f_ibm_vsx(const float *       base,
                                                    gmx_simd_fint32_t   offset,
                                                    gmx_simd_float_t   &v0,
                                                    gmx_simd_float_t   &v1)
{
    gmx_simd_float_t t1, t2;

    /* Use optimized bit-shift multiply for the most common alignments */
    if (align == 2)
    {
        offset = gmx_simd_slli_fi(offset, 1);
    }
    else if (align == 4)
    {
        offset = gmx_simd_slli_fi(offset, 2);
    }
    else if (align == 6)
    {
        /* multiply by 3, then by 2 */
        offset = gmx_simd_add_fi(offset, gmx_simd_slli_fi(offset, 1));
        offset = gmx_simd_slli_fi(offset, 1);
    }
    else if (align == 8)
    {
        offset = gmx_simd_slli_fi(offset, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset = gmx_simd_add_fi(offset, gmx_simd_slli_fi(offset, 1));
        offset = gmx_simd_slli_fi(offset, 2);
    }
    else if (align == 16)
    {
        offset = gmx_simd_slli_fi(offset, 4);
    }

    if (align == 2 || align == 4 || align == 6 ||
        align == 8 || align == 12 || align == 16)
    {
        v0 = (__vector float)vec_splats(*(double *)(base + gmx_simd_extract_fi(offset, 0)));
        v1 = (__vector float)vec_splats(*(double *)(base + gmx_simd_extract_fi(offset, 1)));
        t1 = (__vector float)vec_splats(*(double *)(base + gmx_simd_extract_fi(offset, 2)));
        t2 = (__vector float)vec_splats(*(double *)(base + gmx_simd_extract_fi(offset, 3)));
    }
    else
    {
        v0 = (__vector float)vec_splats(*(double *)(base + align * gmx_simd_extract_fi(offset, 0)));
        v1 = (__vector float)vec_splats(*(double *)(base + align * gmx_simd_extract_fi(offset, 1)));
        t1 = (__vector float)vec_splats(*(double *)(base + align * gmx_simd_extract_fi(offset, 2)));
        t2 = (__vector float)vec_splats(*(double *)(base + align * gmx_simd_extract_fi(offset, 3)));
    }
    t1  = vec_mergeh(v0, t1);
    t2  = vec_mergeh(v1, t2);
    v0  = vec_mergeh(t1, t2);
    v1  = vec_mergel(t1, t2);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_ibm_vsx(const float *       base,
                                                   gmx_simd_fint32_t   offset,
                                                   gmx_simd_float_t   &v0,
                                                   gmx_simd_float_t   &v1)
{
#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)base % 8 == 0);
    assert(align % 2 == 0);
#endif

    gmx_simd_gather_loadu_bysimdint_transpose_f_ibm_vsx<align>(base, offset, v0, v1);
}



static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_f_ibm_vsx(float *           m,
                                            gmx_simd_float_t  v0,
                                            gmx_simd_float_t  v1,
                                            gmx_simd_float_t  v2,
                                            gmx_simd_float_t  v3)
{
#if !defined(__ibmxl__) && !defined(__xlC__)
    assert((size_t)m % 16 == 0);
#endif

    GMX_VSX_TRANSPOSE4(v0, v1, v2, v3);
    v0 = vec_add(v0, v1);
    v2 = vec_add(v2, v3);
    v0 = vec_add(v0, v2);
    v2 = vec_add(v0, gmx_simd_load_f(m));
    gmx_simd_store_f(m, v2);

    return gmx_simd_reduce_f_ibm_vsx(v0);
}
#endif /* __cplusplus */


#endif /* GMX_SIMD_IMPLEMENTATION_IBM_VSX_UTIL_FLOAT_H */
