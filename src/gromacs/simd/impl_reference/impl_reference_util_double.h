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

#ifndef GMX_SIMD_IMPL_REFERENCE_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_REFERENCE_UTIL_DOUBLE_H

/*! \libinternal \file
 *
 * \brief Reference impl., higher-level double prec. SIMD utility functions
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \ingroup module_simd
 */

#include <assert.h>
#include <math.h>

#include "impl_reference_common.h"
#include "impl_reference_simd_double.h"

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/* \name Higher-level SIMD utility functions, double precision.
 *
 * These include generic functions to work with triplets of data, typically
 * coordinates, and a few utility functions to load and update data in the
 * nonbonded kernels. These functions should be available on all implementations.
 *
 * This is part of the new C++ SIMD interface, so these functions are only
 * available when using C++. Since some Gromacs code reliying on the SIMD
 * module is still C (not C++), we have kept the C-style naming for now - this
 * will change once we are entirely C++.
 *
 * Note that we overload function names here for convenience.
 *
 * \{
 */

#ifdef __cplusplus
/*! \brief Load 4 consecutive double from each of GMX_SIMD_DOUBLE_WIDTH offsets,
 *         and transpose into 4 SIMD double variables.
 *
 * \tparam     align  Alignment of the storage, i.e. the distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of output components
 *                    the data is packed without padding.
 * \param      base   Pointer to the start of the memory area
 * \param      offset Array with offsets to the start of each data point.
 * \param[out] v0     1st component of data, base[align*offset[i]] for each i.
 * \param[out] v1     2nd component of data, base[align*offset[i] + 1] for each i.
 * \param[out] v2     3rd component of data, base[align*offset[i] + 2] for each i.
 * \param[out] v3     4th component of data, base[align*offset[i] + 2] for each i.
 *
 * The floating-point memory locations must be aligned, but only to the smaller
 * of four elements and the floating-point SIMD width.
 *
 * The offset memory must be aligned to the corresponding integer SIMD type,
 * i.e. GMX_SIMD_FINT32_WIDTH in single, or GMX_SIMD_DINT32_WIDTH in double.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 */
template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d(const double *        base,
                                 const gmx_int32_t     offset[],
                                 gmx_simd_double_t    &v0,
                                 gmx_simd_double_t    &v1,
                                 gmx_simd_double_t    &v2,
                                 gmx_simd_double_t    &v3)
{
    int i;
#ifndef NDEBUG
    int align_req = (GMX_SIMD_DOUBLE_WIDTH > 4) ? 4 : GMX_SIMD_DOUBLE_WIDTH;
#endif

    /* Make sure:
     * 1) the offset pointer fulfills hardware alignment demand
     * 2) the base pointer fulfills hardware alignment demand
     * 3) the specified alignment is a multiple of hardware alignment demand
     */
    assert((size_t)offset % (GMX_SIMD_DINT32_WIDTH*sizeof(gmx_int32_t)) == 0);
    assert((size_t)base % (align_req*sizeof(double)) == 0);
    assert(align % align_req == 0);

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        v0.r[i] = base[align * offset[i]];
        v1.r[i] = base[align * offset[i] + 1];
        v2.r[i] = base[align * offset[i] + 2];
        v3.r[i] = base[align * offset[i] + 3];
    }
}


/*! \brief Load 2 consecutive double from each of GMX_SIMD_DOUBLE_WIDTH offsets,
 *         and transpose into 2 SIMD double variables.
 *
 * \tparam     align  Alignment of the storage, i.e. the distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of output components
 *                    the data is packed without padding.
 * \param      base   Pointer to the start of the memory area
 * \param      offset Array with offsets to the start of each data point.
 * \param[out] v0     1st component of data, base[align*offset[i]] for each i.
 * \param[out] v1     2nd component of data, base[align*offset[i] + 1] for each i.
 *
 * The floating-point memory locations must be aligned, but only to the smaller
 * of two elements and the floating-point SIMD width.
 *
 * The offset memory must be aligned to the corresponding integer SIMD type,
 * i.e. GMX_SIMD_FINT32_WIDTH in single, or GMX_SIMD_DINT32_WIDTH in double.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 */
template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_d(const double *        base,
                                 const gmx_int32_t     offset[],
                                 gmx_simd_double_t    &v0,
                                 gmx_simd_double_t    &v1)
{
    int i;
#ifndef NDEBUG
    int align_req = (GMX_SIMD_DOUBLE_WIDTH > 2) ? 2 : GMX_SIMD_DOUBLE_WIDTH;
#endif

    /* Make sure:
     * 1) the offset pointer fulfills hardware alignment demand
     * 2) the base pointer fulfills hardware alignment demand
     * 3) the specified alignment is a multiple of hardware alignment demand
     */
    assert((size_t)offset % (GMX_SIMD_DINT32_WIDTH*sizeof(gmx_int32_t)) == 0);
    assert((size_t)base % (align_req*sizeof(double)) == 0);
    assert(align % align_req == 0);

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        v0.r[i] = base[align * offset[i]];
        v1.r[i] = base[align * offset[i] + 1];
    }
}


/*! \brief Best alignment to use for aligned pairs of double data.
 *
 * \copydetails gmx_simd_best_pair_alignment_f
 */
static const int gmx_simd_best_pair_alignment_d = 2;


/*! \brief Load 3 consecutive doubles from each of GMX_SIMD_DOUBLE_WIDTH offsets,
 *         and transpose into 3 SIMD double variables.
 *
 * \copydetails gmx_simd_gather_loadu_transpose_f
 */
template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_d(const double *        base,
                                  const gmx_int32_t     offset[],
                                  gmx_simd_double_t    &v0,
                                  gmx_simd_double_t    &v1,
                                  gmx_simd_double_t    &v2)
{
    int i;

    /* Make sure the offset pointer fulfills hardware alignment demand */
    assert((size_t)offset % (GMX_SIMD_DINT32_WIDTH*sizeof(gmx_int32_t)) == 0);

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        v0.r[i] = base[align * offset[i]];
        v1.r[i] = base[align * offset[i] + 1];
        v2.r[i] = base[align * offset[i] + 2];
    }
}


/*! \brief Transpose and store 3 SIMD doubles to 3 consecutive addresses at
 *         GMX_SIMD_DOUBLE_WIDTH offsets.
 *
 * \copydetails gmx_simd_transpose_scatter_storeu_f
 */
template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_d(double *              base,
                                    const gmx_int32_t     offset[],
                                    gmx_simd_double_t     v0,
                                    gmx_simd_double_t     v1,
                                    gmx_simd_double_t     v2)
{
    int i;

    /* Make sure the offset pointer fulfills hardware alignment demand */
    assert((size_t)offset % (GMX_SIMD_DINT32_WIDTH*sizeof(gmx_int32_t)) == 0);

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        base[align * offset[i]]     = v0.r[i];
        base[align * offset[i] + 1] = v1.r[i];
        base[align * offset[i] + 2] = v2.r[i];
    }
}


/*! \brief Transpose and add 3 SIMD doubles to 3 consecutive addresses at
 *         GMX_SIMD_DOUBLE_WIDTH offsets.
 *
 * \copydetails gmx_simd_transpose_scatter_incru_f
 */
template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_d(double *              base,
                                   const gmx_int32_t     offset[],
                                   gmx_simd_double_t     v0,
                                   gmx_simd_double_t     v1,
                                   gmx_simd_double_t     v2)
{
    int i;

    /* Make sure the offset pointer fulfills hardware alignment demand */
    assert((size_t)offset % (GMX_SIMD_DINT32_WIDTH*sizeof(gmx_int32_t)) == 0);

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        base[align * offset[i]]     += v0.r[i];
        base[align * offset[i] + 1] += v1.r[i];
        base[align * offset[i] + 2] += v2.r[i];
    }
}

/*! \brief Transpose and subtract 3 SIMD doubles to 3 consecutive addresses at
 *         GMX_SIMD_DOUBLE_WIDTH offsets.
 *
 * \copydetails gmx_simd_transpose_scatter_decru_f
 */
template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_d(double *              base,
                                   const gmx_int32_t     offset[],
                                   gmx_simd_double_t     v0,
                                   gmx_simd_double_t     v1,
                                   gmx_simd_double_t     v2)
{
    int i;

    /* Make sure the offset pointer fulfills hardware alignment demand */
    assert((size_t)offset % (GMX_SIMD_DINT32_WIDTH*sizeof(gmx_int32_t)) == 0);

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        base[align * offset[i]]     -= v0.r[i];
        base[align * offset[i] + 1] -= v1.r[i];
        base[align * offset[i] + 2] -= v2.r[i];
    }
}


/*! \brief Expand each element of double SIMD variable into three identical
 *         consecutive elements in three SIMD outputs.
 *
 * \copydetails gmx_simd_expand_scalars_to_triplets_f
 */
static gmx_inline void
gmx_simd_expand_scalars_to_triplets_d(gmx_simd_double_t    scalar,
                                      gmx_simd_double_t   &triplets0,
                                      gmx_simd_double_t   &triplets1,
                                      gmx_simd_double_t   &triplets2)
{
    int i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        triplets0.r[i] = scalar.r[i / 3];
        triplets1.r[i] = scalar.r[(i + GMX_SIMD_DOUBLE_WIDTH) / 3];
        triplets2.r[i] = scalar.r[(i + 2 * GMX_SIMD_DOUBLE_WIDTH) / 3];
    }
}


/*! \brief Load 4 consecutive doubles from each of GMX_SIMD_DOUBLE_WIDTH offsets
 *         specified by a SIMD integer, transpose into 4 SIMD double variables.
 *
 * \tparam     align  Alignment of the storage, i.e. the distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of output components
 *                    the data is packed without padding. This must be a
 *                    multiple of the alignment to keep all data aligned.
 * \param      base   Aligned pointer to the start of the memory.
 * \param      offset SIMD integer type with offsets to the start of each triplet.
 * \param[out] v0     First component, base[align*offset[i]] for each i.
 * \param[out] v1     Second component, base[align*offset[i] + 1] for each i.
 * \param[out] v2     Third component, base[align*offset[i] + 2] for each i.
 * \param[out] v3     Fourth component, base[align*offset[i] + 3] for each i.
 *
 * The memory locations must be aligned, but only to four elements even if the
 * SIMD implementation is even wider.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This is a special routine primarily intended for loading Gromacs
 *       table data as efficiently as possible - this is the reason for using
 *       a SIMD offset index, since the result of the  real-to-integer conversion
 *       is present in a SIMD register just before calling this routine.
 */
template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d(const double *       base,
                                           gmx_simd_dint32_t    offset,
                                           gmx_simd_double_t   &v0,
                                           gmx_simd_double_t   &v1,
                                           gmx_simd_double_t   &v2,
                                           gmx_simd_double_t   &v3)
{
    int i;
#ifndef NDEBUG
    int align_req = (GMX_SIMD_DOUBLE_WIDTH > 4) ? 4 : GMX_SIMD_DOUBLE_WIDTH;
#endif

    /* Make sure:
     * 1) the base pointer fulfills hardware alignment demand
     * 2) the specified alignment is a multiple of hardware alignment demand
     */
    assert((size_t)base % (align_req*sizeof(double)) == 0);
    assert(align % align_req == 0);

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        v0.r[i] = base[align * offset.i[i]];
        v1.r[i] = base[align * offset.i[i] + 1];
        v2.r[i] = base[align * offset.i[i] + 2];
        v3.r[i] = base[align * offset.i[i] + 3];
    }
}


/*! \brief Load 2 consecutive doubles from each of GMX_SIMD_DOUBLE_WIDTH offsets
 *         (unaligned) specified by SIMD integer, transpose into 2 SIMD doubles.
 *
 * \copydetails gmx_simd_gather_loadu_bysimdint_transpose_f
 */
template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_d(const double *       base,
                                            gmx_simd_dint32_t    offset,
                                            gmx_simd_double_t   &v0,
                                            gmx_simd_double_t   &v1)
{
    int i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        v0.r[i] = base[align * offset.i[i]];
        v1.r[i] = base[align * offset.i[i] + 1];
    }
}

/*! \brief Load 2 consecutive doubles from each of GMX_SIMD_DOUBLE_WIDTH offsets
 *         specified by a SIMD integer, transpose into 2 SIMD double variables.
 *
 * \tparam     align  Alignment of the storage, i.e. the distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of output components
 *                    the data is packed without padding. This must be a
 *                    multiple of the alignment to keep all data aligned.
 * \param      base   Aligned pointer to the start of the memory.
 * \param      offset SIMD integer type with offsets to the start of each triplet.
 * \param[out] v0     First component, base[align*offset[i]] for each i.
 * \param[out] v1     Second component, base[align*offset[i] + 1] for each i.
 *
 * The memory locations must be aligned, but only to four elements even if the
 * SIMD implementation is even wider.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This is a special routine primarily intended for loading Gromacs
 *       table data as efficiently as possible - this is the reason for using
 *       a SIMD offset index, since the result of the  real-to-integer conversion
 *       is present in a SIMD register just before calling this routine.
 */
template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_d(const double *       base,
                                           gmx_simd_dint32_t    offset,
                                           gmx_simd_double_t   &v0,
                                           gmx_simd_double_t   &v1)
{
#ifndef NDEBUG
    int align_req = (GMX_SIMD_DOUBLE_WIDTH > 2) ? 2 : GMX_SIMD_DOUBLE_WIDTH;
#endif

    /* Make sure:
     * 1) the base pointer fulfills hardware alignment demand
     * 2) the specified alignment is a multiple of hardware alignment demand
     */
    assert((size_t)base % (align_req*sizeof(double)) == 0);
    assert(align % align_req == 0);

    gmx_simd_gather_loadu_bysimdint_transpose_d<align>(base, offset, v0, v1);
}


/*! \brief Reduce each of four SIMD doubles, add those values to four consecutive
 *         doubles in memory, return sum.
 *
 * \copydetails gmx_simd_reduce_incr_4_return_sum_f
 */
static gmx_inline double
gmx_simd_reduce_incr_4_return_sum_d(double *           m,
                                    gmx_simd_double_t  v0,
                                    gmx_simd_double_t  v1,
                                    gmx_simd_double_t  v2,
                                    gmx_simd_double_t  v3)
{
    /* Note that the 4 here corresponds to the 4 m-elements, not any SIMD width */
    double sum[4];
#ifndef NDEBUG
    int    align_req = (GMX_SIMD_DOUBLE_WIDTH > 4) ? 4 : GMX_SIMD_DOUBLE_WIDTH;
#endif

    /* Make sure the memory pointer is aligned */
    assert((size_t)m % (align_req*sizeof(double)) == 0);

    sum[0] = gmx_simd_reduce_d(v0);
    sum[1] = gmx_simd_reduce_d(v1);
    sum[2] = gmx_simd_reduce_d(v2);
    sum[3] = gmx_simd_reduce_d(v3);

    m[0] += sum[0];
    m[1] += sum[1];
    m[2] += sum[2];
    m[3] += sum[3];

    return sum[0] + sum[1] + sum[2] + sum[3];
}
#endif /* __cplusplus */


/*! \}
 *
 * \name Higher-level SIMD utilities accessing partial (half-width) SIMD doubles.
 *
 * See the single-precision versions for documentation. Since double precision
 * is typically half the width of single, this double version is likely only
 * useful with 512-bit and larger implementations.
 *
 * \{
 */

#ifdef __cplusplus
/*! \brief Load low & high parts of SIMD float from different locations.
 *
 * \copydetails gmx_simd_load_dual_hsimd_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_load_dual_hsimd_d(const double *  m0,
                           const double *  m1)
{
    gmx_simd_double_t a;
    int               i;

    /* Make sure the memory pointers are aligned */
    assert((size_t)m0 % (GMX_SIMD_DOUBLE_WIDTH/2*sizeof(double)) == 0);
    assert((size_t)m1 % (GMX_SIMD_DOUBLE_WIDTH/2*sizeof(double)) == 0);

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        a.r[i]                           = m0[i];
        a.r[GMX_SIMD_DOUBLE_WIDTH/2 + i] = m1[i];
    }
    return a;
}

/*! \brief Load half-SIMD-width float data, spread to both halves.
 *
 * \copydetails gmx_simd_loaddup_hsimd_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_loaddup_hsimd_d(const double *  m)
{
    gmx_simd_double_t a;
    int               i;

    /* Make sure the memory pointer is aligned */
    assert((size_t)m % (GMX_SIMD_DOUBLE_WIDTH/2*sizeof(double)) == 0);

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        a.r[i]                           = m[i];
        a.r[GMX_SIMD_DOUBLE_WIDTH/2 + i] = a.r[i];
    }
    return a;
}

/*! \brief Load two floats, spread 1st in low half, 2nd in high half.
 *
 * \copydetails gmx_simd_load1_dual_hsimd_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_load1_dual_hsimd_d(const double *  m)
{
    gmx_simd_double_t a;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        a.r[i]                           = m[0];
        a.r[GMX_SIMD_DOUBLE_WIDTH/2 + i] = m[1];
    }
    return a;
}


/*! \brief Store low & high parts of SIMD float to different locations.
 *
 * \copydetails gmx_simd_store_dual_hsimd_f
 */
static gmx_inline void
gmx_simd_store_dual_hsimd_d(double *           m0,
                            double *           m1,
                            gmx_simd_double_t  a)
{
    int i;

    /* Make sure the memory pointers are aligned */
    assert((size_t)m0 % (GMX_SIMD_DOUBLE_WIDTH/2*sizeof(double)) == 0);
    assert((size_t)m1 % (GMX_SIMD_DOUBLE_WIDTH/2*sizeof(double)) == 0);

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        m0[i] = a.r[i];
        m1[i] = a.r[GMX_SIMD_DOUBLE_WIDTH/2 + i];
    }
}

/*! \brief Add the two halves of a SIMD double, subtract the sum from
 *         half-SIMD-width consecutive doubles in memory.
 *
 * \copydetails gmx_simd_decr_hsimd_f
 */
static gmx_inline void
gmx_simd_decr_hsimd_d(double *           m,
                      gmx_simd_double_t  a)
{
    int i;

    /* Make sure the memory pointer is aligned */
    assert((size_t)m % (GMX_SIMD_DOUBLE_WIDTH/2*sizeof(double)) == 0);

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        m[i] -= a.r[i] + a.r[GMX_SIMD_DOUBLE_WIDTH/2 + i];
    }
}


/*! \brief Load 2 consecutive doubles from each of GMX_SIMD_DOUBLE_WIDTH/2 offsets,
 *         transpose into SIMD double (low half from base0, high from base1).
 *
 * \copydetails gmx_simd_gather_load_transpose_hsimd_f
 */
template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_hsimd_d(const double *       base0,
                                       const double *       base1,
                                       gmx_int32_t          offset[],
                                       gmx_simd_double_t   &v0,
                                       gmx_simd_double_t   &v1)
{
    int i;
#ifndef NDEBUG
    int align_req = (GMX_SIMD_DOUBLE_WIDTH > 2) ? 2 : GMX_SIMD_DOUBLE_WIDTH;
#endif

    /* Make sure:
     * 1) the offset pointer fulfills hardware alignment demand
     * 2) Both base pointers fulfill hardware alignment demand
     * 3) the specified alignment is a multiple of hardware alignment demand
     */
    assert((size_t)offset % (GMX_SIMD_DINT32_WIDTH/2*sizeof(gmx_int32_t)) == 0);
    assert((size_t)base0 % (align_req*sizeof(double)) == 0);
    assert((size_t)base1 % (align_req*sizeof(double)) == 0);
    assert(align % align_req == 0);

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        v0.r[i] = base0[align * offset[i]];
        v1.r[i] = base0[align * offset[i] + 1];
        v0.r[GMX_SIMD_DOUBLE_WIDTH/2 + i] = base1[align * offset[i]];
        v1.r[GMX_SIMD_DOUBLE_WIDTH/2 + i] = base1[align * offset[i] + 1];
    }
}


/*! \brief Reduce the 4 half-SIMD-with doubles in 2 SIMD variables (sum halves),
 *         increment four consecutive doubles in memory, return sum.
 *
 * \copydetails gmx_simd_reduce_incr_4_return_sum_hsimd_f
 */
static gmx_inline double
gmx_simd_reduce_incr_4_return_sum_hsimd_d(double *           m,
                                          gmx_simd_double_t  v0,
                                          gmx_simd_double_t  v1)
{
    /* The 4 here corresponds to the 4 elements, not any SIMD width */
    double sum[4];
    int    i;
#ifndef NDEBUG
    int    align_req = (GMX_SIMD_DOUBLE_WIDTH > 4) ? 4 : GMX_SIMD_DOUBLE_WIDTH;
#endif

    /* Make sure the memory pointer is aligned */
    assert((size_t)m % (align_req*sizeof(double)) == 0);

    for (i = 0; i < 4; i++)
    {
        sum[i] = 0;
    }

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        sum[0] += v0.r[i];
        sum[1] += v0.r[GMX_SIMD_DOUBLE_WIDTH/2 + i];
        sum[2] += v1.r[i];
        sum[3] += v1.r[GMX_SIMD_DOUBLE_WIDTH/2 + i];
    }

    m[0] += sum[0];
    m[1] += sum[1];
    m[2] += sum[2];
    m[3] += sum[3];

    return sum[0] + sum[1] + sum[2] + sum[3];
}
#endif /* __cplusplus */

/*! \} */

/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_IMPL_REFERENCE_UTIL_DOUBLE_H */
