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

#ifndef GMX_SIMD_IMPL_REFERENCE_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_REFERENCE_UTIL_FLOAT_H

/*! \libinternal \file
 *
 * \brief Reference impl., higher-level single prec. SIMD utility functions
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \ingroup module_simd
 */

/* Avoid adding dependencies on the rest of GROMACS here (e.g. gmxassert.h)
 * since we want to be able run the low-level SIMD implementations independently
 * in simulators for new hardware.
 */


#include <cassert>
#include <cstddef>
#include <cstdint>

#include <algorithm>

#include "impl_reference_definitions.h"
#include "impl_reference_simd_float.h"

namespace gmx
{

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/* \name Higher-level SIMD utility functions, single precision.
 *
 * These include generic functions to work with triplets of data, typically
 * coordinates, and a few utility functions to load and update data in the
 * nonbonded kernels.
 * These functions should be available on all implementations, although
 * some wide SIMD implementations (width>=8) also provide special optional
 * versions to work with half or quarter registers to improve the performance
 * in the nonbonded kernels.
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

/*! \brief Load 4 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH offsets,
 *         and transpose into 4 SIMD float variables.
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
static inline void
simdGatherLoadTransposeF(const float  *        base,
                         const std::int32_t    offset[],
                         SimdFloat *           v0,
                         SimdFloat *           v1,
                         SimdFloat *           v2,
                         SimdFloat *           v3)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH*sizeof(std::int32_t)) == 0);
    // Base pointer must be aligned to the smaller of 4 elements and float SIMD width
    assert(std::size_t(base) % (std::min(GMX_SIMD_FLOAT_WIDTH, 4)*sizeof(float)) == 0);
    // align parameter must also be a multiple of the above alignment requirement
    assert(align % std::min(GMX_SIMD_FLOAT_WIDTH, 4) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        v0->r[i] = base[align * offset[i]];
        v1->r[i] = base[align * offset[i] + 1];
        v2->r[i] = base[align * offset[i] + 2];
        v3->r[i] = base[align * offset[i] + 3];
    }
}

/*! \brief Load 2 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH offsets,
 *         and transpose into 2 SIMD float variables.
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
 * To achieve the best possible performance, you should store your data with
 * alignment \ref c_simdBestPairAlignmentF in single,
 * \ref c_simdBestPairAlignmentD in double, or use
 * \ref c_simdBestPairAlignment for default Gromacs precision.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 */
template <int align>
static inline void
simdGatherLoadTransposeF(const float  *        base,
                         const std::int32_t    offset[],
                         SimdFloat *           v0,
                         SimdFloat *           v1)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH*sizeof(std::int32_t)) == 0);
    // Base pointer must be aligned to the smaller of 2 elements and float SIMD width
    assert(std::size_t(base) % (std::min(GMX_SIMD_FLOAT_WIDTH, 2)*sizeof(float)) == 0);
    // align parameter must also be a multiple of the above alignment requirement
    assert(align % std::min(GMX_SIMD_FLOAT_WIDTH, 2) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        v0->r[i] = base[align * offset[i]];
        v1->r[i] = base[align * offset[i] + 1];
    }
}


/*! \brief Best alignment to use for aligned pairs of float data.
 *
 *  The routines to load and transpose data will work with a wide range of
 *  alignments, but some might be faster than others, depending on the load
 *  instructions available in the hardware. This specifies the best
 *  alignment for each implementation when working with pairs of data.
 *
 *  To allow each architecture to use the most optimal form, we use a constant
 *  that code outside the SIMD module should use to store things properly. It
 *  must be at least 2. For example, a value of 2 means the two parameters A & B
 *  are stored as [A0 B0 A1 B1] while align-4 means [A0 B0 - - A1 B1 - -].
 *
 *  This alignment depends on the efficiency of partial-register load/store
 *  operations, and will depend on the architecture.
 */
static const int c_simdBestPairAlignmentF = 2;


/*! \brief Load 3 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH offsets,
 *         and transpose into 3 SIMD float variables.
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
 *
 * The floating-point memory does not have to be aligned. The offset memory
 * must be aligned to the corresponding integer SIMD type, i.e.
 * GMX_SIMD_FINT32_WIDTH in single, or GMX_SIMD_DINT32_WIDTH in double precision.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This routine uses a normal array for the offsets, since we typically
 *       load this data from memory. On the architectures we have tested this
 *       is faster even when a SIMD integer datatype is present.
 * \note If your alignment is known to be always 4, you will likely lose
 *       performance by using this function. In that case, a better alternative
 *       is to make sure your memory (base pointer) is aligned,
 *       and use the routines for loading four values instead, such as
 *       simdGatherLoadTranspose for real-value data. The extra
 *       dummy variable is optimized away by all compilers we have tried.
 */
template <int align>
static inline void
simdGatherLoadUTransposeF(const float  *        base,
                          const std::int32_t    offset[],
                          SimdFloat *           v0,
                          SimdFloat *           v1,
                          SimdFloat *           v2)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH*sizeof(std::int32_t)) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        v0->r[i] = base[align * offset[i]];
        v1->r[i] = base[align * offset[i] + 1];
        v2->r[i] = base[align * offset[i] + 2];
    }
}


/*! \brief Transpose and store 3 SIMD floats to 3 consecutive addresses at
 *         GMX_SIMD_FLOAT_WIDTH offsets.
 *
 * \tparam     align  Alignment of the storage, i.e. the distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of output components
 *                    the data is packed without padding.
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Aligned array with offsets to the start of each triplet.
 * \param      v0     1st component of triplets, written to base[align*offset[i]].
 * \param      v1     2nd component of triplets, written to base[align*offset[i] + 1].
 * \param      v2     3rd component of triplets, written to base[align*offset[i] + 2].
 *
 * The floating-point memory does not have to be aligned. The offset memory
 * must be aligned to the corresponding integer SIMD type, i.e.
 * GMX_SIMD_FINT32_WIDTH in single, or GMX_SIMD_DINT32_WIDTH in double precision.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This routine uses a normal array for the offsets, since we typically
 *       load the data from memory. On the architectures we have tested this
 *       is faster even when a SIMD integer datatype is present.
 */
template <int align>
static inline void
simdTransposeScatterStoreUF(float  *              base,
                            const std::int32_t    offset[],
                            SimdFloat             v0,
                            SimdFloat             v1,
                            SimdFloat             v2)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH*sizeof(std::int32_t)) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        base[align * offset[i]]     = v0.r[i];
        base[align * offset[i] + 1] = v1.r[i];
        base[align * offset[i] + 2] = v2.r[i];
    }
}


/*! \brief Transpose and add 3 SIMD floats to 3 consecutive addresses at
 *         GMX_SIMD_FLOAT_WIDTH offsets.
 *
 * \tparam     align  Alignment of the storage, i.e. the distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of output components
 *                    the data is packed without padding.
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Aligned array with offsets to the start of each triplet.
 * \param      v0     1st component of triplets, added to base[align*offset[i]].
 * \param      v1     2nd component of triplets, added to base[align*offset[i] + 1].
 * \param      v2     3rd component of triplets, added to base[align*offset[i] + 2].
 *
 * The floating-point memory does not have to be aligned. The offset memory
 * must be aligned to the corresponding integer SIMD type, i.e.
 * GMX_SIMD_FINT32_WIDTH in single, or GMX_SIMD_DINT32_WIDTH in double precision.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This routine uses a normal array for the offsets, since we typically
 *       load the data from memory. On the architectures we have tested this
 *       is faster even when a SIMD integer datatype is present.
 */
template <int align>
static inline void
simdTransposeScatterIncrUF(float  *              base,
                           const std::int32_t    offset[],
                           SimdFloat             v0,
                           SimdFloat             v1,
                           SimdFloat             v2)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH*sizeof(std::int32_t)) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        base[align * offset[i]]     += v0.r[i];
        base[align * offset[i] + 1] += v1.r[i];
        base[align * offset[i] + 2] += v2.r[i];
    }
}


/*! \brief Transpose and subtract 3 SIMD floats to 3 consecutive addresses at
 *         GMX_SIMD_FLOAT_WIDTH offsets.
 *
 * \tparam     align   Alignment of storage, i.e. distance between index points
 *                     measured in elements (not bytes). When this is identical
 *                     to the number of output components the data is packed.
 * \param[out] base    Pointer to start of memory.
 * \param      offset  Aligned array with offsets to the start of each triplet.
 * \param      v0      1st component, subtracted from base[align*offset[i]]
 * \param      v1      2nd component, subtracted from base[align*offset[i]+1]
 * \param      v2      3rd component, subtracted from base[align*offset[i]+2]
 *
 * The floating-point memory does not have to be aligned. The offset memory
 * must be aligned to the corresponding integer SIMD type, i.e.
 * GMX_SIMD_FINT32_WIDTH in single, or GMX_SIMD_DINT32_WIDTH in double precision.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This routine uses a normal array for the offsets, since we typically
 *       load the data from memory. On the architectures we have tested this
 *       is faster even when a SIMD integer datatype is present.
 */
template <int align>
static inline void
simdTransposeScatterDecrUF(float  *              base,
                           const std::int32_t    offset[],
                           SimdFloat             v0,
                           SimdFloat             v1,
                           SimdFloat             v2)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH*sizeof(std::int32_t)) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        base[align * offset[i]]     -= v0.r[i];
        base[align * offset[i] + 1] -= v1.r[i];
        base[align * offset[i] + 2] -= v2.r[i];
    }
}


/*! \brief Expand each element of float SIMD variable into three identical
 *         consecutive elements in three SIMD outputs.
 *
 * \param      scalar    Floating-point input, e.g. [s0 s1 s2 s3] if width=4.
 * \param[out] triplets0 First output, e.g. [s0 s0 s0 s1] if width=4.
 * \param[out] triplets1 Second output, e.g. [s1 s1 s2 s2] if width=4.
 * \param[out] triplets2 Third output, e.g. [s2 s3 s3 s3] if width=4.
 *
 * This routine is meant to use for things like scalar-vector multiplication,
 * where the vectors are stored in a merged format like [x0 y0 z0 x1 y1 z1 ...],
 * while the scalars are stored as [s0 s1 s2...], and the data cannot easily
 * be changed to SIMD-friendly layout.
 *
 * In this case, load 3 full-width SIMD variables from the vector array (This
 * will always correspond to GMX_SIMD_FLOAT_WIDTH/GMX_SIMD_DOUBLE_WIDTH
 * triplets), load a single full-width variable from the scalar array, and
 * call this routine to expand the data. You can then simply multiply the
 * first, second and third pair of SIMD variables, and store the three
 * results back into a suitable vector-format array.
 */
static inline void
simdExpandScalarsToTripletsF(SimdFloat    scalar,
                             SimdFloat *  triplets0,
                             SimdFloat *  triplets1,
                             SimdFloat *  triplets2)
{
    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        triplets0->r[i] = scalar.r[i / 3];
        triplets1->r[i] = scalar.r[(i + GMX_SIMD_FLOAT_WIDTH) / 3];
        triplets2->r[i] = scalar.r[(i + 2 * GMX_SIMD_FLOAT_WIDTH) / 3];
    }
}

/*! \brief Load 4 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH offsets
 *         specified by a SIMD integer, transpose into 4 SIMD float variables.
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
 * The floating-point memory locations must be aligned, but only to the smaller
 * of four elements and the floating-point SIMD width.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This is a special routine primarily intended for loading Gromacs
 *       table data as efficiently as possible - this is the reason for using
 *       a SIMD offset index, since the result of the  real-to-integer conversion
 *       is present in a SIMD register just before calling this routine.
 */
template <int align>
static inline void
simdGatherLoadBySimdIntTransposeF(const float *       base,
                                  SimdFInt32          offset,
                                  SimdFloat *         v0,
                                  SimdFloat *         v1,
                                  SimdFloat *         v2,
                                  SimdFloat *         v3)
{
    // Base pointer must be aligned to the smaller of 4 elements and float SIMD width
    assert(std::size_t(base) % (std::min(GMX_SIMD_FLOAT_WIDTH, 4)*sizeof(float)) == 0);
    // align parameter must also be a multiple of the above alignment requirement
    assert(align % std::min(GMX_SIMD_FLOAT_WIDTH, 4) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        v0->r[i] = base[align * offset.i[i]];
        v1->r[i] = base[align * offset.i[i] + 1];
        v2->r[i] = base[align * offset.i[i] + 2];
        v3->r[i] = base[align * offset.i[i] + 3];
    }
}


/*! \brief Load 2 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH offsets
 *         (unaligned) specified by SIMD integer, transpose into 2 SIMD floats.
 *
 * \tparam     align  Alignment of the storage, i.e. the distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of output components
 *                    the data is packed without padding.
 * \param      base   Pointer to the start of the memory.
 * \param      offset SIMD integer type with offsets to the start of each triplet.
 * \param[out] v0     First component, base[align*offset[i]] for each i.
 * \param[out] v1     Second component, base[align*offset[i] + 1] for each i.
 *
 * This routine is only available when
 * GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE is defined.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This is a special routine primarily intended for loading Gromacs
 *       table data as efficiently as possible - this is the reason for using
 *       a SIMD offset index, since the result of the  real-to-integer conversion
 *       is present in a SIMD register just before calling this routine.
 */
template <int align>
static inline void
simdGatherLoadUBySimdIntTransposeF(const float *       base,
                                   SimdFInt32          offset,
                                   SimdFloat *         v0,
                                   SimdFloat *         v1)
{
    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        v0->r[i] = base[align * offset.i[i]];
        v1->r[i] = base[align * offset.i[i] + 1];
    }
}


/*! \brief Load 2 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH offsets
 *         specified by a SIMD integer, transpose into 2 SIMD float variables.
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
 * The floating-point memory locations must be aligned, but only to the smaller
 * of two elements and the floating-point SIMD width.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This is a special routine primarily intended for loading Gromacs
 *       table data as efficiently as possible - this is the reason for using
 *       a SIMD offset index, since the result of the  real-to-integer conversion
 *       is present in a SIMD register just before calling this routine.
 */
template <int align>
static inline void
simdGatherLoadBySimdIntTransposeF(const float *       base,
                                  SimdFInt32          offset,
                                  SimdFloat *         v0,
                                  SimdFloat *         v1)
{
    // Base pointer must be aligned to the smaller of 2 elements and float SIMD width
    assert(std::size_t(base) % (std::min(GMX_SIMD_FLOAT_WIDTH, 2)*sizeof(float)) == 0);
    // align parameter must also be a multiple of the above alignment requirement
    assert(align % std::min(GMX_SIMD_FLOAT_WIDTH, 2) == 0);

    simdGatherLoadUBySimdIntTransposeF<align>(base, offset, v0, v1);
}


/*! \brief Reduce each of four SIMD floats, add those values to four consecutive
 *         floats in memory, return sum.
 *
 * \param m   Pointer to memory where four floats should be incremented
 * \param v0  SIMD variable whose sum should be added to m[0]
 * \param v1  SIMD variable whose sum should be added to m[1]
 * \param v2  SIMD variable whose sum should be added to m[2]
 * \param v3  SIMD variable whose sum should be added to m[3]
 *
 * \return Sum of all elements in the four SIMD variables.
 *
 * The pointer m must be aligned to the smaller of four elements and the
 * floating-point SIMD width.
 *
 * \note This is a special routine intended for the Gromacs nonbonded kernels.
 * It is used in the epilogue of the outer loop, where the variables will
 * contain unrolled forces for one outer-loop-particle each, corresponding to
 * a single coordinate (i.e, say, four x-coordinate force variables). These
 * should be summed and added to the force array in memory. Since we always work
 * with contiguous SIMD-layout , we can use efficient aligned loads/stores.
 * When calculating the virial, we also need the total sum of all forces for
 * each coordinate. This is provided as the return value. For routines that
 * do not need these, this extra code will be optimized away completely if you
 * just ignore the return value (Checked with gcc-4.9.1 and clang-3.6 for AVX).
 */
static inline float
simdReduceIncr4ReturnSumF(float *           m,
                          SimdFloat         v0,
                          SimdFloat         v1,
                          SimdFloat         v2,
                          SimdFloat         v3)
{
    float sum[4]; // Note that the 4 here corresponds to the 4 m-elements, not any SIMD width

    // Make sure the memory pointer is aligned to the smaller of 4 elements and float SIMD width
    assert(std::size_t(m) % (std::min(GMX_SIMD_FLOAT_WIDTH, 4)*sizeof(float)) == 0);

    sum[0] = simdReduceF(v0);
    sum[1] = simdReduceF(v1);
    sum[2] = simdReduceF(v2);
    sum[3] = simdReduceF(v3);

    m[0] += sum[0];
    m[1] += sum[1];
    m[2] += sum[2];
    m[3] += sum[3];

    return sum[0] + sum[1] + sum[2] + sum[3];
}


/*! \}
 *
 * \name Higher-level SIMD utilities accessing partial (half-width) SIMD floats.
 *
 * These functions are optional. The are only useful for SIMD implementation
 * where the width is 8 or larger, and where it would be inefficient
 * to process 4*8, 8*8, or more, interactions in parallel.
 *
 * Currently, only Intel provides very wide SIMD implementations, but these
 * also come with excellent support for loading, storing, accessing and
 * shuffling parts of the register in so-called 'lanes' of 4 bytes each.
 * We can use this to load separate parts into the low/high halves of the
 * register in the inner loop of the nonbonded kernel, which e.g. makes it
 * possible to process 4*4 nonbonded interactions as a pattern of 2*8. We
 * can also use implementations with width 16 or greater.
 *
 * To make this more generic, when \ref GMX_SIMD_HAVE_HSIMD_UTIL_REAL is
 * defined, the SIMD implementation provides seven special routines that:
 *
 * - Load the low/high parts of a SIMD variable from different pointers
 * - Load half the SIMD width from one pointer, and duplicate in low/high parts
 * - Load two reals, put 1st one in all low elements, and 2nd in all high ones.
 * - Store the low/high parts of a SIMD variable to different pointers
 * - Subtract both SIMD halves from a single half-SIMD-width memory location.
 * - Load aligned pairs (LJ parameters) from two base pointers, with a common
 *   offset list, and put these in the low/high SIMD halves.
 * - Reduce each half of two SIMD registers (i.e., 4 parts in total), increment
 *   four adjacent memory positions, and return the total sum.
 *
 * Remember: this is ONLY used when the native SIMD width is large. You will
 * just waste time if you implement it for normal 16-byte SIMD architectures.
 *
 * This is part of the new C++ SIMD interface, so these functions are only
 * available when using C++. Since some Gromacs code reliying on the SIMD
 * module is still C (not C++), we have kept the C-style naming for now - this
 * will change once we are entirely C++.
 *
 * \{
 */

/*! \brief Load low & high parts of SIMD float from different locations.
 *
 * \param m0 Pointer to memory aligned to half SIMD width.
 * \param m1 Pointer to memory aligned to half SIMD width.
 *
 * \return SIMD variable with low part loaded from m0, high from m1.
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE for double.
 */
static inline SimdFloat
simdLoadDualHsimdF(const float *  m0,
                   const float *  m1)
{
    SimdFloat        a;

    // Make sure the memory pointers are aligned to half float SIMD width
    assert(std::size_t(m0) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);
    assert(std::size_t(m1) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        a.r[i]                          = m0[i];
        a.r[GMX_SIMD_FLOAT_WIDTH/2 + i] = m1[i];
    }
    return a;
}

/*! \brief Load half-SIMD-width float data, spread to both halves.
 *
 * \param m Pointer to memory aligned to half SIMD width.
 *
 * \return SIMD variable with both halves loaded from m..
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE for double.
 */
static inline SimdFloat
simdLoadDupHsimdF(const float *  m)
{
    SimdFloat        a;

    // Make sure the memory pointer is aligned
    assert(std::size_t(m) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        a.r[i]                          = m[i];
        a.r[GMX_SIMD_FLOAT_WIDTH/2 + i] = a.r[i];
    }
    return a;
}

/*! \brief Load two floats, spread 1st in low half, 2nd in high half.
 *
 * \param m Pointer to two adjacent floating-point values.
 *
 * \return SIMD variable where all elements in the low half have been set
 *         to m[0], and all elements in high half to m[1].
 *
 * \note This routine always loads two values and sets the halves separately.
 *       If you want to set all elements to the same value, simply use
 *       the standard \ref simdLoad1().
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE for double.
 */
static inline SimdFloat
simdLoad1DualHsimdF(const float *  m)
{
    SimdFloat        a;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        a.r[i]                          = m[0];
        a.r[GMX_SIMD_FLOAT_WIDTH/2 + i] = m[1];
    }
    return a;
}


/*! \brief Store low & high parts of SIMD float to different locations.
 *
 * \param m0 Pointer to memory aligned to half SIMD width.
 * \param m1 Pointer to memory aligned to half SIMD width.
 * \param a  SIMD variable. Low half should be stored to m0, high to m1.
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE for double.
 */
static inline void
simdStoreDualHsimdF(float *           m0,
                    float *           m1,
                    SimdFloat         a)
{
    // Make sure the memory pointers are aligned to half float SIMD width
    assert(std::size_t(m0) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);
    assert(std::size_t(m1) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        m0[i] = a.r[i];
        m1[i] = a.r[GMX_SIMD_FLOAT_WIDTH/2 + i];
    }
}

/*! \brief Add the two halves of a SIMD float, subtract the sum from
 *         half-SIMD-width consecutive floats in memory.
 *
 * \param m  half-width aligned memory, from which sum of the halves will be subtracted.
 * \param a  SIMD variable. Upper & lower halves will first be added.
 *
 * If the SIMD width is 8 and contains [a b c d e f g h], the
 * memory will be modified to [m[0]-(a+e) m[1]-(b+f) m[2]-(c+g) m[3]-(d+h)].
 *
 * The memory must be aligned to half SIMD width.
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE for double.
 */
static inline void
simdDecrHsimdF(float *           m,
               SimdFloat         a)
{
    // Make sure the memory pointer is aligned to half float SIMD width
    assert(std::size_t(m) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        m[i] -= a.r[i] + a.r[GMX_SIMD_FLOAT_WIDTH/2 + i];
    }
}

/*! \brief Load 2 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH/2 offsets,
 *         transpose into SIMD float (low half from base0, high from base1).
 *
 * \tparam     align  Alignment of the storage, i.e. the distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of output components
 *                    the data is packed without padding. This must be a
 *                    multiple of the alignment to keep all data aligned.
 * \param      base0  Pointer to base of first aligned memory
 * \param      base1  Pointer to base of second aligned memory
 * \param      offset Offset to the start of each pair
 * \param[out] v0     1st element in each pair, base0 in low and base1 in high half.
 * \param[out] v1     2nd element in each pair, base0 in low and base1 in high half.
 *
 * The offset array should be of half the SIMD width length, so it corresponds
 * to the half-SIMD-register operations. This also means it must be aligned
 * to half the integer SIMD width (i.e., GMX_SIMD_FINT32_WIDTH/2).
 *
 * The floating-point memory locations must be aligned, but only to the smaller
 * of two elements and the floating-point SIMD width.
 *
 * This routine is primarily designed to load nonbonded parameters in the
 * kernels. It is the equivalent of the full-width routine
 * simdGatherLoadTranspose(), but just
 * as the other hsimd routines it will pick half-SIMD-width data from base0
 * and put in the lower half, while the upper half comes from base1.
 *
 * For an example, assume the SIMD width is 8, align is 2, that
 * base0 is [A0 A1 B0 B1 C0 C1 D0 D1 ...], and base1 [E0 E1 F0 F1 G0 G1 H0 H1...].
 *
 * Then we will get v0 as [A0 B0 C0 D0 E0 F0 G0 H0] and v1 as [A1 B1 C1 D1 E1 F1 G1 H1].
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE for double.
 */
template <int align>
static inline void
simdGatherLoadTransposeHsimdF(const float  *       base0,
                              const float  *       base1,
                              std::int32_t         offset[],
                              SimdFloat *          v0,
                              SimdFloat *          v1)
{
    // Offset list must be aligned for half SIMD FINT32 width
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH/2*sizeof(std::int32_t)) == 0);
    // base pointers must be aligned to the smaller of 2 elements and float SIMD width
    assert(std::size_t(base0) % (std::min(GMX_SIMD_FLOAT_WIDTH, 2)*sizeof(float)) == 0);
    assert(std::size_t(base1) % (std::min(GMX_SIMD_FLOAT_WIDTH, 2)*sizeof(float)) == 0);
    // alignment parameter must be also be multiple of the above required alignment
    assert(align % std::min(GMX_SIMD_FLOAT_WIDTH, 2) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        v0->r[i] = base0[align * offset[i]];
        v1->r[i] = base0[align * offset[i] + 1];
        v0->r[GMX_SIMD_FLOAT_WIDTH/2 + i] = base1[align * offset[i]];
        v1->r[GMX_SIMD_FLOAT_WIDTH/2 + i] = base1[align * offset[i] + 1];
    }
}

/*! \brief Reduce the 4 half-SIMD-with floats in 2 SIMD variables (sum halves),
 *         increment four consecutive floats in memory, return sum.
 *
 * \param m    Pointer to memory where the four values should be incremented
 * \param v0   Variable whose half-SIMD sums should be added to m[0]/m[1], respectively.
 * \param v1   Variable whose half-SIMD sums should be added to m[2]/m[3], respectively.
 *
 * \return Sum of all elements in the four SIMD variables.
 *
 * The pointer m must be aligned, but only to the smaller
 * of four elements and the floating-point SIMD width.
 *
 * \note This is the half-SIMD-width version of
 *      simdReduceIncr4ReturnSum(). The only difference is that the
 *      four half-SIMD inputs needed are present in the low/high halves of the
 *      two SIMD arguments.
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE for double.
 */
static inline float
simdReduceIncr4ReturnSumHsimdF(float *            m,
                               SimdFloat          v0,
                               SimdFloat          v1)
{
    // The 4 here corresponds to the 4 elements in memory, not any SIMD width
    float sum[4] = { 0.0f, 0.0f, 0.0f, 0.0f };

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        sum[0] += v0.r[i];
        sum[1] += v0.r[GMX_SIMD_FLOAT_WIDTH/2 + i];
        sum[2] += v1.r[i];
        sum[3] += v1.r[GMX_SIMD_FLOAT_WIDTH/2 + i];
    }

    // Make sure the memory pointer is aligned to the smaller of 4 elements and float SIMD width
    assert(std::size_t(m) % (std::min(GMX_SIMD_FLOAT_WIDTH, 4)*sizeof(float)) == 0);

    m[0] += sum[0];
    m[1] += sum[1];
    m[2] += sum[2];
    m[3] += sum[3];

    return sum[0] + sum[1] + sum[2] + sum[3];
}

/*! \} */

/*! \} */
/*! \endcond */

}      // namespace gmx

#endif // GMX_SIMD_IMPL_REFERENCE_UTIL_FLOAT_H
