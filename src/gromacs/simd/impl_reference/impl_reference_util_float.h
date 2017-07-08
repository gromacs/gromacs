/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
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

#include "config.h"

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
 * \{
 */

/*! \brief Load 4 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH offsets,
 *         and transpose into 4 SIMD float variables.
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of SIMD variables
 *                    (i.e., 4 for this routine) the input data is packed without
 *                    padding in memory. See the SIMD parameters for exactly
 *                    what memory positions are loaded.
 * \param      base   Pointer to the start of the memory area
 * \param      offset Array with offsets to the start of each data point.
 * \param[out] v0     1st component of data, base[align*offset[i]] for each i.
 * \param[out] v1     2nd component of data, base[align*offset[i] + 1] for each i.
 * \param[out] v2     3rd component of data, base[align*offset[i] + 2] for each i.
 * \param[out] v3     4th component of data, base[align*offset[i] + 3] for each i.
 *
 * The floating-point memory locations must be aligned, but only to the smaller
 * of four elements and the floating-point SIMD width.
 *
 * The offset memory must be aligned to GMX_SIMD_DINT32_WIDTH.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 */
template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const float  *        base,
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

    for (std::size_t i = 0; i < v0->simdInternal_.size(); i++)
    {
        v0->simdInternal_[i] = base[align * offset[i]];
        v1->simdInternal_[i] = base[align * offset[i] + 1];
        v2->simdInternal_[i] = base[align * offset[i] + 2];
        v3->simdInternal_[i] = base[align * offset[i] + 3];
    }
}

/*! \brief Load 2 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH offsets,
 *         and transpose into 2 SIMD float variables.
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of SIMD variables
 *                    (i.e., 2 for this routine) the input data is packed without
 *                    padding in memory. See the SIMD parameters for exactly
 *                    what memory positions are loaded.
 * \param      base   Pointer to the start of the memory area
 * \param      offset Array with offsets to the start of each data point.
 * \param[out] v0     1st component of data, base[align*offset[i]] for each i.
 * \param[out] v1     2nd component of data, base[align*offset[i] + 1] for each i.
 *
 * The floating-point memory locations must be aligned, but only to the smaller
 * of two elements and the floating-point SIMD width.
 *
 * The offset memory must be aligned to GMX_SIMD_FINT32_WIDTH.
 *
 * To achieve the best possible performance, you should store your data with
 * alignment \ref c_simdBestPairAlignmentFloat in single, or
 * \ref c_simdBestPairAlignmentDouble in double.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 */
template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const float  *        base,
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

    for (std::size_t i = 0; i < v0->simdInternal_.size(); i++)
    {
        v0->simdInternal_[i] = base[align * offset[i]];
        v1->simdInternal_[i] = base[align * offset[i] + 1];
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
static const int c_simdBestPairAlignmentFloat = 2;


/*! \brief Load 3 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH offsets,
 *         and transpose into 3 SIMD float variables.
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of SIMD variables
 *                    (i.e., 3 for this routine) the input data is packed without
 *                    padding in memory. See the SIMD parameters for exactly
 *                    what memory positions are loaded.
 * \param      base   Pointer to the start of the memory area
 * \param      offset Array with offsets to the start of each data point.
 * \param[out] v0     1st component of data, base[align*offset[i]] for each i.
 * \param[out] v1     2nd component of data, base[align*offset[i] + 1] for each i.
 * \param[out] v2     3rd component of data, base[align*offset[i] + 2] for each i.
 *
 * This function can work with both aligned (better performance) and unaligned
 * memory. When the align parameter is not a power-of-two (align==3 would be normal
 * for packed atomic coordinates) the memory obviously cannot be aligned, and
 * we account for this.
 * However, in the case where align is a power-of-two, we assume the base pointer
 * also has the same alignment, which will enable many platforms to use faster
 * aligned memory load operations.
 * An easy way to think of this is that each triplet of data in memory must be
 * aligned to the align parameter you specify when it's a power-of-two.
 *
 * The offset memory must always be aligned to GMX_SIMD_FINT32_WIDTH, since this
 * enables us to use SIMD loads and gather operations on platforms that support it.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This routine uses a normal array for the offsets, since we typically
 *       load this data from memory. On the architectures we have tested this
 *       is faster even when a SIMD integer datatype is present.
 * \note To improve performance, this function might use full-SIMD-width
 *       unaligned loads. This means you need to ensure the memory is padded
 *       at the end, so we always can load GMX_SIMD_REAL_WIDTH elements
 *       starting at the last offset. If you use the Gromacs aligned memory
 *       allocation routines this will always be the case.
 */
template <int align>
static inline void gmx_simdcall
gatherLoadUTranspose(const float  *        base,
                     const std::int32_t    offset[],
                     SimdFloat *           v0,
                     SimdFloat *           v1,
                     SimdFloat *           v2)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH*sizeof(std::int32_t)) == 0);

    for (std::size_t i = 0; i < v0->simdInternal_.size(); i++)
    {
        v0->simdInternal_[i] = base[align * offset[i]];
        v1->simdInternal_[i] = base[align * offset[i] + 1];
        v2->simdInternal_[i] = base[align * offset[i] + 2];
    }
}


/*! \brief Transpose and store 3 SIMD floats to 3 consecutive addresses at
 *         GMX_SIMD_FLOAT_WIDTH offsets.
 *
 * \tparam     align  Alignment of the memory to which we write, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of SIMD variables
 *                    (i.e., 3 for this routine) the output data is packed without
 *                    padding in memory. See the SIMD parameters for exactly
 *                    what memory positions are written.
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Aligned array with offsets to the start of each triplet.
 * \param      v0     1st component of triplets, written to base[align*offset[i]].
 * \param      v1     2nd component of triplets, written to base[align*offset[i] + 1].
 * \param      v2     3rd component of triplets, written to base[align*offset[i] + 2].
 *
 * This function can work with both aligned (better performance) and unaligned
 * memory. When the align parameter is not a power-of-two (align==3 would be normal
 * for packed atomic coordinates) the memory obviously cannot be aligned, and
 * we account for this.
 * However, in the case where align is a power-of-two, we assume the base pointer
 * also has the same alignment, which will enable many platforms to use faster
 * aligned memory store operations.
 * An easy way to think of this is that each triplet of data in memory must be
 * aligned to the align parameter you specify when it's a power-of-two.
 *
 * The offset memory must always be aligned to GMX_SIMD_FINT32_WIDTH, since this
 * enables us to use SIMD loads and gather operations on platforms that support it.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This routine uses a normal array for the offsets, since we typically
 *       load the data from memory. On the architectures we have tested this
 *       is faster even when a SIMD integer datatype is present.
 */
template <int align>
static inline void gmx_simdcall
transposeScatterStoreU(float  *              base,
                       const std::int32_t    offset[],
                       SimdFloat             v0,
                       SimdFloat             v1,
                       SimdFloat             v2)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH*sizeof(std::int32_t)) == 0);

    for (std::size_t i = 0; i < v0.simdInternal_.size(); i++)
    {
        base[align * offset[i]]     = v0.simdInternal_[i];
        base[align * offset[i] + 1] = v1.simdInternal_[i];
        base[align * offset[i] + 2] = v2.simdInternal_[i];
    }
}


/*! \brief Transpose and add 3 SIMD floats to 3 consecutive addresses at
 *         GMX_SIMD_FLOAT_WIDTH offsets.
 *
 * \tparam     align  Alignment of the memory to which we write, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of SIMD variables
 *                    (i.e., 3 for this routine) the output data is packed without
 *                    padding in memory. See the SIMD parameters for exactly
 *                    what memory positions are incremented.
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Aligned array with offsets to the start of each triplet.
 * \param      v0     1st component of triplets, added to base[align*offset[i]].
 * \param      v1     2nd component of triplets, added to base[align*offset[i] + 1].
 * \param      v2     3rd component of triplets, added to base[align*offset[i] + 2].
 *
 * This function can work with both aligned (better performance) and unaligned
 * memory. When the align parameter is not a power-of-two (align==3 would be normal
 * for packed atomic coordinates) the memory obviously cannot be aligned, and
 * we account for this.
 * However, in the case where align is a power-of-two, we assume the base pointer
 * also has the same alignment, which will enable many platforms to use faster
 * aligned memory load/store operations.
 * An easy way to think of this is that each triplet of data in memory must be
 * aligned to the align parameter you specify when it's a power-of-two.
 *
 * The offset memory must always be aligned to GMX_SIMD_FINT32_WIDTH, since this
 * enables us to use SIMD loads and gather operations on platforms that support it.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This routine uses a normal array for the offsets, since we typically
 *       load the data from memory. On the architectures we have tested this
 *       is faster even when a SIMD integer datatype is present.
 * \note To improve performance, this function might use full-SIMD-width
 *       unaligned load/store, and add 0.0 to the extra elements.
 *       This means you need to ensure the memory is padded
 *       at the end, so we always can load GMX_SIMD_REAL_WIDTH elements
 *       starting at the last offset. If you use the Gromacs aligned memory
 *       allocation routines this will always be the case.
 */
template <int align>
static inline void gmx_simdcall
transposeScatterIncrU(float  *              base,
                      const std::int32_t    offset[],
                      SimdFloat             v0,
                      SimdFloat             v1,
                      SimdFloat             v2)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH*sizeof(std::int32_t)) == 0);

    for (std::size_t i = 0; i < v0.simdInternal_.size(); i++)
    {
        base[align * offset[i]]     += v0.simdInternal_[i];
        base[align * offset[i] + 1] += v1.simdInternal_[i];
        base[align * offset[i] + 2] += v2.simdInternal_[i];
    }
}


/*! \brief Transpose and subtract 3 SIMD floats to 3 consecutive addresses at
 *         GMX_SIMD_FLOAT_WIDTH offsets.
 *
 * \tparam     align  Alignment of the memory to which we write, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of SIMD variables
 *                    (i.e., 3 for this routine) the output data is packed without
 *                    padding in memory. See the SIMD parameters for exactly
 *                    what memory positions are decremented.
 * \param[out] base    Pointer to start of memory.
 * \param      offset  Aligned array with offsets to the start of each triplet.
 * \param      v0      1st component, subtracted from base[align*offset[i]]
 * \param      v1      2nd component, subtracted from base[align*offset[i]+1]
 * \param      v2      3rd component, subtracted from base[align*offset[i]+2]
 *
 * This function can work with both aligned (better performance) and unaligned
 * memory. When the align parameter is not a power-of-two (align==3 would be normal
 * for packed atomic coordinates) the memory obviously cannot be aligned, and
 * we account for this.
 * However, in the case where align is a power-of-two, we assume the base pointer
 * also has the same alignment, which will enable many platforms to use faster
 * aligned memory load/store operations.
 * An easy way to think of this is that each triplet of data in memory must be
 * aligned to the align parameter you specify when it's a power-of-two.
 *
 * The offset memory must always be aligned to GMX_SIMD_FINT32_WIDTH, since this
 * enables us to use SIMD loads and gather operations on platforms that support it.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This routine uses a normal array for the offsets, since we typically
 *       load the data from memory. On the architectures we have tested this
 *       is faster even when a SIMD integer datatype is present.
 * \note To improve performance, this function might use full-SIMD-width
 *       unaligned load/store, and subtract 0.0 from the extra elements.
 *       This means you need to ensure the memory is padded
 *       at the end, so we always can load GMX_SIMD_REAL_WIDTH elements
 *       starting at the last offset. If you use the Gromacs aligned memory
 *       allocation routines this will always be the case.
 */
template <int align>
static inline void gmx_simdcall
transposeScatterDecrU(float  *              base,
                      const std::int32_t    offset[],
                      SimdFloat             v0,
                      SimdFloat             v1,
                      SimdFloat             v2)
{
    // Offset list must be aligned for SIMD FINT32
    assert(std::size_t(offset) % (GMX_SIMD_FINT32_WIDTH*sizeof(std::int32_t)) == 0);

    for (std::size_t i = 0; i < v0.simdInternal_.size(); i++)
    {
        base[align * offset[i]]     -= v0.simdInternal_[i];
        base[align * offset[i] + 1] -= v1.simdInternal_[i];
        base[align * offset[i] + 2] -= v2.simdInternal_[i];
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
 * will always correspond to GMX_SIMD_FLOAT_WIDTH triplets),
 * load a single full-width variable from the scalar array, and
 * call this routine to expand the data. You can then simply multiply the
 * first, second and third pair of SIMD variables, and store the three
 * results back into a suitable vector-format array.
 */
static inline void gmx_simdcall
expandScalarsToTriplets(SimdFloat    scalar,
                        SimdFloat *  triplets0,
                        SimdFloat *  triplets1,
                        SimdFloat *  triplets2)
{
    for (std::size_t i = 0; i < scalar.simdInternal_.size(); i++)
    {
        triplets0->simdInternal_[i] = scalar.simdInternal_[i / 3];
        triplets1->simdInternal_[i] = scalar.simdInternal_[(i + scalar.simdInternal_.size()) / 3];
        triplets2->simdInternal_[i] = scalar.simdInternal_[(i + 2 * scalar.simdInternal_.size()) / 3];
    }
}

/*! \brief Load 4 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH offsets
 *         specified by a SIMD integer, transpose into 4 SIMD float variables.
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of SIMD variables
 *                    (i.e., 4 for this routine) the input data is packed without
 *                    padding in memory. See the SIMD parameters for exactly
 *                    what memory positions are loaded.
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
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *       base,
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

    for (std::size_t i = 0; i < v0->simdInternal_.size(); i++)
    {
        v0->simdInternal_[i] = base[align * offset.simdInternal_[i]];
        v1->simdInternal_[i] = base[align * offset.simdInternal_[i] + 1];
        v2->simdInternal_[i] = base[align * offset.simdInternal_[i] + 2];
        v3->simdInternal_[i] = base[align * offset.simdInternal_[i] + 3];
    }
}


/*! \brief Load 2 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH offsets
 *         (unaligned) specified by SIMD integer, transpose into 2 SIMD floats.
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of SIMD variables
 *                    (i.e., 2 for this routine) the input data is packed without
 *                    padding in memory. See the SIMD parameters for exactly
 *                    what memory positions are loaded.
 * \param      base   Pointer to the start of the memory.
 * \param      offset SIMD integer type with offsets to the start of each triplet.
 * \param[out] v0     First component, base[align*offset[i]] for each i.
 * \param[out] v1     Second component, base[align*offset[i] + 1] for each i.
 *
 * Since some SIMD architectures cannot handle any unaligned loads, this routine
 * is only available if GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE is 1.
 *
 * \note You should NOT scale offsets before calling this routine; it is
 *       done internally by using the alignment template parameter instead.
 * \note This is a special routine primarily intended for loading Gromacs
 *       table data as efficiently as possible - this is the reason for using
 *       a SIMD offset index, since the result of the  real-to-integer conversion
 *       is present in a SIMD register just before calling this routine.
 */
template <int align>
static inline void gmx_simdcall
gatherLoadUBySimdIntTranspose(const float *       base,
                              SimdFInt32          offset,
                              SimdFloat *         v0,
                              SimdFloat *         v1)
{
    for (std::size_t i = 0; i < v0->simdInternal_.size(); i++)
    {
        v0->simdInternal_[i] = base[align * offset.simdInternal_[i]];
        v1->simdInternal_[i] = base[align * offset.simdInternal_[i] + 1];
    }
}

/*! \brief Load 2 consecutive floats from each of GMX_SIMD_FLOAT_WIDTH offsets
 *         specified by a SIMD integer, transpose into 2 SIMD float variables.
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 *                    When this is identical to the number of SIMD variables
 *                    (i.e., 2 for this routine) the input data is packed without
 *                    padding in memory. See the SIMD parameters for exactly
 *                    what memory positions are loaded.
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
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *       base,
                             SimdFInt32          offset,
                             SimdFloat *         v0,
                             SimdFloat *         v1)
{
    // Base pointer must be aligned to the smaller of 2 elements and float SIMD width
    assert(std::size_t(base) % (std::min(GMX_SIMD_FLOAT_WIDTH, 2)*sizeof(float)) == 0);
    // align parameter must also be a multiple of the above alignment requirement
    assert(align % std::min(GMX_SIMD_FLOAT_WIDTH, 2) == 0);

    for (std::size_t i = 0; i < v0->simdInternal_.size(); i++)
    {
        v0->simdInternal_[i] = base[align * offset.simdInternal_[i]];
        v1->simdInternal_[i] = base[align * offset.simdInternal_[i] + 1];
    }
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
static inline float gmx_simdcall
reduceIncr4ReturnSum(float *           m,
                     SimdFloat         v0,
                     SimdFloat         v1,
                     SimdFloat         v2,
                     SimdFloat         v3)
{
    float sum[4]; // Note that the 4 here corresponds to the 4 m-elements, not any SIMD width

    // Make sure the memory pointer is aligned to the smaller of 4 elements and float SIMD width
    assert(std::size_t(m) % (std::min(GMX_SIMD_FLOAT_WIDTH, 4)*sizeof(float)) == 0);

    sum[0] = reduce(v0);
    sum[1] = reduce(v1);
    sum[2] = reduce(v2);
    sum[3] = reduce(v3);

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
 * To make this more generic, when \ref GMX_SIMD_HAVE_HSIMD_UTIL_REAL is 1,
 * the SIMD implementation provides seven special routines that:
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
 * Available if \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is 1.
 */
static inline SimdFloat gmx_simdcall
loadDualHsimd(const float *  m0,
              const float *  m1)
{
    SimdFloat        a;

    // Make sure the memory pointers are aligned to half float SIMD width
    assert(std::size_t(m0) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);
    assert(std::size_t(m1) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);

    for (std::size_t i = 0; i < a.simdInternal_.size()/2; i++)
    {
        a.simdInternal_[i]                            = m0[i];
        a.simdInternal_[a.simdInternal_.size()/2 + i] = m1[i];
    }
    return a;
}

/*! \brief Load half-SIMD-width float data, spread to both halves.
 *
 * \param m Pointer to memory aligned to half SIMD width.
 *
 * \return SIMD variable with both halves loaded from m..
 *
 * Available if \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is 1.
 */
static inline SimdFloat gmx_simdcall
loadDuplicateHsimd(const float *  m)
{
    SimdFloat        a;

    // Make sure the memory pointer is aligned
    assert(std::size_t(m) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);

    for (std::size_t i = 0; i < a.simdInternal_.size()/2; i++)
    {
        a.simdInternal_[i]                            = m[i];
        a.simdInternal_[a.simdInternal_.size()/2 + i] = a.simdInternal_[i];
    }
    return a;
}

/*! \brief Load two floats, spread 1st in low half, 2nd in high half.
 *
 * \param m Pointer to two adjacent float values.
 *
 * \return SIMD variable where all elements in the low half have been set
 *         to m[0], and all elements in high half to m[1].
 *
 * \note This routine always loads two values and sets the halves separately.
 *       If you want to set all elements to the same value, simply use
 *       the standard (non-half-SIMD) operations.
 *
 * Available if \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is 1.
 */
static inline SimdFloat gmx_simdcall
loadU1DualHsimd(const float *  m)
{
    SimdFloat        a;

    for (std::size_t i = 0; i < a.simdInternal_.size()/2; i++)
    {
        a.simdInternal_[i]                            = m[0];
        a.simdInternal_[a.simdInternal_.size()/2 + i] = m[1];
    }
    return a;
}


/*! \brief Store low & high parts of SIMD float to different locations.
 *
 * \param m0 Pointer to memory aligned to half SIMD width.
 * \param m1 Pointer to memory aligned to half SIMD width.
 * \param a  SIMD variable. Low half should be stored to m0, high to m1.
 *
 * Available if \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is 1.
 */
static inline void gmx_simdcall
storeDualHsimd(float *           m0,
               float *           m1,
               SimdFloat         a)
{
    // Make sure the memory pointers are aligned to half float SIMD width
    assert(std::size_t(m0) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);
    assert(std::size_t(m1) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);

    for (std::size_t i = 0; i < a.simdInternal_.size()/2; i++)
    {
        m0[i] = a.simdInternal_[i];
        m1[i] = a.simdInternal_[a.simdInternal_.size()/2 + i];
    }
}

/*! \brief Add each half of SIMD variable to separate memory adresses
 *
 * \param m0 Pointer to memory aligned to half SIMD width.
 * \param m1 Pointer to memory aligned to half SIMD width.
 * \param a  SIMD variable. Lower half will be added to m0, upper half to m1.
 *
 * The memory must be aligned to half SIMD width.
 *
 * \note The updated m0 value is written before m1 is read from memory, so
 *       the result will be correct even if the memory regions overlap.
 *
 * Available if \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is 1.
 */
static inline void gmx_simdcall
incrDualHsimd(float *           m0,
              float *           m1,
              SimdFloat         a)
{
    // Make sure the memory pointer is aligned to half float SIMD width
    assert(std::size_t(m0) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);
    assert(std::size_t(m1) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);

    for (std::size_t i = 0; i < a.simdInternal_.size()/2; i++)
    {
        m0[i] += a.simdInternal_[i];
    }
    for (std::size_t i = 0; i < a.simdInternal_.size()/2; i++)
    {
        m1[i] += a.simdInternal_[a.simdInternal_.size()/2 + i];
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
 * Available if \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is 1.
 */
static inline void gmx_simdcall
decrHsimd(float *           m,
          SimdFloat         a)
{
    // Make sure the memory pointer is aligned to half float SIMD width
    assert(std::size_t(m) % (GMX_SIMD_FLOAT_WIDTH/2*sizeof(float)) == 0);

    for (std::size_t i = 0; i < a.simdInternal_.size()/2; i++)
    {
        m[i] -= a.simdInternal_[i] + a.simdInternal_[a.simdInternal_.size()/2 + i];
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
 * gatherLoadTranspose(), but just
 * as the other hsimd routines it will pick half-SIMD-width data from base0
 * and put in the lower half, while the upper half comes from base1.
 *
 * For an example, assume the SIMD width is 8, align is 2, that
 * base0 is [A0 A1 B0 B1 C0 C1 D0 D1 ...], and base1 [E0 E1 F0 F1 G0 G1 H0 H1...].
 *
 * Then we will get v0 as [A0 B0 C0 D0 E0 F0 G0 H0] and v1 as [A1 B1 C1 D1 E1 F1 G1 H1].
 *
 * Available if \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is 1.
 */
template <int align>
static inline void gmx_simdcall
gatherLoadTransposeHsimd(const float  *       base0,
                         const float  *       base1,
                         const std::int32_t   offset[],
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

    for (std::size_t i = 0; i < v0->simdInternal_.size()/2; i++)
    {
        v0->simdInternal_[i] = base0[align * offset[i]];
        v1->simdInternal_[i] = base0[align * offset[i] + 1];
        v0->simdInternal_[v0->simdInternal_.size()/2 + i] = base1[align * offset[i]];
        v1->simdInternal_[v1->simdInternal_.size()/2 + i] = base1[align * offset[i] + 1];
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
 *      reduceIncr4ReturnSum(). The only difference is that the
 *      four half-SIMD inputs needed are present in the low/high halves of the
 *      two SIMD arguments.
 *
 * Available if \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT is 1.
 */
static inline float gmx_simdcall
reduceIncr4ReturnSumHsimd(float *            m,
                          SimdFloat          v0,
                          SimdFloat          v1)
{
    // The 4 here corresponds to the 4 elements in memory, not any SIMD width
    float sum[4] = { 0.0f, 0.0f, 0.0f, 0.0f };

    for (std::size_t i = 0; i < v0.simdInternal_.size()/2; i++)
    {
        sum[0] += v0.simdInternal_[i];
        sum[1] += v0.simdInternal_[v0.simdInternal_.size()/2 + i];
        sum[2] += v1.simdInternal_[i];
        sum[3] += v1.simdInternal_[v1.simdInternal_.size()/2 + i];
    }

    // Make sure the memory pointer is aligned to the smaller of 4 elements and float SIMD width
    assert(std::size_t(m) % (std::min(GMX_SIMD_FLOAT_WIDTH, 4)*sizeof(float)) == 0);

    m[0] += sum[0];
    m[1] += sum[1];
    m[2] += sum[2];
    m[3] += sum[3];

    return sum[0] + sum[1] + sum[2] + sum[3];
}

#if GMX_SIMD_FLOAT_WIDTH > 8 || defined DOXYGEN
/*! \brief Load N floats and duplicate them 4 times each.
 *
 * \param m Pointer to unaligned memory
 *
 * \return SIMD variable with N floats from m duplicated 4x.
 *
 * Available if \ref GMX_SIMD_HAVE_4NSIMD_UTIL_FLOAT is 1.
 * N is GMX_SIMD_FLOAT_WIDTH/4. Duplicated values are
 * contigous and different values are 4 positions in SIMD
 * apart.
 */
static inline SimdFloat gmx_simdcall
loadUNDuplicate4(const float* m)
{
    SimdFloat        a;
    for (std::size_t i = 0; i < a.simdInternal_.size()/4; i++)
    {
        a.simdInternal_[i*4]   = m[i];
        a.simdInternal_[i*4+1] = m[i];
        a.simdInternal_[i*4+2] = m[i];
        a.simdInternal_[i*4+3] = m[i];
    }
    return a;
}

/*! \brief Load 4 floats and duplicate them N times each.
 *
 * \param m Pointer to memory aligned to 4 floats
 *
 * \return SIMD variable with 4 floats from m duplicated Nx.
 *
 * Available if \ref GMX_SIMD_HAVE_4NSIMD_UTIL_FLOAT is 1.
 * N is GMX_SIMD_FLOAT_WIDTH/4. Different values are
 * contigous and same values are 4 positions in SIMD
 * apart.
 */
static inline SimdFloat gmx_simdcall
load4DuplicateN(const float* m)
{
    SimdFloat        a;
    for (std::size_t i = 0; i < a.simdInternal_.size()/4; i++)
    {
        a.simdInternal_[i*4]   = m[0];
        a.simdInternal_[i*4+1] = m[1];
        a.simdInternal_[i*4+2] = m[2];
        a.simdInternal_[i*4+3] = m[3];
    }
    return a;
}
#endif

#if GMX_SIMD_FLOAT_WIDTH >= 8 || defined DOXYGEN
/*! \brief Load floats in blocks of 4 at fixed offsets
 *
 * \param m Pointer to unaligned memory
 * \param offset Offset in memory between input blocks of 4
 *
 * \return SIMD variable with floats from m.
 *
 * Available if \ref GMX_SIMD_HAVE_4NSIMD_UTIL_FLOAT is 1.
 * Blocks of 4 floats are loaded from m+n*offset where n
 * is the n-th block of 4 floats.
 */
static inline SimdFloat gmx_simdcall
loadU4NOffset(const float* m, int offset)
{
    SimdFloat        a;
    for (std::size_t i = 0; i < a.simdInternal_.size()/4; i++)
    {
        a.simdInternal_[i*4]   = m[offset*i + 0];
        a.simdInternal_[i*4+1] = m[offset*i + 1];
        a.simdInternal_[i*4+2] = m[offset*i + 2];
        a.simdInternal_[i*4+3] = m[offset*i + 3];
    }
    return a;
}
#endif

/*! \} */

/*! \} */
/*! \endcond */

}      // namespace gmx

#endif // GMX_SIMD_IMPL_REFERENCE_UTIL_FLOAT_H
