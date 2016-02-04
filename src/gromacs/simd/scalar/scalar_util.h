/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
#ifndef GMX_SIMD_SCALAR_UTIL_H
#define GMX_SIMD_SCALAR_UTIL_H

#include <cmath>

/*! \libinternal \file
 *
 * \brief Scalar utility functions mimicking GROMACS SIMD utility functions
 *
 * These versions make it possible to write functions that are templated with
 * either a SIMD or scalar type. While some of these functions might not appear
 * SIMD-specific, we have placed them here because the only reason to use these
 * instead of generic function is in templated combined SIMD/non-SIMD code.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_simd
 */

namespace gmx
{

/*****************************************************************************
 *   Single-precision utility load/store functions mimicking SIMD versions   *
 *****************************************************************************/

/*! \brief Load 4 consecutive floats from base/offset into four variables
 *
 * \tparam     align  Alignment of the memory from which we read.
 * \param      base   Pointer to the start of the memory area
 * \param      offset Index to data.
 * \param[out] v0     1st float, base[align*offset[0]].
 * \param[out] v1     2nd float, base[align*offset[0] + 1].
 * \param[out] v2     3rd float, base[align*offset[0] + 2].
 * \param[out] v3     4th float, base[align*offset[0] + 3].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
gatherLoadTranspose(const float  *        base,
                    const std::int32_t    offset[],
                    float *               v0,
                    float *               v1,
                    float *               v2,
                    float *               v3)
{
    *v0 = base[align*offset[0]];
    *v1 = base[align*offset[0]+1];
    *v2 = base[align*offset[0]+2];
    *v3 = base[align*offset[0]+3];
}

/*! \brief Load 2 consecutive floats from base/offset into four variables
 *
 * \tparam     align  Alignment of the memory from which we read.
 * \param      base   Pointer to the start of the memory area
 * \param      offset Index to data.
 * \param[out] v0     1st float, base[align*offset[0]].
 * \param[out] v1     2nd float, base[align*offset[0] + 1].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
gatherLoadTranspose(const float  *        base,
                    const std::int32_t    offset[],
                    float *               v0,
                    float *               v1)
{
    *v0 = base[align*offset[0]];
    *v1 = base[align*offset[0]+1];
}


/*! \brief Load 3 consecutive floats from base/offsets, store into three vars.
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param      base   Pointer to the start of the memory area
 * \param      offset Offset to the start of data.
 * \param[out] v0     1st value, base[align*offset[0]].
 * \param[out] v1     2nd value, base[align*offset[0] + 1].
 * \param[out] v2     3rd value, base[align*offset[0] + 2].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
gatherLoadUTranspose(const float  *        base,
                     const std::int32_t    offset[],
                     float *               v0,
                     float *               v1,
                     float *               v2)
{
    *v0 = base[align*offset[0]];
    *v1 = base[align*offset[0]+1];
    *v2 = base[align*offset[0]+2];
}

/*! \brief Store 3 floats to 3 to base/offset.
 *
 * \tparam     align  Alignment of the memory to which we write, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Offset to the start of triplet.
 * \param      v0     1st value, written to base[align*offset[0]].
 * \param      v1     2nd value, written to base[align*offset[0] + 1].
 * \param      v2     3rd value, written to base[align*offset[0] + 2].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
transposeScatterStoreU(float  *              base,
                       const std::int32_t    offset[],
                       float                 v0,
                       float                 v1,
                       float                 v2)
{
    base[align*offset[0]]   = v0;
    base[align*offset[0]+1] = v1;
    base[align*offset[0]+2] = v2;
}

/*! \brief Add 3 floats to base/offset.
 *
 * \tparam     align  Alignment of the memory to which we write, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Offset to the start of triplet.
 * \param      v0     1st value, added to base[align*offset[0]].
 * \param      v1     2nd value, added to base[align*offset[0] + 1].
 * \param      v2     3rd value, added to base[align*offset[0] + 2].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
transposeScatterIncrU(float  *              base,
                      const std::int32_t    offset[],
                      float                 v0,
                      float                 v1,
                      float                 v2)
{
    base[align*offset[0]]   += v0;
    base[align*offset[0]+1] += v1;
    base[align*offset[0]+2] += v2;
}

/*! \brief Subtract 3 floats from base/offset.
 *
 * \tparam     align  Alignment of the memory to which we write, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Offset to the start of triplet.
 * \param      v0     1st value, subtracted from base[align*offset[0]].
 * \param      v1     2nd value, subtracted from base[align*offset[0] + 1].
 * \param      v2     3rd value, subtracted from base[align*offset[0] + 2].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
transposeScatterDecrU(float  *              base,
                      const std::int32_t    offset[],
                      float                 v0,
                      float                 v1,
                      float                 v2)
{
    base[align*offset[0]]   -= v0;
    base[align*offset[0]+1] -= v1;
    base[align*offset[0]+2] -= v2;
}

/*! \brief Copy single float to three variables.
 *
 * \param      scalar    Floating-point input.
 * \param[out] triplets0 Copy 1.
 * \param[out] triplets1 Copy 2.
 * \param[out] triplets2 Copy 3.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
expandScalarsToTriplets(float    scalar,
                        float *  triplets0,
                        float *  triplets1,
                        float *  triplets2)
{
    *triplets0 = scalar;
    *triplets1 = scalar;
    *triplets2 = scalar;
}

/*! \brief Load 4 floats from base/offsets and store into variables.
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param      base   Aligned pointer to the start of the memory.
 * \param      offset Integer type with offset to the start of each triplet.
 * \param[out] v0     First float, base[align*offset[0]].
 * \param[out] v1     Second float, base[align*offset[0] + 1].
 * \param[out] v2     Third float, base[align*offset[0] + 2].
 * \param[out] v3     Fourth float, base[align*offset[0] + 3].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
gatherLoadBySimdIntTranspose(const float *       base,
                             std::int32_t        offset,
                             float *             v0,
                             float *             v1,
                             float *             v2,
                             float *             v3)
{
    *v0 = base[align*offset];
    *v1 = base[align*offset+1];
    *v2 = base[align*offset+2];
    *v3 = base[align*offset+3];
}

/*! \brief Load 2 floats from base/offsets and store into variables (unaligned).
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param      base   Aligned pointer to the start of the memory.
 * \param      offset Integer type with offset to the start of each triplet.
 * \param[out] v0     First float, base[align*offset[0]].
 * \param[out] v1     Second float, base[align*offset[0] + 1].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
gatherLoadUBySimdIntTranspose(const float *       base,
                              std::int32_t        offset,
                              float *             v0,
                              float *             v1)
{
    *v0 = base[align*offset];
    *v1 = base[align*offset+1];
}

/*! \brief Load 2 floats from base/offsets and store into variables (aligned).
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param      base   Aligned pointer to the start of the memory.
 * \param      offset Integer type with offset to the start of each triplet.
 * \param[out] v0     First float, base[align*offset[0]].
 * \param[out] v1     Second float, base[align*offset[0] + 1].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
gatherLoadBySimdIntTranspose(const float *       base,
                             std::int32_t        offset,
                             float *             v0,
                             float *             v1)
{
    *v0 = base[align*offset];
    *v1 = base[align*offset+1];
}

/*! \brief Add each float to four consecutive memory locations, return sum.
 *
 * \param m   Pointer to memory where four floats should be incremented
 * \param v0  float to be added to m[0]
 * \param v1  float to be added to m[1]
 * \param v2  float to be added to m[2]
 * \param v3  float to be added to m[3]
 *
 * \return v0+v1+v2+v3.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
reduceIncr4ReturnSum(float *           m,
                     float             v0,
                     float             v1,
                     float             v2,
                     float             v3)
{
    m[0] += v0;
    m[1] += v1;
    m[2] += v2;
    m[3] += v3;

    return v0 + v1 + v2 + v3;
}


/*****************************************************************************
 *   Double-precision utility load/store functions mimicking SIMD versions   *
 *****************************************************************************/

/*! \brief Load 4 consecutive doubles from base/offset into four variables
 *
 * \tparam     align  Alignment of the memory from which we read.
 * \param      base   Pointer to the start of the memory area
 * \param      offset Index to data.
 * \param[out] v0     1st double, base[align*offset[0]].
 * \param[out] v1     2nd double, base[align*offset[0] + 1].
 * \param[out] v2     3rd double, base[align*offset[0] + 2].
 * \param[out] v3     4th double, base[align*offset[0] + 3].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
gatherLoadTranspose(const double  *        base,
                    const std::int32_t     offset[],
                    double *               v0,
                    double *               v1,
                    double *               v2,
                    double *               v3)
{
    *v0 = base[align*offset[0]];
    *v1 = base[align*offset[0]+1];
    *v2 = base[align*offset[0]+2];
    *v3 = base[align*offset[0]+3];
}

/*! \brief Load 2 consecutive doubles from base/offset into four variables
 *
 * \tparam     align  Alignment of the memory from which we read.
 * \param      base   Pointer to the start of the memory area
 * \param      offset Index to data.
 * \param[out] v0     1st double, base[align*offset[0]].
 * \param[out] v1     2nd double, base[align*offset[0] + 1].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
gatherLoadTranspose(const double  *        base,
                    const std::int32_t     offset[],
                    double *               v0,
                    double *               v1)
{
    *v0 = base[align*offset[0]];
    *v1 = base[align*offset[0]+1];
}


/*! \brief Load 3 consecutive doubles from base/offsets, store into three vars.
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param      base   Pointer to the start of the memory area
 * \param      offset Offset to the start of data.
 * \param[out] v0     1st double, base[align*offset[0]].
 * \param[out] v1     2nd double, base[align*offset[0] + 1].
 * \param[out] v2     3rd double, base[align*offset[0] + 2].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
gatherLoadUTranspose(const double  *        base,
                     const std::int32_t     offset[],
                     double *               v0,
                     double *               v1,
                     double *               v2)
{
    *v0 = base[align*offset[0]];
    *v1 = base[align*offset[0]+1];
    *v2 = base[align*offset[0]+2];
}

/*! \brief Store 3 doubles to 3 to base/offset.
 *
 * \tparam     align  Alignment of the memory to which we write, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Offset to the start of triplet.
 * \param      v0     1st value, written to base[align*offset[0]].
 * \param      v1     2nd value, written to base[align*offset[0] + 1].
 * \param      v2     3rd value, written to base[align*offset[0] + 2].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
transposeScatterStoreU(double  *             base,
                       const std::int32_t    offset[],
                       double                v0,
                       double                v1,
                       double                v2)
{
    base[align*offset[0]]   = v0;
    base[align*offset[0]+1] = v1;
    base[align*offset[0]+2] = v2;
}

/*! \brief Add 3 doubles to base/offset.
 *
 * \tparam     align  Alignment of the memory to which we write, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Offset to the start of triplet.
 * \param      v0     1st value, added to base[align*offset[0]].
 * \param      v1     2nd value, added to base[align*offset[0] + 1].
 * \param      v2     3rd value, added to base[align*offset[0] + 2].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
transposeScatterIncrU(double  *              base,
                      const std::int32_t     offset[],
                      double                 v0,
                      double                 v1,
                      double                 v2)
{
    base[align*offset[0]]   += v0;
    base[align*offset[0]+1] += v1;
    base[align*offset[0]+2] += v2;
}

/*! \brief Subtract 3 doubles from base/offset.
 *
 * \tparam     align  Alignment of the memory to which we write, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Offset to the start of triplet.
 * \param      v0     1st value, subtracted from base[align*offset[0]].
 * \param      v1     2nd value, subtracted from base[align*offset[0] + 1].
 * \param      v2     3rd value, subtracted from base[align*offset[0] + 2].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
transposeScatterDecrU(double  *              base,
                      const std::int32_t     offset[],
                      double                 v0,
                      double                 v1,
                      double                 v2)
{
    base[align*offset[0]]   -= v0;
    base[align*offset[0]+1] -= v1;
    base[align*offset[0]+2] -= v2;
}

/*! \brief Copy single double to three variables.
 *
 * \param      scalar    Floating-point input.
 * \param[out] triplets0 Copy 1.
 * \param[out] triplets1 Copy 2.
 * \param[out] triplets2 Copy 3.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
expandScalarsToTriplets(double    scalar,
                        double *  triplets0,
                        double *  triplets1,
                        double *  triplets2)
{
    *triplets0 = scalar;
    *triplets1 = scalar;
    *triplets2 = scalar;
}

/*! \brief Load 4 doubles from base/offsets and store into variables.
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param      base   Aligned pointer to the start of the memory.
 * \param      offset Integer type with offset to the start of each triplet.
 * \param[out] v0     First double, base[align*offset[0]].
 * \param[out] v1     Second double, base[align*offset[0] + 1].
 * \param[out] v2     Third double, base[align*offset[0] + 2].
 * \param[out] v3     Fourth double, base[align*offset[0] + 3].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
gatherLoadBySimdIntTranspose(const double *       base,
                             std::int32_t         offset,
                             double *             v0,
                             double *             v1,
                             double *             v2,
                             double *             v3)
{
    *v0 = base[align*offset];
    *v1 = base[align*offset+1];
    *v2 = base[align*offset+2];
    *v3 = base[align*offset+3];
}

/*! \brief Load 2 doubles from base/offsets and store into variables (unaligned).
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param      base   Aligned pointer to the start of the memory.
 * \param      offset Integer type with offset to the start of each triplet.
 * \param[out] v0     First double, base[align*offset[0]].
 * \param[out] v1     Second double, base[align*offset[0] + 1].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
gatherLoadUBySimdIntTranspose(const double *       base,
                              std::int32_t         offset,
                              double *             v0,
                              double *             v1)
{
    *v0 = base[align*offset];
    *v1 = base[align*offset+1];
}

/*! \brief Load 2 doubles from base/offsets and store into variables (aligned).
 *
 * \tparam     align  Alignment of the memory from which we read, i.e. distance
 *                    (measured in elements, not bytes) between index points.
 * \param      base   Aligned pointer to the start of the memory.
 * \param      offset Integer type with offset to the start of each triplet.
 * \param[out] v0     First double, base[align*offset[0]].
 * \param[out] v1     Second double, base[align*offset[0] + 1].
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <int align>
static inline void
gatherLoadBySimdIntTranspose(const double *      base,
                             std::int32_t        offset,
                             double *            v0,
                             double *            v1)
{
    *v0 = base[align*offset];
    *v1 = base[align*offset+1];
}

/*! \brief Add each double to four consecutive memory locations, return sum.
 *
 * \param m   Pointer to memory where four floats should be incremented
 * \param v0  double to be added to m[0]
 * \param v1  double to be added to m[1]
 * \param v2  double to be added to m[2]
 * \param v3  double to be added to m[3]
 *
 * \return v0+v1+v2+v3.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
reduceIncr4ReturnSum(double *           m,
                     double             v0,
                     double             v1,
                     double             v2,
                     double             v3)
{
    m[0] += v0;
    m[1] += v1;
    m[2] += v2;
    m[3] += v3;

    return v0 + v1 + v2 + v3;
}

} // namespace gmx


#endif // GMX_SIMD_SCALAR_UTIL_H
