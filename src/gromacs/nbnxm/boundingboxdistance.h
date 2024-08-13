/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \internal \file
 *
 * \brief
 * Declares and defines functions for computing distances between bounding boxes
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_BOUNDINGBOXDISTANCE_H
#define GMX_NBNXM_BOUNDINGBOXDISTANCE_H

#include "config.h"

#include "gromacs/simd/simd.h"
#include "gromacs/utility/gmxassert.h"

#include "boundingbox.h"
#include "boundingbox_simd.h"

namespace gmx
{

//! Loads a corner of a bounding box into a float vector
template<typename T>
std::enable_if_t<std::is_same_v<T, gmx::BasicVector<float>>, gmx::BasicVector<float>> static inline loadBoundingBoxCorner(
        const BoundingBox::Corner& corner)
{
    return { corner.x, corner.y, corner.z };
}

#if NBNXN_SEARCH_BB_SIMD4

/*! Loads a corner of a bounding box into a 4-wide SIMD register
 *
 * \param[in] corner  A bounding box corner, should be SIMD4 aligned
 */
template<typename T>
std::enable_if_t<std::is_same_v<T, gmx::Simd4Float>, gmx::Simd4Float> static inline gmx_simdcall
loadBoundingBoxCorner(const BoundingBox::Corner& corner)
{
    static_assert(sizeof(BoundingBox::Corner) == 4 * sizeof(float));

    GMX_ASSERT(std::size_t(corner.ptr()) % (4 * sizeof(*corner.ptr())) == 0,
               "Bounding box should be SIMD4 aligned");

    return gmx::load4(corner.ptr());
}

#endif

//! Return the element-wise max of two 3-float vectors, needed to share code with SIMD
static inline gmx::BasicVector<float> max(const gmx::BasicVector<float>& v1,
                                          const gmx::BasicVector<float>& v2)
{
    return gmx::elementWiseMax(v1, v2);
}

//! Return the dot product of two 3-float vectors, needed to share code with SIMD
static inline float dotProduct(const gmx::BasicVector<float>& v1, const gmx::BasicVector<float>& v2)
{
    return gmx::dot(v1, v2);
}

/*! \brief Returns the distance^2 between two bounding boxes
 *
 * Uses 4-wide SIMD operations when available.
 *
 * \param[in] bb_i  First bounding box, has to be aligned for 4-wide SIMD
 * \param[in] bb_j  Second bounding box, has to be aligned for 4-wide SIMD
 */
static inline float clusterBoundingBoxDistance2(const BoundingBox& bb_i, const BoundingBox& bb_j)
{
#if NBNXN_SEARCH_BB_SIMD4
    using T = Simd4Float;

    const T zero = simd4SetZeroF();
#else
    using T = BasicVector<float>;

    const T zero({ 0.0f, 0.0f, 0.0f });
#endif

    const T iLowerCorner = loadBoundingBoxCorner<T>(bb_i.lower);
    const T iUpperCorner = loadBoundingBoxCorner<T>(bb_i.upper);
    const T jLowerCorner = loadBoundingBoxCorner<T>(bb_j.lower);
    const T jUpperCorner = loadBoundingBoxCorner<T>(bb_j.upper);

    const T difference0 = iLowerCorner - jUpperCorner;
    const T difference1 = jLowerCorner - iUpperCorner;

    const T maxDifference     = max(difference0, difference1);
    const T maxDifferenceZero = max(maxDifference, zero);

    return dotProduct(maxDifferenceZero, maxDifferenceZero);
}

#if NBNXN_SEARCH_BB_SIMD4

/*! Calculates bounding box distances^2 of four i-bounding-boxes with one j-bounding-box
 *
 * \tparam     boundingBoxStart  The offset for reading bounding boxes and storing distances
 * \param[in]  iBoundingBoxes    Bounding boxes, entries \p boundingBoxStart to \p boundingBoxStart+4 are used
 * \param[in]  jLowerCorner      Lower X/Y/Z corners of the j-bounding box
 * \param[in]  jUpperCorner      Upper X/Y/Z corners of the j-bounding box
 * \param[out] distancesSquared  The four squared distances are returned in this list at offset \p boundingBoxStart
 */
template<int boundingBoxStart>
static inline void gmx_simdcall
clusterBoundingBoxDistance2_xxxx_simd4_inner(const float*                            iBoundingBoxes,
                                             const std::array<gmx::Simd4Float, DIM>& jLowerCorner,
                                             const std::array<gmx::Simd4Float, DIM>& jUpperCorner,
                                             float* distancesSquared)
{
    constexpr int stride = c_packedBoundingBoxesDimSize;

    const int shi = boundingBoxStart * c_numBoundingBoxBounds1D * DIM;

    const Simd4Float zero = setZero();

    const Simd4Float iLowerCornerX = load4(iBoundingBoxes + shi + 0 * stride);
    const Simd4Float iLowerCornerY = load4(iBoundingBoxes + shi + 1 * stride);
    const Simd4Float iLowerCornerZ = load4(iBoundingBoxes + shi + 2 * stride);
    const Simd4Float iUpperCornerX = load4(iBoundingBoxes + shi + 3 * stride);
    const Simd4Float iUpperCornerY = load4(iBoundingBoxes + shi + 4 * stride);
    const Simd4Float iUpperCornerZ = load4(iBoundingBoxes + shi + 5 * stride);

    const Simd4Float dx_0 = iLowerCornerX - jUpperCorner[XX];
    const Simd4Float dy_0 = iLowerCornerY - jUpperCorner[YY];
    const Simd4Float dz_0 = iLowerCornerZ - jUpperCorner[ZZ];

    const Simd4Float dx_1 = jLowerCorner[XX] - iUpperCornerX;
    const Simd4Float dy_1 = jLowerCorner[YY] - iUpperCornerY;
    const Simd4Float dz_1 = jLowerCorner[ZZ] - iUpperCornerZ;

    const Simd4Float mx = max(dx_0, dx_1);
    const Simd4Float my = max(dy_0, dy_1);
    const Simd4Float mz = max(dz_0, dz_1);

    const Simd4Float m0x = max(mx, zero);
    const Simd4Float m0y = max(my, zero);
    const Simd4Float m0z = max(mz, zero);

    const Simd4Float d2x = m0x * m0x;
    const Simd4Float d2y = m0y * m0y;
    const Simd4Float d2z = m0z * m0z;

    const Simd4Float d2s = d2x + d2y;
    const Simd4Float d2t = d2s + d2z;

    store4(distancesSquared + boundingBoxStart, d2t);
}

//! 4-wide SIMD code for nsi bb distances for bb format xxxxyyyyzzzz
gmx_unused static void clusterBoundingBoxDistance2_xxxx_simd4(const float* bb_j,
                                                              const int    nsi,
                                                              const float* bb_i,
                                                              float*       d2)
{
    constexpr int stride = c_packedBoundingBoxesDimSize;

    const std::array<Simd4Float, DIM> jLowerCorner = { Simd4Float(bb_j[0 * stride]),
                                                       Simd4Float(bb_j[1 * stride]),
                                                       Simd4Float(bb_j[2 * stride]) };
    const std::array<Simd4Float, DIM> jUpperCorner = { Simd4Float(bb_j[3 * stride]),
                                                       Simd4Float(bb_j[4 * stride]),
                                                       Simd4Float(bb_j[5 * stride]) };

    /* Here we "loop" over si (0,stride) from 0 to nsi with step stride.
     * But as we know the number of iterations is 1 or 2, we unroll manually.
     */
    clusterBoundingBoxDistance2_xxxx_simd4_inner<0>(bb_i, jLowerCorner, jUpperCorner, d2);
    if (stride < nsi)
    {
        clusterBoundingBoxDistance2_xxxx_simd4_inner<stride>(bb_i, jLowerCorner, jUpperCorner, d2);
    }
}

#endif /* NBNXN_SEARCH_BB_SIMD4 */

} // namespace gmx

#endif // GMX_NBNXM_BOUNDINGBOXDISTANCE_H
