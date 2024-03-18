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

#include "boundingboxes.h"
#include "grid.h"

namespace Nbnxm
{

#if !NBNXN_SEARCH_BB_SIMD4

/*! \brief Plain C code calculating the distance^2 between two bounding boxes in xyz0 format
 *
 * \param[in] bb_i  First bounding box
 * \param[in] bb_j  Second bounding box
 */
gmx_unused static float clusterBoundingBoxDistance2(const BoundingBox& bb_i, const BoundingBox& bb_j)
{
    float dl  = bb_i.lower.x - bb_j.upper.x;
    float dh  = bb_j.lower.x - bb_i.upper.x;
    float dm  = std::max(dl, dh);
    float dm0 = std::max(dm, 0.0F);
    float d2  = dm0 * dm0;

    dl  = bb_i.lower.y - bb_j.upper.y;
    dh  = bb_j.lower.y - bb_i.upper.y;
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0F);
    d2 += dm0 * dm0;

    dl  = bb_i.lower.z - bb_j.upper.z;
    dh  = bb_j.lower.z - bb_i.upper.z;
    dm  = std::max(dl, dh);
    dm0 = std::max(dm, 0.0F);
    d2 += dm0 * dm0;

    return d2;
}

#else /* NBNXN_SEARCH_BB_SIMD4 */

/*! \brief 4-wide SIMD code calculating the distance^2 between two bounding boxes in xyz0 format
 *
 * \param[in] bb_i  First bounding box, should be aligned for 4-wide SIMD
 * \param[in] bb_j  Second bounding box, should be aligned for 4-wide SIMD
 */
static float clusterBoundingBoxDistance2(const BoundingBox& bb_i, const BoundingBox& bb_j)
{
    using namespace gmx;

    const Simd4Float bb_i_S0 = load4(bb_i.lower.ptr());
    const Simd4Float bb_i_S1 = load4(bb_i.upper.ptr());
    const Simd4Float bb_j_S0 = load4(bb_j.lower.ptr());
    const Simd4Float bb_j_S1 = load4(bb_j.upper.ptr());

    const Simd4Float dl_S = bb_i_S0 - bb_j_S1;
    const Simd4Float dh_S = bb_j_S0 - bb_i_S1;

    const Simd4Float dm_S  = max(dl_S, dh_S);
    const Simd4Float dm0_S = max(dm_S, simd4SetZeroF());

    return dotProduct(dm0_S, dm0_S);
}

/* Calculate bb bounding distances of bb_i[si,...,si+3] and store them in d2 */
template<int boundingBoxStart>
static inline void gmx_simdcall clusterBoundingBoxDistance2_xxxx_simd4_inner(const float* bb_i,
                                                                             float*       d2,
                                                                             const gmx::Simd4Float xj_l,
                                                                             const gmx::Simd4Float yj_l,
                                                                             const gmx::Simd4Float zj_l,
                                                                             const gmx::Simd4Float xj_h,
                                                                             const gmx::Simd4Float yj_h,
                                                                             const gmx::Simd4Float zj_h)
{
    using namespace gmx;

    constexpr int stride = c_packedBoundingBoxesDimSize;

    const int shi = boundingBoxStart * c_numBoundingBoxBounds1D * DIM;

    const Simd4Float zero = setZero();

    const Simd4Float xi_l = load4(bb_i + shi + 0 * stride);
    const Simd4Float yi_l = load4(bb_i + shi + 1 * stride);
    const Simd4Float zi_l = load4(bb_i + shi + 2 * stride);
    const Simd4Float xi_h = load4(bb_i + shi + 3 * stride);
    const Simd4Float yi_h = load4(bb_i + shi + 4 * stride);
    const Simd4Float zi_h = load4(bb_i + shi + 5 * stride);

    const Simd4Float dx_0 = xi_l - xj_h;
    const Simd4Float dy_0 = yi_l - yj_h;
    const Simd4Float dz_0 = zi_l - zj_h;

    const Simd4Float dx_1 = xj_l - xi_h;
    const Simd4Float dy_1 = yj_l - yi_h;
    const Simd4Float dz_1 = zj_l - zi_h;

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

    store4(d2 + boundingBoxStart, d2t);
}

/* 4-wide SIMD code for nsi bb distances for bb format xxxxyyyyzzzz */
gmx_unused static void clusterBoundingBoxDistance2_xxxx_simd4(const float* bb_j,
                                                              const int    nsi,
                                                              const float* bb_i,
                                                              float*       d2)
{
    using namespace gmx;

    constexpr int stride = c_packedBoundingBoxesDimSize;

    const Simd4Float xj_l = Simd4Float(bb_j[0 * stride]);
    const Simd4Float yj_l = Simd4Float(bb_j[1 * stride]);
    const Simd4Float zj_l = Simd4Float(bb_j[2 * stride]);
    const Simd4Float xj_h = Simd4Float(bb_j[3 * stride]);
    const Simd4Float yj_h = Simd4Float(bb_j[4 * stride]);
    const Simd4Float zj_h = Simd4Float(bb_j[5 * stride]);

    /* Here we "loop" over si (0,stride) from 0 to nsi with step stride.
     * But as we know the number of iterations is 1 or 2, we unroll manually.
     */
    clusterBoundingBoxDistance2_xxxx_simd4_inner<0>(bb_i, d2, xj_l, yj_l, zj_l, xj_h, yj_h, zj_h);
    if (stride < nsi)
    {
        clusterBoundingBoxDistance2_xxxx_simd4_inner<stride>(bb_i, d2, xj_l, yj_l, zj_l, xj_h, yj_h, zj_h);
    }
}

#endif /* NBNXN_SEARCH_BB_SIMD4 */

} // namespace Nbnxm

#endif // GMX_NBNXM_BOUNDINGBOXDISTANCE_H
