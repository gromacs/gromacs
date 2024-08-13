/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 * Declares the SIMD cluster pair distance kernel and helpers
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_PAIRLIST_SIMD_KERNEL_H
#define GMX_NBNXM_PAIRLIST_SIMD_KERNEL_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "pairlistwork.h"

namespace gmx
{
struct NbnxmPairlistCpu;
struct NbnxmPairlistCpuWork;
class Grid;

//! Copies PBC shifted i-cell packed atom coordinates to working array for the 4xM layout
void setICellCoordinatesSimd4xM(int                   ci,
                                const RVec&           shift,
                                int gmx_unused        stride,
                                const real*           x,
                                NbnxmPairlistCpuWork* work);

//! Copies PBC shifted i-cell packed atom coordinates to working array for the 2xMM layout
void setICellCoordinatesSimd2xMM(int                   ci,
                                 const RVec&           shift,
                                 int gmx_unused        stride,
                                 const real*           x,
                                 NbnxmPairlistCpuWork* work);

/*! \brief SIMD code for checking and adding cluster-pairs to the list using the 4xM layout.
 *
 * Checks bounding box distances and possibly atom pair distances.
 * This is an accelerated version of make_cluster_list_simple.
 *
 * \param[in]     jGrid               The j-grid
 * \param[in,out] nbl                 The pair-list to store the cluster pairs in
 * \param[in]     icluster            The index of the i-cluster
 * \param[in]     firstCell           The first cluster in the j-range, using i-cluster size indexing
 * \param[in]     lastCell            The last cluster in the j-range, using i-cluster size indexing
 * \param[in]     excludeSubDiagonal  Exclude atom pairs with i-index > j-index
 * \param[in]     x_j                 Coordinates for the j-atom, in SIMD packed format
 * \param[in]     rlist2              The squared list cut-off
 * \param[in]     rbb2                The squared cut-off for putting cluster-pairs in the list based on bounding box distance only
 * \param[in,out] numDistanceChecks   The number of distance checks performed
 */
void makeClusterListSimd4xM(const Grid&              jGrid,
                            NbnxnPairlistCpu*        nbl,
                            int                      icluster,
                            int                      firstCell,
                            int                      lastCell,
                            bool                     excludeSubDiagonal,
                            const real* gmx_restrict x_j,
                            real                     rlist2,
                            float                    rbb2,
                            int* gmx_restrict        numDistanceChecks);

/*! \brief SIMD code for checking and adding cluster-pairs to the list using the 2xMM layout.
 *
 * Checks bounding box distances and possibly atom pair distances.
 * This is an accelerated version of make_cluster_list_simple.
 *
 * \param[in]     jGrid               The j-grid
 * \param[in,out] nbl                 The pair-list to store the cluster pairs in
 * \param[in]     icluster            The index of the i-cluster
 * \param[in]     firstCell           The first cluster in the j-range, using i-cluster size indexing
 * \param[in]     lastCell            The last cluster in the j-range, using i-cluster size indexing
 * \param[in]     excludeSubDiagonal  Exclude atom pairs with i-index > j-index
 * \param[in]     x_j                 Coordinates for the j-atom, in SIMD packed format
 * \param[in]     rlist2              The squared list cut-off
 * \param[in]     rbb2                The squared cut-off for putting cluster-pairs in the list based on bounding box distance only
 * \param[in,out] numDistanceChecks   The number of distance checks performed
 */
void makeClusterListSimd2xMM(const Grid&              jGrid,
                             NbnxnPairlistCpu*        nbl,
                             int                      icluster,
                             int                      firstCell,
                             int                      lastCell,
                             bool                     excludeSubDiagonal,
                             const real* gmx_restrict x_j,
                             real                     rlist2,
                             float                    rbb2,
                             int* gmx_restrict        numDistanceChecks);

} // namespace gmx

#endif // GMX_NBNXM_PAIRLIST_SIMD_KERNEL_H
