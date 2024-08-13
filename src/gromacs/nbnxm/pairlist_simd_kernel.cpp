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
 * Defines the SIMD cluster pair distance kernel and helpers
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#include "pairlist_simd_kernel.h"

#include "config.h"

#include <algorithm>
#include <array>
#include <memory>
#include <type_traits>

#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/nbnxm/boundingbox.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

#include "boundingboxdistance.h"
#include "clusterdistancekerneltype.h"
#include "grid.h"
#include "nbnxm_simd.h"
#include "pairlist_imask.h"
#include "pairlistwork.h"

namespace gmx
{

template<ClusterDistanceKernelType kernelType>
gmx_unused static constexpr int sc_iClusterSizeSimd()
{
    static_assert(kernelType == ClusterDistanceKernelType::CpuSimd_4xM
                  || kernelType == ClusterDistanceKernelType::CpuSimd_2xMM);

    return 4;
}

template<ClusterDistanceKernelType kernelType>
gmx_unused static constexpr int sc_jClusterSizeSimd()
{
    static_assert(kernelType == ClusterDistanceKernelType::CpuSimd_4xM
                  || kernelType == ClusterDistanceKernelType::CpuSimd_2xMM);

#if GMX_SIMD && GMX_SIMD_HAVE_REAL
    if constexpr (kernelType == ClusterDistanceKernelType::CpuSimd_4xM)
    {
        return GMX_SIMD_REAL_WIDTH;
    }
    else
    {
        return GMX_SIMD_REAL_WIDTH / 2;
    }
#else
    return 0;
#endif
}

//! Stride of the packed x coordinate array
template<ClusterDistanceKernelType kernelType>
gmx_unused static constexpr int sc_xStride()
{
    return std::max(sc_iClusterSizeSimd<kernelType>(), sc_jClusterSizeSimd<kernelType>());
}

//! Returns the nbnxn coordinate data index given the i-cluster index
template<ClusterDistanceKernelType kernelType>
gmx_unused static inline int xIndexFromCi(int ci)
{
    constexpr int c_iClusterSize = sc_iClusterSizeSimd<kernelType>();
    constexpr int c_jClusterSize = sc_jClusterSizeSimd<kernelType>();

    static_assert(c_jClusterSize == c_iClusterSize / 2 || c_jClusterSize == c_iClusterSize
                          || c_jClusterSize == c_iClusterSize * 2,
                  "Only j-cluster sizes 2, 4 and 8 are currently implemented");

    constexpr int c_stride = sc_xStride<kernelType>();

    if constexpr (c_iClusterSize == c_stride)
    {
        return ci * DIM * c_stride;
    }
    else
    {
        static_assert(2 * c_iClusterSize == c_stride);

        // i-clusters are half the stride
        return (ci >> 1) * DIM * c_stride + (ci & 1) * (c_stride >> 1);
    }
}

//! Returns the nbnxn coordinate data index given the j-cluster index
template<ClusterDistanceKernelType kernelType>
gmx_unused static inline int xIndexFromCj(int cj)
{
    constexpr int c_iClusterSize = sc_iClusterSizeSimd<kernelType>();
    constexpr int c_jClusterSize = sc_jClusterSizeSimd<kernelType>();

    static_assert(c_jClusterSize == c_iClusterSize / 2 || c_jClusterSize == c_iClusterSize
                          || c_jClusterSize == c_iClusterSize * 2,
                  "Only j-cluster sizes 2, 4 and 8 are currently implemented");

    constexpr int c_stride = sc_xStride<kernelType>();

    if constexpr (c_jClusterSize == c_stride)
    {
        return cj * DIM * c_stride;
    }
    else
    {
        static_assert(2 * c_jClusterSize == c_stride);

        // j-clusters are half the stride
        return (cj >> 1) * DIM * c_stride + (cj & 1) * (c_stride >> 1);
    }
}

/*! \brief Returns the j-cluster index given the i-cluster index.
 *
 * \tparam    kernelType        The kernel type
 * \tparam    jSubClusterIndex  The j-sub-cluster index (0/1), used when size(j-cluster) <
 *                              size(i-cluster)
 * \param[in] ci                The i-cluster index
 */
template<ClusterDistanceKernelType kernelType, int jSubClusterIndex>
gmx_unused static inline int cjFromCi(int ci)
{
    constexpr int c_iClusterSize = sc_iClusterSizeSimd<kernelType>();
    constexpr int c_jClusterSize = sc_jClusterSizeSimd<kernelType>();

    static_assert(c_jClusterSize == c_iClusterSize / 2 || c_jClusterSize == c_iClusterSize
                          || c_jClusterSize == c_iClusterSize * 2,
                  "Only j-cluster sizes 2, 4 and 8 are currently implemented");

    static_assert(jSubClusterIndex == 0 || jSubClusterIndex == 1,
                  "Only sub-cluster indices 0 and 1 are supported");

    if constexpr (c_jClusterSize == c_iClusterSize / 2)
    {
        if (jSubClusterIndex == 0)
        {
            return ci << 1;
        }
        else
        {
            return ((ci + 1) << 1) - 1;
        }
    }
    else if constexpr (c_jClusterSize == c_iClusterSize)
    {
        return ci;
    }
    else
    {
        return ci >> 1;
    }
}

#if GMX_SIMD && GMX_USE_SIMD_KERNELS
//! Copies PBC shifted i-cell packed atom coordinates to working array
template<ClusterDistanceKernelType kernelType>
static inline void setICellCoordinatesSimd(int                   ci,
                                           const RVec&           shift,
                                           int gmx_unused        stride,
                                           const real*           x,
                                           NbnxmPairlistCpuWork* work)
{
    using namespace gmx;

    constexpr int c_iClusterSize = sc_iClusterSizeSimd<kernelType>();
    constexpr int c_jClusterSize = sc_jClusterSizeSimd<kernelType>();

    static_assert(c_jClusterSize == GMX_SIMD_REAL_WIDTH || 2 * c_jClusterSize == GMX_SIMD_REAL_WIDTH);

    constexpr int c_xStride = sc_xStride<kernelType>();

    real* x_ci_simd = work->iClusterData.xSimd.data();

    const int ia = xIndexFromCi<kernelType>(ci);

    if constexpr (c_jClusterSize == GMX_SIMD_REAL_WIDTH)
    {
        for (int i = 0; i < c_iClusterSize; i++)
        {
            for (int d = 0; d < DIM; d++)
            {
                store(x_ci_simd + (i * DIM + d) * GMX_SIMD_REAL_WIDTH,
                      SimdReal(x[ia + d * c_xStride + i] + shift[d]));
            }
        }
    }
    else
    {
#    if GMX_SIMD_HAVE_HSIMD_UTIL_REAL
        for (int i = 0; i < c_iClusterSize / 2; i++)
        {
            for (int d = 0; d < DIM; d++)
            {
                store(x_ci_simd + (i * DIM + d) * GMX_SIMD_REAL_WIDTH,
                      loadU1DualHsimd(x + ia + d * c_xStride + i * 2) + SimdReal(shift[d]));
            }
        }
#    else
        GMX_RELEASE_ASSERT(false, "Function called that is not supported in this build");
#    endif
    }
}
#endif // GMX_SIMD && GMX_USE_SIMD_KERNELS

void setICellCoordinatesSimd4xM(int gmx_unused ci,
                                const RVec gmx_unused& shift,
                                int gmx_unused         stride,
                                const real gmx_unused* x,
                                NbnxmPairlistCpuWork gmx_unused* work)
{
#if GMX_HAVE_NBNXM_SIMD_4XM
    setICellCoordinatesSimd<ClusterDistanceKernelType::CpuSimd_4xM>(ci, shift, stride, x, work);
#else
    GMX_RELEASE_ASSERT(false, "Function called that is not supported in this build");
#endif
}

void setICellCoordinatesSimd2xMM(int gmx_unused ci,
                                 const RVec gmx_unused& shift,
                                 int gmx_unused         stride,
                                 const real gmx_unused* x,
                                 NbnxmPairlistCpuWork gmx_unused* work)
{
#if GMX_HAVE_NBNXM_SIMD_2XMM
    setICellCoordinatesSimd<ClusterDistanceKernelType::CpuSimd_2xMM>(ci, shift, stride, x, work);
#else
    GMX_RELEASE_ASSERT(false, "Function called that is not supported in this build");
#endif
}

#if GMX_SIMD && GMX_USE_SIMD_KERNELS

template<ClusterDistanceKernelType kernelType>
static inline SimdReal loadJData(const real* x)
{
    constexpr int c_jClusterSize = sc_jClusterSizeSimd<kernelType>();

    static_assert(c_jClusterSize == GMX_SIMD_REAL_WIDTH || 2 * c_jClusterSize == GMX_SIMD_REAL_WIDTH);

    if constexpr (c_jClusterSize == GMX_SIMD_REAL_WIDTH)
    {
        return load<SimdReal>(x);
    }
    else
    {
#    if GMX_SIMD_HAVE_HSIMD_UTIL_REAL
        return loadDuplicateHsimd(x);
#    else
        GMX_RELEASE_ASSERT(false, "Function called that is not supported in this build");
        return setZero();
#    endif
    }
}

/*! \brief SIMD code for checking and adding cluster-pairs to the list using coordinates in packed format.
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
template<ClusterDistanceKernelType kernelType>
static inline void makeClusterListSimd(const Grid&              jGrid,
                                       NbnxnPairlistCpu*        nbl,
                                       int                      icluster,
                                       int                      firstCell,
                                       int                      lastCell,
                                       bool                     excludeSubDiagonal,
                                       const real* gmx_restrict x_j,
                                       real                     rlist2,
                                       float                    rbb2,
                                       int* gmx_restrict        numDistanceChecks)
{
    using namespace gmx;

    constexpr int c_iClusterSize = sc_iClusterSizeSimd<kernelType>();
    constexpr int c_jClusterSize = sc_jClusterSizeSimd<kernelType>();

    constexpr int c_numIRegisters = (c_iClusterSize * c_jClusterSize) / GMX_SIMD_REAL_WIDTH;

    const real* gmx_restrict        x_ci_simd = nbl->work->iClusterData.xSimd.data();
    const BoundingBox* gmx_restrict bb_ci     = nbl->work->iClusterData.bb.data();

    /* Convert the j-range from i-cluster size indexing to j-cluster indexing */
    int jclusterFirst = cjFromCi<kernelType, 0>(firstCell);
    int jclusterLast  = cjFromCi<kernelType, 1>(lastCell);
    GMX_ASSERT(jclusterLast >= jclusterFirst,
               "We should have a non-empty j-cluster range, since the calling code should have "
               "ensured a non-empty cell range");

    const SimdReal rc2_S = SimdReal(rlist2);

    bool InRange = false;
    while (!InRange && jclusterFirst <= jclusterLast)
    {
        const float d2 = clusterBoundingBoxDistance2(bb_ci[0], jGrid.jBoundingBoxes()[jclusterFirst]);
        *numDistanceChecks += 2;

        /* Check if the distance is within the distance where
         * we use only the bounding box distance rbb,
         * or within the cut-off and there is at least one atom pair
         * within the cut-off.
         */
        if (d2 < rbb2)
        {
            InRange = true;
        }
        else if (d2 < rlist2)
        {
            const int xind_f =
                    xIndexFromCj<kernelType>(cjFromCi<kernelType, 0>(jGrid.cellOffset()) + jclusterFirst);

            const SimdReal jx_S = loadJData<kernelType>(x_j + xind_f + 0 * sc_xStride<kernelType>());
            const SimdReal jy_S = loadJData<kernelType>(x_j + xind_f + 1 * sc_xStride<kernelType>());
            const SimdReal jz_S = loadJData<kernelType>(x_j + xind_f + 2 * sc_xStride<kernelType>());


            /* Calculate distance */
            std::array<std::array<SimdReal, 3>, c_numIRegisters> d;
            for (int i = 0; i < c_numIRegisters; i++)
            {
                d[i][0] = load<SimdReal>(x_ci_simd + (i * 3 + 0) * GMX_SIMD_REAL_WIDTH) - jx_S;
                d[i][1] = load<SimdReal>(x_ci_simd + (i * 3 + 1) * GMX_SIMD_REAL_WIDTH) - jy_S;
                d[i][2] = load<SimdReal>(x_ci_simd + (i * 3 + 2) * GMX_SIMD_REAL_WIDTH) - jz_S;
            }

            /* rsq = dx*dx+dy*dy+dz*dz */
            std::array<SimdReal, c_numIRegisters> rsq;
            for (int i = 0; i < c_numIRegisters; i++)
            {
                rsq[i] = norm2(d[i][0], d[i][1], d[i][2]);
            }

            std::array<SimdBool, c_numIRegisters> wco;
            for (int i = 0; i < c_numIRegisters; i++)
            {
                wco[i] = (rsq[i] < rc2_S);
            }

            constexpr int numBitShifts = StaticLog2<c_numIRegisters>::value;
            for (int bitShift = 0; bitShift < numBitShifts; bitShift++)
            {
                const int offset = (1 << bitShift);
                for (int i = 0; i < c_numIRegisters; i += 2 * offset)
                {
                    wco[i] = wco[i] || wco[i + offset];
                }
            }

            InRange = anyTrue(wco[0]);

            *numDistanceChecks += c_numIRegisters * GMX_SIMD_REAL_WIDTH;
        }
        if (!InRange)
        {
            jclusterFirst++;
        }
    }
    if (!InRange)
    {
        return;
    }

    InRange = false;
    while (!InRange && jclusterLast > jclusterFirst)
    {
        const float d2 = clusterBoundingBoxDistance2(bb_ci[0], jGrid.jBoundingBoxes()[jclusterLast]);
        *numDistanceChecks += 2;

        /* Check if the distance is within the distance where
         * we use only the bounding box distance rbb,
         * or within the cut-off and there is at least one atom pair
         * within the cut-off.
         */
        if (d2 < rbb2)
        {
            InRange = true;
        }
        else if (d2 < rlist2)
        {
            const int xind_l =
                    xIndexFromCj<kernelType>(cjFromCi<kernelType, 0>(jGrid.cellOffset()) + jclusterLast);

            const SimdReal jx_S = loadJData<kernelType>(x_j + xind_l + 0 * sc_xStride<kernelType>());
            const SimdReal jy_S = loadJData<kernelType>(x_j + xind_l + 1 * sc_xStride<kernelType>());
            const SimdReal jz_S = loadJData<kernelType>(x_j + xind_l + 2 * sc_xStride<kernelType>());

            /* Calculate distance */
            std::array<std::array<SimdReal, 3>, c_numIRegisters> d;
            for (int i = 0; i < c_numIRegisters; i++)
            {
                d[i][0] = load<SimdReal>(x_ci_simd + (i * 3 + 0) * GMX_SIMD_REAL_WIDTH) - jx_S;
                d[i][1] = load<SimdReal>(x_ci_simd + (i * 3 + 1) * GMX_SIMD_REAL_WIDTH) - jy_S;
                d[i][2] = load<SimdReal>(x_ci_simd + (i * 3 + 2) * GMX_SIMD_REAL_WIDTH) - jz_S;
            }

            /* rsq = dx*dx+dy*dy+dz*dz */
            std::array<SimdReal, c_numIRegisters> rsq;
            for (int i = 0; i < c_numIRegisters; i++)
            {
                rsq[i] = norm2(d[i][0], d[i][1], d[i][2]);
            }

            std::array<SimdBool, c_numIRegisters> wco;
            for (int i = 0; i < c_numIRegisters; i++)
            {
                wco[i] = (rsq[i] < rc2_S);
            }

            constexpr int numBitShifts = StaticLog2<c_numIRegisters>::value;
            for (int bitShift = 0; bitShift < numBitShifts; bitShift++)
            {
                const int offset = (1 << bitShift);
                for (int i = 0; i < c_numIRegisters; i += 2 * offset)
                {
                    wco[i] = wco[i] || wco[i + offset];
                }
            }

            InRange = anyTrue(wco[0]);

            *numDistanceChecks += c_numIRegisters * GMX_SIMD_REAL_WIDTH;
        }
        if (!InRange)
        {
            jclusterLast--;
        }
    }

    if (jclusterFirst <= jclusterLast)
    {
        for (int jcluster = jclusterFirst; jcluster <= jclusterLast; jcluster++)
        {
            /* Store cj and the interaction mask */
            nbnxn_cj_t cjEntry;
            cjEntry.cj = cjFromCi<kernelType, 0>(jGrid.cellOffset()) + jcluster;
            cjEntry.excl =
                    getImask<c_iClusterSize, c_jClusterSize>(excludeSubDiagonal, icluster, jcluster);
            nbl->cj.push_back(cjEntry);
        }
        /* Increase the closing index in the i list */
        nbl->ci.back().cj_ind_end = nbl->cj.size();
    }
}

#endif // GMX_SIMD && GMX_USE_SIMD_KERNELS

void makeClusterListSimd4xM(const Grid gmx_unused& jGrid,
                            NbnxnPairlistCpu gmx_unused* nbl,
                            int gmx_unused               icluster,
                            int gmx_unused               firstCell,
                            int gmx_unused               lastCell,
                            bool gmx_unused              excludeSubDiagonal,
                            const real gmx_unused* gmx_restrict x_j,
                            real gmx_unused                     rlist2,
                            float gmx_unused                    rbb2,
                            int gmx_unused* gmx_restrict numDistanceChecks)
{
#if GMX_HAVE_NBNXM_SIMD_4XM
    makeClusterListSimd<ClusterDistanceKernelType::CpuSimd_4xM>(
            jGrid, nbl, icluster, firstCell, lastCell, excludeSubDiagonal, x_j, rlist2, rbb2, numDistanceChecks);
#else
    GMX_RELEASE_ASSERT(false, "Function called that is not supported in this build");
#endif
}

void makeClusterListSimd2xMM(const Grid gmx_unused& jGrid,
                             NbnxnPairlistCpu gmx_unused* nbl,
                             int gmx_unused               icluster,
                             int gmx_unused               firstCell,
                             int gmx_unused               lastCell,
                             bool gmx_unused              excludeSubDiagonal,
                             const real gmx_unused* gmx_restrict x_j,
                             real gmx_unused                     rlist2,
                             float gmx_unused                    rbb2,
                             int gmx_unused* gmx_restrict numDistanceChecks)
{
#if GMX_HAVE_NBNXM_SIMD_2XMM
    makeClusterListSimd<ClusterDistanceKernelType::CpuSimd_2xMM>(
            jGrid, nbl, icluster, firstCell, lastCell, excludeSubDiagonal, x_j, rlist2, rbb2, numDistanceChecks);
#else
    GMX_RELEASE_ASSERT(false, "Function called that is not supported in this build");
#endif
}

} // namespace gmx
