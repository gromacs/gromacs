/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 *  \brief
 *  NBNXM SYCL kernels
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "nbnxm_sycl_kernel_pruneonly.h"

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/utility/template_mp.h"

#include "nbnxm_sycl_kernel_utils.h"
#include "nbnxm_sycl_types.h"

using sycl::access::fence_space;
using mode = sycl::access_mode;

//! \brief Class name for NBNXM prune-only kernel
template<bool haveFreshList>
class NbnxmKernelPruneOnly;

namespace Nbnxm
{

/*! \brief Prune-only kernel for NBNXM.
 *
 */
template<bool haveFreshList>
auto nbnxmKernelPruneOnly(sycl::handler&                                cgh,
                          DeviceAccessor<Float4, mode::read>            a_xq,
                          DeviceAccessor<Float3, mode::read>            a_shiftVec,
                          DeviceAccessor<nbnxn_cj4_t, mode::read_write> a_plistCJ4,
                          DeviceAccessor<nbnxn_sci_t, mode::read>       a_plistSci,
                          DeviceAccessor<unsigned int, haveFreshList ? mode::write : mode::read> a_plistIMask,
                          const float rlistOuterSq,
                          const float rlistInnerSq,
                          const int   numParts,
                          const int   part)
{
    a_xq.bind(cgh);
    a_shiftVec.bind(cgh);
    a_plistCJ4.bind(cgh);
    a_plistSci.bind(cgh);
    a_plistIMask.bind(cgh);

    /* shmem buffer for i x+q pre-loading */
    sycl_2020::local_accessor<Float4, 1> sm_xq(
            sycl::range<1>(c_nbnxnGpuNumClusterPerSupercluster * c_clSize), cgh);

    constexpr int warpSize = c_clSize * c_clSize / 2;

    /* Somewhat weird behavior inherited from OpenCL.
     * With clSize == 4, we use sub_group size of 16 (not enforced in OpenCL implementation, but chosen
     * by the IGC compiler), however for data layout we consider it to be 8.
     * Setting sub_group size to 8 slows down the prune-only kernel 1.5-2 times.
     * For clSize == But we need to set specific sub_group size >= 32 for clSize == 8 for correctness,
     * but it causes very poor performance.
     */
    constexpr int gmx_unused requiredSubGroupSize = (c_clSize == 4) ? 16 : warpSize;

    /* Requirements:
     * Work group (block) must have range (c_clSize, c_clSize, ...) (for itemIdx calculation, easy
     * to change). */
    return [=](sycl::nd_item<3> itemIdx) [[intel::reqd_sub_group_size(requiredSubGroupSize)]]
    {
        // thread/block/warp id-s
        const unsigned tidxi = itemIdx.get_local_id(2);
        const unsigned tidxj = itemIdx.get_local_id(1);
        const int      tidx  = tidxj * c_clSize + tidxi;
        const unsigned tidxz = itemIdx.get_local_id(0);
        const unsigned bidx  = itemIdx.get_group(0);

        const sycl::sub_group sg   = itemIdx.get_sub_group();
        const unsigned        widx = tidx / warpSize;

        // my i super-cluster's index = sciOffset + current bidx * numParts + part
        const nbnxn_sci_t nbSci     = a_plistSci[bidx * numParts + part];
        const int         sci       = nbSci.sci;           /* super-cluster */
        const int         cij4Start = nbSci.cj4_ind_start; /* first ...*/
        const int         cij4End   = nbSci.cj4_ind_end;   /* and last index of j clusters */

        // We may need only a subset of threads active for preloading i-atoms
        // depending on the super-cluster and cluster / thread-block size.
        constexpr bool c_loadUsingAllXYThreads = (c_clSize == c_nbnxnGpuNumClusterPerSupercluster);
        if (tidxz == 0 && (c_loadUsingAllXYThreads || tidxj < c_nbnxnGpuNumClusterPerSupercluster))
        {
            for (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i += c_clSize)
            {
                /* Pre-load i-atom x and q into shared memory */
                const int ci = sci * c_nbnxnGpuNumClusterPerSupercluster + tidxj + i;
                const int ai = ci * c_clSize + tidxi;

                /* We don't need q, but using float4 in shmem avoids bank conflicts.
                   (but it also wastes L2 bandwidth). */
                const Float4 xq    = a_xq[ai];
                const Float3 shift = a_shiftVec[nbSci.shift];
                const Float4 xi(xq[0] + shift[0], xq[1] + shift[1], xq[2] + shift[2], xq[3]);
                sm_xq[(tidxj + i) * c_clSize + tidxi] = xi;
            }
        }
        itemIdx.barrier(fence_space::local_space);

        /* loop over the j clusters = seen by any of the atoms in the current super-cluster.
         * The loop stride c_syclPruneKernelJ4Concurrency ensures that consecutive warps-pairs are
         * assigned consecutive j4's entries. */
        for (int j4 = cij4Start + tidxz; j4 < cij4End; j4 += c_syclPruneKernelJ4Concurrency)
        {
            unsigned imaskFull, imaskCheck, imaskNew;

            if constexpr (haveFreshList)
            {
                /* Read the mask from the list transferred from the CPU */
                imaskFull = a_plistCJ4[j4].imei[widx].imask;
                /* We attempt to prune all pairs present in the original list */
                imaskCheck = imaskFull;
                imaskNew   = 0;
            }
            else
            {
                /* Read the mask from the "warp-pruned" by rlistOuter mask array */
                imaskFull = a_plistIMask[j4 * c_nbnxnGpuClusterpairSplit + widx];
                /* Read the old rolling pruned mask, use as a base for new */
                imaskNew = a_plistCJ4[j4].imei[widx].imask;
                /* We only need to check pairs with different mask */
                imaskCheck = (imaskNew ^ imaskFull);
            }

            if (imaskCheck)
            {
                for (int jm = 0; jm < c_nbnxnGpuJgroupSize; jm++)
                {
                    if (imaskCheck & (superClInteractionMask << (jm * c_nbnxnGpuNumClusterPerSupercluster)))
                    {
                        unsigned mask_ji = (1U << (jm * c_nbnxnGpuNumClusterPerSupercluster));
                        // SYCL-TODO: Reevaluate prefetching methods
                        const int cj = a_plistCJ4[j4].cj[jm];
                        const int aj = cj * c_clSize + tidxj;

                        /* load j atom data */
                        const Float4 tmp = a_xq[aj];
                        const Float3 xj(tmp[0], tmp[1], tmp[2]);

                        for (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
                        {
                            if (imaskCheck & mask_ji)
                            {
                                // load i-cluster coordinates from shmem
                                const Float4 xi = sm_xq[i * c_clSize + tidxi];
                                // distance between i and j atoms
                                Float3 rv(xi[0], xi[1], xi[2]);
                                rv -= xj;
                                const float r2 = norm2(rv);

                                /* If _none_ of the atoms pairs are in rlistOuter
                                 * range, the bit corresponding to the current
                                 * cluster-pair in imask gets set to 0. */
                                if (haveFreshList && !(sycl::any_of_group(sg, r2 < rlistOuterSq)))
                                {
                                    imaskFull &= ~mask_ji;
                                }
                                /* If any atom pair is within range, set the bit
                                 * corresponding to the current cluster-pair. */
                                if (sycl::any_of_group(sg, r2 < rlistInnerSq))
                                {
                                    imaskNew |= mask_ji;
                                }
                            } // (imaskCheck & mask_ji)
                            /* shift the mask bit by 1 */
                            mask_ji += mask_ji;
                        } // (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
                    } // (imaskCheck & (superClInteractionMask << (jm * c_nbnxnGpuNumClusterPerSupercluster)))
                } // for (int jm = 0; jm < c_nbnxnGpuJgroupSize; jm++)

                if constexpr (haveFreshList)
                {
                    /* copy the list pruned to rlistOuter to a separate buffer */
                    a_plistIMask[j4 * c_nbnxnGpuClusterpairSplit + widx] = imaskFull;
                }
                /* update the imask with only the pairs up to rlistInner */
                a_plistCJ4[j4].imei[widx].imask = imaskNew;
            } // (imaskCheck)
        } // for (int j4 = cij4_start + tidxz; j4 < cij4_end; j4 += c_syclPruneKernelJ4Concurrency)
    };
}

//! \brief Leap Frog SYCL prune-only kernel launch code.
template<bool haveFreshList, class... Args>
sycl::event launchNbnxmKernelPruneOnly(const DeviceStream& deviceStream, const int numSciInPart, Args&&... args)
{
    using kernelNameType = NbnxmKernelPruneOnly<haveFreshList>;

    /* Kernel launch config:
     * - The thread block dimensions match the size of i-clusters, j-clusters,
     *   and j-cluster concurrency, in x, y, and z, respectively.
     * - The 1D block-grid contains as many blocks as super-clusters.
     */
    const unsigned long     numBlocks = numSciInPart;
    const sycl::range<3>    blockSize{ c_syclPruneKernelJ4Concurrency, c_clSize, c_clSize };
    const sycl::range<3>    globalSize{ numBlocks * blockSize[0], blockSize[1], blockSize[2] };
    const sycl::nd_range<3> range{ globalSize, blockSize };

    sycl::queue q = deviceStream.stream();

    sycl::event e = q.submit([&](sycl::handler& cgh) {
        auto kernel = nbnxmKernelPruneOnly<haveFreshList>(cgh, std::forward<Args>(args)...);
        cgh.parallel_for<kernelNameType>(range, kernel);
    });

    return e;
}

//! \brief Select templated kernel and launch it.
template<class... Args>
sycl::event chooseAndLaunchNbnxmKernelPruneOnly(bool haveFreshList, Args&&... args)
{
    return gmx::dispatchTemplatedFunction(
            [&](auto haveFreshList_) {
                return launchNbnxmKernelPruneOnly<haveFreshList_>(std::forward<Args>(args)...);
            },
            haveFreshList);
}

void launchNbnxmKernelPruneOnly(NbnxmGpu*                 nb,
                                const InteractionLocality iloc,
                                const int                 numParts,
                                const int                 part,
                                const int                 numSciInPart)
{
    NBAtomDataGpu*      adat          = nb->atdat;
    NBParamGpu*         nbp           = nb->nbparam;
    gpu_plist*          plist         = nb->plist[iloc];
    const bool          haveFreshList = plist->haveFreshList;
    const DeviceStream& deviceStream  = *nb->deviceStreams[iloc];

    sycl::event e = chooseAndLaunchNbnxmKernelPruneOnly(haveFreshList,
                                                        deviceStream,
                                                        numSciInPart,
                                                        adat->xq,
                                                        adat->shiftVec,
                                                        plist->cj4,
                                                        plist->sci,
                                                        plist->imask,
                                                        nbp->rlistOuter_sq,
                                                        nbp->rlistInner_sq,
                                                        numParts,
                                                        part);
}

} // namespace Nbnxm
