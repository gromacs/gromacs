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
 *  Data management and kernel launch functions for nbnxm sycl.
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "gromacs/nbnxm/gpu_common.h"
#include "gromacs/utility/exceptions.h"

#include "nbnxm_sycl_kernel.h"
#include "nbnxm_sycl_kernel_pruneonly.h"
#include "nbnxm_sycl_types.h"

namespace gmx
{

static void launchSciSortOnGpu(GpuPairlist* plist, const DeviceStream& deviceStream);

void gpu_launch_kernel_pruneonly(NbnxmGpu* nb, const InteractionLocality iloc, const int numParts)
{
    auto* plist = nb->plist[iloc].get();

    if (plist->haveFreshList)
    {
        GMX_ASSERT(numParts == 1, "With first pruning we expect 1 part");

        /* Set rollingPruningNumParts to signal that it is not set */
        plist->rollingPruningNumParts = 0;
    }
    else
    {
        if (plist->rollingPruningNumParts == 0)
        {
            plist->rollingPruningNumParts = numParts;
        }
        else
        {
            GMX_ASSERT(numParts == plist->rollingPruningNumParts,
                       "It is not allowed to change numParts in between list generation steps");
        }
    }

    /* Compute the max number of list entries to prune across all passes
     * Note that the actual number for a specific pass will be computed inside the kernel.
     * Also note that this implementation (parts tracking on device) differs from some
     * other backends (parts tracking on host, passed as kernel argument).
     */

    /* Compute the max number of list entries to prune in this pass */
    const int numSciInPartMax = (plist->numSci) / numParts;

    /* Don't launch the kernel if there is no work to do */
    if (numSciInPartMax <= 0)
    {
        plist->haveFreshList = false;
        return;
    }

    launchNbnxmKernelPruneOnly(nb, iloc, numParts, numSciInPartMax);
    if (plist->haveFreshList && nbnxmSortListsOnGpu())
    {
        launchSciSortOnGpu(plist, *nb->deviceStreams[iloc]);
    }

    if (plist->haveFreshList)
    {
        plist->haveFreshList = false;
        nb->didPrune[iloc]   = true; // Mark that pruning has been done
    }
    else
    {
        nb->didRollingPrune[iloc] = true; // Mark that rolling pruning has been done
    }
}


void gpu_launch_kernel(NbnxmGpu* nb, const gmx::StepWorkload& stepWork, const InteractionLocality iloc)
{
    const NBParamGpu* nbp   = nb->nbparam;
    auto*             plist = nb->plist[iloc].get();

    if (canSkipNonbondedWork(*nb, iloc))
    {
        plist->haveFreshList = false;
        return;
    }

    if (nbp->useDynamicPruning && plist->haveFreshList)
    {
        // Prunes for rlistOuter and rlistInner, sets plist->haveFreshList=false
        gpu_launch_kernel_pruneonly(nb, iloc, 1);
    }

    if (plist->numSci == 0)
    {
        /* Don't launch an empty local kernel */
        return;
    }

    /* Whether we need to call a combined prune and interaction kernel or just an interaction
     * kernel. doPrune being true implies we are not using dynamic pruning and are in the first
     * call to the interaction kernel after a neighbour list step */
    bool doPrune = (plist->haveFreshList && !nb->timers->interaction[iloc].didPrune);

    launchNbnxmKernel(nb, stepWork, iloc, doPrune);
    if (doPrune && nbnxmSortListsOnGpu())
    {
        launchSciSortOnGpu(plist, *nb->deviceStreams[iloc]);
    }
}

/*! \brief SYCL exclusive prefix sum kernel for list sorting.
 *
 * As of oneAPI 2024.1, \c oneapi::dpl::experimental::exclusive_scan_async for inputs <= 16384
 * elements simply launches a single work-group and uses sycl::joint_exclusive_scan. We have,
 * somewhat arbitrary, input size of 8192, so we're fine replicating the same approach.
 *
 * NVIDIA's CUB uses fancier approach ("Single-pass Parallel Prefix Scan with Decoupled Look-back"),
 * but we are unlikely to need it here since this kernel is very small anyway.
 *
 * \tparam workGroupSize Size of the (only) work-group.
 * \tparam nElements Input array size; should be a multiple of workGroupSize.
 * \param gm_input Input data buffer, should contain nElements elements of type \c int.
 * \param gm_output Output data buffer, should have enough space for nElements elements of type \c int.
 *
 * \warning This kernel should be launched with only a single work-group of size workGroupSize.
 * \warning This kernel is inefficient for large inputs (oneDPL uses it with <= 16384 elements).
 */
template<int workGroupSize, int nElements>
static auto nbnxnKernelExclusivePrefixSum(const int* __restrict__ gm_input, int* __restrict__ gm_output)
{
    static_assert(nElements % workGroupSize == 0, "This simple scan kernel does not handle padding");
    return [=](sycl::nd_item<1> itemIdx) {
        const sycl::group<1> workGroup = itemIdx.get_group();
        sycl::joint_exclusive_scan(
                workGroup, gm_input, gm_input + nElements, gm_output, 0, sycl::plus<int>{});
    };
}
//! SYCL kernel name
class ExclusivePrefixSum;

/*! \brief SYCL bucket sci sort kernel.
 *
 *  Sorts sci in order from most to least neighbours, using the count sort algorithm
 *
 *  Unlike the cpu version of sci sort, this kernel uses counts which only contain pairs which have
 *  not been masked out, giving an ordering which more accurately represents the work which will be
 *  done in the non bonded force kernel. The counts themselves are generated in the prune kernel.
 *
 * \param gm_sci Unsorted pair list.
 * \param gm_sciCount Total number of sci with exactly \c i neighbours
 * \param gm_sciOffset Exclusive prefix sum of gm_sciCount. \c gm_sciOffset[i] is the offset that
 *                     the first sci with \c i neighbours will have in the sorted sci
 *                     list. All other sci with i neighbours will be placed randomly in
 *                     positions \c gm_sciOffset[i] to \c gm_sciOffset[i+1] exclusive.
 * \param gm_sciSorted Sorted pair list.
 */
static auto nbnxnKernelBucketSciSort(const nbnxn_sci_t* __restrict__ gm_sci,
                                     const int* __restrict__ gm_sciCount,
                                     int* __restrict__ gm_sciOffset,
                                     nbnxn_sci_t* __restrict__ gm_sciSorted)
{
    return [=](sycl::id<1> itemIdx) {
        using sycl::memory_order, sycl::memory_scope, sycl::access::address_space;

        const int         idx      = itemIdx[0];
        const nbnxn_sci_t sci      = gm_sci[idx];
        const int         sciCount = gm_sciCount[idx];
        // Choose an order in the sorted list for sci with the same number of neighbours as each other.
        // We want each sci with the same count to go somewhere in the list between sciOffset[count]
        // to sciOffset[count+1] exclusive.
        // As the amount of work for each of these in the non bonded force kernel will
        // be equivalent, the order within these bounds does not matter.
        sycl::atomic_ref<int, memory_order::relaxed, memory_scope::device, address_space::global_space> offsetRef(
                gm_sciOffset[sciCount]);

        const int sciOffset = offsetRef.fetch_add(1);
        // Insert the sci into the sorted list at the chosen index
        gm_sciSorted[sciOffset] = sci;
    };
}
//! SYCL kernel name
class BucketSciSort;

template<int workGroupSize>
static void launchPrefixSumKernel(sycl::queue& q, GpuPairlistSorting* sorting)
{
    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        cgh.parallel_for<ExclusivePrefixSum>(
                sycl::nd_range<1>{ workGroupSize, workGroupSize },
                nbnxnKernelExclusivePrefixSum<workGroupSize, c_sciHistogramSize>(
                        sorting->sciHistogram.get_pointer(), sorting->sciOffset.get_pointer()));
    });
}

static void launchBucketSortKernel(sycl::queue& q, GpuPairlist* plist)
{
    const size_t size = plist->numSci;
    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        cgh.parallel_for<BucketSciSort>(
                sycl::range<1>{ size },
                nbnxnKernelBucketSciSort(plist->sci.get_pointer(),
                                         plist->sorting.sciCount.get_pointer(),
                                         plist->sorting.sciOffset.get_pointer(),
                                         plist->sorting.sciSorted.get_pointer()));
    });
}
static void launchSciSortOnGpu(GpuPairlist* plist, const DeviceStream& deviceStream)
{
    sycl::queue q = deviceStream.stream();

    /* We are launching a single work-group, and it should be, in principle, as large as possible.
     * E.g., on PVC1100, wgSizeScan=1024 is ~1.2 times faster than wgSizeScan=512, and on
     * MI250X ~1.7 times faster. But Intel iGPUs only handle 512.
     * It should be autodetected, but for now we just use the universally supported value.
     * The kernel takes 10-100µs and is run rarely, so it is of minor concern. */
    constexpr int wgSizeScan = 512;
    launchPrefixSumKernel<wgSizeScan>(q, &plist->sorting);
    launchBucketSortKernel(q, plist);
}

} // namespace gmx
