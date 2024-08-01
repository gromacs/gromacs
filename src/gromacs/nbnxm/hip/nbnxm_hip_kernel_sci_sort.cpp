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
 *  \brief
 *  Explicitly instantiate NBNXM HIP kernel for sci sorting
 *
 *  \ingroup module_nbnxm
 */

#include "gromacs/nbnxm/hip/nbnxm_hip_kernel_sci_sort.h"

#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/gpu_types_common.h"

#include "nbnxm_hip_types.h"

namespace gmx
{

namespace
{

/*! \brief HIP bucket sci sort kernel.
 *
 *  Sorts sci in order from most to least neighbours, using the count sort algorithm
 *
 *  Unlike the cpu version of sci sort, this kernel uses counts which only contain pairs which have
 *  not been masked out, giving an ordering which more accurately represents the work which will be
 *  done in the non bonded force kernel. The counts themselves are generated in the prune kernel.
 *
 *  Inputs:
 *   - plist.sci = unsorted pair list
 *   - plist.sorting.sciCount[i] = total number of sci with exactly i neighbours
 *   - plist.sorting.sciOffset = exclusive prefix sum of sciCount.
 *     plist.sorting.sciOffset[i] represents the offset that the first sci with i neighbours will
 *     have in the sorted sci list. All other sci with i neighbours will be placed randomly in
 * positions plist.sorting.sciOffset[i] to plist.sorting.sciOffset[i+1] exclusive
 */
template<PairlistType pairlistType, int threadsPerBlock>
__launch_bounds__(threadsPerBlock) __global__
        void nbnxmKernelBucketSciSort(GpuPairlist<pairlistType> plist)
{
    int size = plist.numSci;

    const int tid         = threadIdx.x;
    const int blockOffset = blockIdx.x * threadsPerBlock;

    if (size > (blockOffset + tid))
    {
        nbnxn_sci_t sci      = plist.sci[blockOffset + tid];
        int         sciCount = plist.sorting.sciCount[blockOffset + tid];

        // Choose an order in the sorted list for sci with the same number of neighbours as each other.
        // We want each sci with the same count to go somewhere in the list between sciOffset[count]
        // to sciOffset[count+1] exclusive.
        // As the amount of work for each of these in the non bonded force kernel will
        // be equivalent, the order within these bounds does not matter.
        int sciOffset = atomicAdd(&plist.sorting.sciOffset[sciCount], 1);
        // Insert the sci into the sorted list at the chosen index
        plist.sorting.sciSorted[sciOffset] = sci;
    }
}

//! \brief NBNXM kernel launch code.
template<PairlistType pairlistType>
void launchNbnxmKernelHelperSciSort(const DeviceStream& deviceStream, GpuPairlist<pairlistType>* plist)
{
    performExclusiveScan(plist->sorting.nscanTemporary, plist->sorting.scanTemporary, plist, deviceStream);

    KernelLaunchConfig configSortSci;
    configSortSci.blockSize[0]     = c_sciSortingThreadsPerBlock;
    configSortSci.blockSize[1]     = 1;
    configSortSci.blockSize[2]     = 1;
    configSortSci.gridSize[0]      = divideRoundUp(plist->numSci, c_sciSortingThreadsPerBlock);
    configSortSci.sharedMemorySize = 0;
    const auto kernelSciSort = nbnxmKernelBucketSciSort<pairlistType, c_sciSortingThreadsPerBlock>;
    const auto kernelSciSortArgs = prepareGpuKernelArguments(kernelSciSort, configSortSci, plist);

    launchGpuKernel(kernelSciSort, configSortSci, deviceStream, nullptr, "nbnxn_kernel_sci_sort", kernelSciSortArgs);
}

} // namespace

void launchNbnxmKernelSciSort(NbnxmGpu* nb, InteractionLocality iloc)
{
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];
    std::visit(
            [&](auto&& pairlists)
            {
                auto* plist                 = pairlists[iloc].get();
                using T                     = std::decay_t<decltype(*plist)>;
                constexpr auto pairlistType = gmx::getPairlistTypeFromPairlist<T>();
                if constexpr (isGpuSpecificPairlist(pairlistType))
                {
                    launchNbnxmKernelHelperSciSort<pairlistType>(deviceStream, plist);
                }
            },
            nb->plist);
}

} // namespace gmx
