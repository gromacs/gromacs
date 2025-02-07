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
 *  NBNXM HIP prune only kernels
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "nbnxm_hip_kernel_pruneonly.h"

#include <type_traits>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/devicebuffer_hip.h"
#include "gromacs/gpu_utils/gpu_kernel_utils.h"
#include "gromacs/gpu_utils/hiputils.h"
#include "gromacs/nbnxm/gpu_types_common.h"
#include "gromacs/nbnxm/nbnxm_enums.h"
#include "gromacs/utility/template_mp.h"

#include "nbnxm_hip_kernel_utils.h"
#include "nbnxm_hip_types.h"

namespace gmx
{
//! Helper method to calculate launch bounds
template<bool hasLargeRegisterPool>
constexpr int minBlocksPp = hasLargeRegisterPool ? 8 : 1;

#if defined(__gfx90a__)
constexpr bool c_deviceLargeRegisterPool = true;
#else
constexpr bool c_deviceLargeRegisterPool = false;
#endif

/*! \brief Prune-only kernel for NBNXM.
 *
 * The number of threads per block is always c_clSizeSq, so we can use it here without having
 * to pass another template argument. The remaining multipliers are architecture and kernel
 * flavor specific and need to be passed in.
 *
 */
template<PairlistType pairlistType, bool hasLargeRegisterPool, bool haveFreshList, int threadZ>
__launch_bounds__(c_clSizeSq<pairlistType>* threadZ, minBlocksPp<hasLargeRegisterPool>) __global__
        static void nbnxmKernelPruneOnly(NBAtomDataGpu             atdat,
                                         NBParamGpu                nbparam,
                                         GpuPairlist<pairlistType> plist,
                                         int                       numParts)
{
    // we only want to actually compile the kernel for archs that have matching execution widths
    constexpr bool c_matchingExecutionWidth = sc_gpuParallelExecutionWidth(pairlistType) == warpSize;
    constexpr bool c_deviceRegisterPoolMatch = hasLargeRegisterPool == c_deviceLargeRegisterPool;
    constexpr bool c_doGenerateKernel = c_deviceRegisterPoolMatch && c_matchingExecutionWidth;
    if constexpr (c_doGenerateKernel)
    {
        constexpr int  c_clSize                 = sc_gpuClusterSize(pairlistType);
        constexpr int  c_clusterPerSuperCluster = sc_gpuClusterPerSuperCluster(pairlistType);
        constexpr int  c_gpuJGroupSize          = sc_gpuJgroupSize(pairlistType);
        constexpr int  c_parallelExecutionWidth = sc_gpuParallelExecutionWidth(pairlistType);
        const unsigned superClInteractionMask   = ((1U << c_clusterPerSuperCluster) - 1U);

        // thread/block/warp id-s
        const int tidxi      = threadIdx.x;
        const int tidxj      = threadIdx.y;
        const int tidxz      = threadZ == 1 ? 0 : threadIdx.z;
        const int bidx       = blockIdx.x;
        const int tidx       = tidxi + c_clSize * tidxj;
        const int tidxInWarp = tidx & (c_parallelExecutionWidth - 1);
        const int widx       = (tidxj * c_clSize) / c_subWarp<pairlistType>;
        // Get part for this kernel from global memory. Each block has its own copy to allow asynchronous incrementation.
        int part = plist.d_rollingPruningPart[bidx];
        __syncthreads();

        // Single thread per block increments the block's copy of the part index
        if (tidxi == 0 && tidxj == 0 && tidxz == 0)
        {
            plist.d_rollingPruningPart[bidx] = (part + 1) % numParts;
        }

        // Kernel has been launched with max number of blocks across all passes (plist.nsci/numParts),
        // but the last pass will require 1 less block, so extra block should return early.
        int numSciInPart = divideRoundUp((plist.numSci - part), numParts);
        if (bidx >= numSciInPart)
        {
            return;
        }

        const float3*                    shiftVec         = asFloat3(atdat.shiftVec);
        nbnxn_cj_packed_t<pairlistType>* gm_plistCJPacked = plist.cjPacked;
        const nbnxn_sci_t* plistSci = haveFreshList ? plist.sci : plist.sorting.sciSorted;
        int*               gm_plistSciHistogram = plist.sorting.sciHistogram;
        int*               gm_sciCount          = plist.sorting.sciCount;
        unsigned int*      gm_plistIMask        = plist.imask;
        const float        rlistOuterSq         = nbparam.rlistOuter_sq;
        const float        rlistInnerSq         = nbparam.rlistInner_sq;
        /* shmem buffer for i x+q pre-loading */
        extern __shared__ char sm_dynamicSharedMemory[];
        char*                  sm_nextSlotPtr = sm_dynamicSharedMemory;
        float4*                sm_xqBufferPtr = reinterpret_cast<float4*>(sm_nextSlotPtr);
        sm_nextSlotPtr += incrementSharedMemorySlotPtr<pairlistType, float4>();

        AmdFastBuffer<const float4> xq{ atdat.xq };

        // my i super-cluster's index = sciOffset + current bidx * numParts + part
        const nbnxn_sci_t nbSci          = plistSci[bidx * numParts + part];
        const int         sci            = nbSci.sci;           /* super-cluster */
        const int         cijPackedBegin = nbSci.cjPackedBegin; /* first ...*/
        const int         cijPackedEnd   = nbSci.cjPackedEnd;   /* and last index of j clusters */

        // We may need only a subset of threads active for preloading i-atoms
        // depending on the super-cluster and cluster / thread-block size.
        constexpr bool c_loadUsingAllXYThreads = (c_clSize == c_clusterPerSuperCluster);
        if (tidxz == 0 && (c_loadUsingAllXYThreads || tidxj < c_clusterPerSuperCluster))
        {
            /* Pre-load i-atom x and q into shared memory */
            const int ci = sci * c_clusterPerSuperCluster + tidxj;
            const int ai = ci * c_clSize + tidxi;

            /* We don't need q, but using float4 in shmem avoids bank conflicts.
               (but it also wastes L2 bandwidth). */
            const float4 local = xq[ai];
            const float3 shift = shiftVec[nbSci.shift];
            float4       xi;
            xi.x = -(local.x + shift.x);
            xi.y = -(local.y + shift.y);
            xi.z = -(local.z + shift.z);

            sm_xqBufferPtr[tidxj * c_clSize + tidxi] = xi;
        }
        __syncthreads();

        /* loop over the j clusters = seen by any of the atoms in the current super-cluster.
         * The loop stride c_syclPruneKernelJPackedConcurrency ensures that consecutive warps-pairs
         * are assigned consecutive jPacked's entries. */
        int prunedPairCount = 0;
        for (int jPacked = cijPackedBegin + tidxz; jPacked < cijPackedEnd; jPacked += threadZ)
        {
            int       imaskFull, imaskCheck, imaskNew;
            const int imask = gm_plistCJPacked[jPacked].imei[widx].imask;

            if constexpr (haveFreshList)
            {
                /* Read the mask from the list transferred from the CPU */
                imaskFull = (c_clSizeSq<pairlistType> == c_parallelExecutionWidth)
                                    ? __builtin_amdgcn_readfirstlane(imask)
                                    : imask;
                /* We attempt to prune all pairs present in the original list */
                imaskCheck = imaskFull;
                imaskNew   = 0;
            }
            else
            {
                /* Read the mask from the "warp-pruned" by rlistOuter mask array */
                imaskFull = gm_plistIMask[jPacked * sc_gpuClusterPairSplit(pairlistType) + widx];
                /* "Scalarize" imaskFull when possible, compiler always generates vector loads for that
                 * This means that imaskFull is now stored in a vector register, code is simpler this way
                 */
                imaskFull = (c_clSizeSq<pairlistType> == c_parallelExecutionWidth)
                                    ? __builtin_amdgcn_readfirstlane(imaskFull)
                                    : imaskFull;
                /* Read the old rolling pruned mask, use as a base for new */
                imaskNew = imask;
                /* We only need to check pairs with different mask */
                imaskCheck = (imaskNew ^ imaskFull);
            }

            if (imaskCheck)
            {
#pragma unroll
                for (int jm = 0; jm < c_gpuJGroupSize; jm++)
                {
                    if (imaskCheck & (superClInteractionMask << (jm * c_clusterPerSuperCluster)))
                    {
                        int       mask_ji = (1U << (jm * c_clusterPerSuperCluster));
                        const int cj      = gm_plistCJPacked[jPacked].cj[jm];
                        const int aj      = cj * c_clSize + tidxj;

                        /* load j atom data */
                        const float4          tmp = xq[aj];
                        const AmdPackedFloat3 xj(tmp.x, tmp.y, tmp.z);

                        for (int i = 0; i < c_clusterPerSuperCluster; i++)
                        {
                            if (imaskCheck & mask_ji)
                            {
                                // load i-cluster coordinates from shmem
                                const float4 local = sm_xqBufferPtr[i * c_clSize + tidxi];
                                // distance between i and j atoms
                                AmdPackedFloat3 xi(local.x, local.y, local.z);
                                AmdPackedFloat3 rv = xi + xj;
                                const float     r2 = rv.norm2();

                                /* If _none_ of the atoms pairs are in rlistOuter
                                 * range, the bit corresponding to the current
                                 * cluster-pair in imask gets set to 0. */
                                if (haveFreshList
                                    && !(nb_any_internal<pairlistType>(r2 < rlistOuterSq, widx)))
                                {
                                    imaskFull &= ~mask_ji;
                                }
                                /* If any atom pair is within range, set the bit
                                 * corresponding to the current cluster-pair. */
                                if (nb_any_internal<pairlistType>(r2 < rlistInnerSq, widx))
                                {
                                    imaskNew |= mask_ji;
                                }
                            } // (imaskCheck & mask_ji)
                            /* shift the mask bit by 1 */
                            mask_ji += mask_ji;
                        } // (int i = 0; i < c_clusterPerSuperCluster; i++)
                    } // (imaskCheck & (superClInteractionMask << (jm * c_clusterPerSuperCluster)))
                } // for (int jm = 0; jm < c_gpuJGroupSize; jm++)

                if constexpr (haveFreshList)
                {
                    /* copy the list pruned to rlistOuter to a separate buffer */
                    gm_plistIMask[jPacked * gmx::sc_gpuClusterPairSplit(pairlistType) + widx] = imaskFull;
                    prunedPairCount += __popc(imaskNew);
                }
                /* update the imask with only the pairs up to rlistInner */
                gm_plistCJPacked[jPacked].imei[widx].imask = imaskNew;
            } // (imaskCheck)
            __builtin_amdgcn_wave_barrier();
        }
        if constexpr (haveFreshList)
        {
            if (tidxInWarp == 0)
            {
                if constexpr (threadZ > 1)
                {
                    __syncthreads();
                    char* sm_prunedPairCountBuffer = sm_dynamicSharedMemory;
                    int*  sm_prunedPairCount = reinterpret_cast<int*>(sm_prunedPairCountBuffer);

                    sm_prunedPairCount[tidxz] = prunedPairCount;
                    __syncthreads();

                    for (int indexThreadZ = 1; indexThreadZ < threadZ; indexThreadZ++)
                    {
                        prunedPairCount += sm_prunedPairCount[indexThreadZ];
                    }
                    __syncthreads();
                }

                if (tidxi == 0 && tidxj == 0 && tidxz == 0)
                {
                    int index = max(c_sciHistogramSize - prunedPairCount - 1, 0);
                    atomicAdd(gm_plistSciHistogram + index, 1);
                    gm_sciCount[bidx * numParts + part] = index;
                }
            }
        }
    }
    else
    {
        GMX_UNUSED_VALUE(atdat);
        GMX_UNUSED_VALUE(nbparam);
        GMX_UNUSED_VALUE(plist);
        GMX_UNUSED_VALUE(numParts);
        GMX_DEVICE_ASSERT(false);
    }
}

template<bool haveFreshList, PairlistType pairlistType>
static std::string getPruneKernelName()
{
    static constexpr std::string_view baseName = "pruneOnly";

    static constexpr std::string_view freshListName        = "_freshlist";
    static constexpr std::string_view dummyValue           = "";
    static constexpr std::string_view executionWidth64Name = "_wave64";
    static constexpr std::string_view executionWidth32Name = "_wave32";

    constexpr bool isWave64 = sc_gpuParallelExecutionWidth(pairlistType) == 64;

#if !defined(_MSC_VER)
    return std::string(CompileTimeStringJoin_v<baseName,
                                               (haveFreshList ? freshListName : dummyValue),
                                               (isWave64 ? executionWidth64Name : executionWidth32Name)>);
#else
    std::string returnValue;
    returnValue.reserve(1024);
    return returnValue.append(baseName)
            .append(haveFreshList ? freshListName : dummyValue)
            .append(isWave64 ? executionWidth64Name : executionWidth32Name);
#endif
}

//! \brief Leap Frog HIP prune-only kernel launch code.
template<PairlistType pairlistType, bool haveFreshList, bool hasLargeRegisterPool, class... Args>
void launchNbnxmKernelPruneOnly(const DeviceInformation& deviceInfo,
                                const DeviceStream&      deviceStream,
                                const int                numSciInPart,
                                Args*... args)
{


    constexpr int numThreadZ =
            haveFreshList ? c_pruneKernelJPackedConcurrency / 2 : c_pruneKernelJPackedConcurrency;
    constexpr int c_clSize = sc_gpuClusterSize(pairlistType);
    /* Kernel launch config:
     * - The thread block dimensions match the size of i-clusters, j-clusters,
     *   and j-cluster concurrency, in x, y, and z, respectively.
     * - The 1D block-grid contains as many blocks as super-clusters.
     */
    KernelLaunchConfig config;
    config.blockSize[0] = c_clSize;
    config.blockSize[1] = c_clSize;
    config.blockSize[2] = numThreadZ;
    config.gridSize[0]  = numberOfKernelBlocksSanityCheck(numSciInPart, deviceInfo);
    config.sharedMemorySize = requiredSharedMemorySize<true, numThreadZ, VdwType::Count, pairlistType>();


    auto kernel = nbnxmKernelPruneOnly<pairlistType, hasLargeRegisterPool, haveFreshList, numThreadZ>;
    std::string kernelName = getPruneKernelName<haveFreshList, pairlistType>();

    const auto kernelArgs = prepareGpuKernelArguments(kernel, config, args...);

    launchGpuKernel(kernel, config, deviceStream, nullptr, kernelName.c_str(), kernelArgs);
}

//! \brief Select templated kernel and launch it.
template<PairlistType pairlistType, class... Args>
void chooseAndLaunchNbnxmKernelPruneOnly(bool                     haveFreshList,
                                         const DeviceInformation& deviceInfo,
                                         const DeviceStream&      deviceStream,
                                         const int                numSciInPart,
                                         Args*... args)
{
    const bool hasLargeRegisterPool = targetHasLargeRegisterPool(deviceInfo);


    gmx::dispatchTemplatedFunction(
            [&](auto haveFreshList_, auto hasLargeRegisterPool_)
            {
                launchNbnxmKernelPruneOnly<pairlistType, haveFreshList_, hasLargeRegisterPool_>(
                        deviceInfo, deviceStream, numSciInPart, args...);
            },
            haveFreshList,
            hasLargeRegisterPool);
}

void launchNbnxmKernelPruneOnly(NbnxmGpu* nb, const InteractionLocality iloc, const int* numParts, const int numSciInPart)
{
    NBAtomDataGpu*           adat         = nb->atdat;
    NBParamGpu*              nbp          = nb->nbparam;
    const DeviceStream&      deviceStream = *nb->deviceStreams[iloc];
    const DeviceInformation& deviceInfo   = nb->deviceContext_->deviceInfo();

    std::visit(
            [&](auto&& pairlists)
            {
                auto*      plist            = pairlists[iloc].get();
                const bool haveFreshList    = plist->haveFreshList;
                using T                     = std::decay_t<decltype(*plist)>;
                constexpr auto pairlistType = getPairlistTypeFromPairlist<T>();

                if constexpr (isGpuSpecificPairlist(pairlistType))
                {
                    chooseAndLaunchNbnxmKernelPruneOnly<pairlistType>(
                            haveFreshList, deviceInfo, deviceStream, numSciInPart, adat, nbp, plist, numParts);
                }
            },
            nb->plist);
}

} // namespace gmx
