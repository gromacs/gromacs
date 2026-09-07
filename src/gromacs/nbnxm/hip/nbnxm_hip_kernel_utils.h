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
 *  Utility constant and function declaration for the HIP non-bonded kernels.
 *  This header should be included once at the top level, just before the
 *  kernels are included (has to be preceded by nbnxm_hip_types.h).
 *
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *  \ingroup module_nbnxm
 */

#ifndef NBNXM_HIP_KERNEL_UTILS_H
#define NBNXM_HIP_KERNEL_UTILS_H

#include <assert.h>

#include <type_traits>

/* Note that floating-point constants in HIP code should be suffixed
 * with f (e.g. 0.5f), to stop the compiler producing intermediate
 * code that is in double precision.
 */

#include "gromacs/gpu_utils/gpu_kernel_utils.h"
#include "gromacs/gpu_utils/hip_kernel_utils.h"
#include "gromacs/gpu_utils/hip_sycl_kernel_utils.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/gpu_utils/vectype_ops_hip.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/math/functions.h"
#include "gromacs/nbnxm/nbnxm_enums.h"
#include "gromacs/pbcutil/ishift.h"

#include "nbnxm_hip_types.h"

namespace gmx
{

template<PairlistType pairlistType>
__device__ constexpr int c_subWarp = sc_gpuParallelExecutionWidth(pairlistType);

/*! \brief Log of the i and j cluster size.
 *  change this together with c_clSize !*/
template<PairlistType pairlistType>
__device__ constexpr int c_clSizeLog2 = StaticLog2<sc_gpuClusterSize(pairlistType)>::value;

/*! \brief Square of cluster size. */
template<PairlistType pairlistType>
__device__ constexpr int c_clSizeSq = sc_gpuClusterSize(pairlistType) * sc_gpuClusterSize(pairlistType);

template<PairlistType pairlistType>
__device__ static constexpr int c_fbufStride = c_clSizeSq<pairlistType>;

template<PairlistType pairlistType>
__device__ __forceinline__ int nb_any_internal(int predicate, int widx)
{
    if constexpr (c_subWarp<pairlistType> == deviceWavefrontSize())
    {
        return __any(predicate);
    }
    else
    {
        return static_cast<int>(__ballot(predicate) >> (widx * c_subWarp<pairlistType>));
    }
}

/*! \brief Increment the pointer into shared memory.
 *
 * \tparam T which type we use to calculate the new offset
 */
template<PairlistType pairlistType, typename T>
__device__ inline size_t incrementSharedMemorySlotPtr()
{
    constexpr int offset = sc_gpuClusterPerSuperCluster(pairlistType) * sc_gpuClusterSize(pairlistType);
    return offset * sizeof(T);
}

inline size_t numberOfKernelBlocksSanityCheck(int numSci, const DeviceInformation& deviceInfo)
{
    GMX_ASSERT(numSci > 0, "Grid dimensions of zero are not accepted in HIP");
    const int maximumGridSize = deviceInfo.prop.maxGridSize[0];
    if (numSci > maximumGridSize)
    {
        auto message = formatString(
                "The number of nonbonded work units (number of super clusters) %d exceeds the "
                "maximum grid size in the x dimension %d!",
                numSci,
                maximumGridSize);
        GMX_THROW(InternalError(message));
    }
    return numSci;
}

template<bool isPruneKernel, int numThreadZ, VdwType vdwType, PairlistType pairlistType>
constexpr size_t requiredSharedMemorySize()
{
    constexpr int offset =
            numThreadZ * sc_gpuClusterPerSuperCluster(pairlistType) * sc_gpuClusterSize(pairlistType);
    size_t shmemSize = offset * sizeof(float4);

    if constexpr (!isPruneKernel)
    {
        if constexpr (vdwType == VdwType::CutCombGeom || vdwType == VdwType::CutCombLB)
        {
            shmemSize += offset * sizeof(float2);
        }
        else
        {
            shmemSize += offset * sizeof(int);
        }
        // Need additional storage for pruning data
        shmemSize += 1 * sizeof(int);
    }
    else
    {
        constexpr int pruneKernelOffset = numThreadZ * gmx::sc_gpuClusterPairSplit(pairlistType)
                                          * gmx::sc_gpuJgroupSize(pairlistType);
        shmemSize += pruneKernelOffset * sizeof(int);
    }
    return shmemSize;
}

//! Find out if the target device has a large enough register pool (MI2xx and later)
inline bool targetHasLargeRegisterPool(const DeviceInformation& deviceInfo)
{
    return deviceInfo.deviceHasLargeRegisterPool;
}

__device__ static inline float2 fetchNbfpC6C12(const float2* nbfpComb, int type)
{
    return *indexedAddress(nbfpComb, type);
}

/*! \brief Reduce c_clSize j-force components using AMD DPP instruction.
 *
 * c_clSize consecutive threads hold the force components of a j-atom which we
 * reduced in log2(cl_Size) steps using shift
 *
 * Note: This causes massive amount of spills with the tabulated kernel on gfx803 using ROCm 5.3.
 * We don't disable it only for the tabulated kernel as the analytical is the default anyway.
 */
template<PairlistType pairlistType>
__device__ static inline float reduceForceJWarpShuffle(AmdPackedFloat3 f, const int tidxi)
{
    static_assert(isPowerOfTwo(sc_gpuClusterPerSuperCluster(pairlistType)));
    static_assert(sc_gpuClusterSize(pairlistType) == 8);
    f[0] += amdDppUpdateShfl<float, /* row_shl:1 */ 0x101>(f[0]);
    f[1] += amdDppUpdateShfl<float, /* row_shr:1 */ 0x111>(f[1]);
    f[2] += amdDppUpdateShfl<float, /* row_shl:1 */ 0x101>(f[2]);

    if (tidxi & 1)
    {
        f[0] = f[1];
    }

    f[0] += amdDppUpdateShfl<float, /* row_shl:2 */ 0x102>(f[0]);
    f[2] += amdDppUpdateShfl<float, /* row_shr:2 */ 0x112>(f[2]);

    if (tidxi & 2)
    {
        f[0] = f[2];
    }

    f[0] += amdDppUpdateShfl<float, /* row_shl:4 */ 0x104>(f[0]);
    return f[0];
}

/*! \brief Lowest level i force reduction.
 *
 * Only works for array sizes that are power of 2.
 * Uses AMD DPP instructions to avoid use of atomic operations.
 */
template<PairlistType pairlistType>
__device__ static inline float reduceForceIWarpShuffle(AmdPackedFloat3 f, const int tidxi, const int tidxj)
{
    static_assert(isPowerOfTwo(sc_gpuClusterPerSuperCluster(pairlistType)));
    static_assert(sc_gpuClusterSize(pairlistType) == 8);
    constexpr int c_clSize = sc_gpuClusterSize(pairlistType);

    // transpose first to enable later DDP based reduction
    f[0] = __shfl(f[0], tidxi * c_clSize + tidxj);
    f[1] = __shfl(f[1], tidxi * c_clSize + tidxj);
    f[2] = __shfl(f[2], tidxi * c_clSize + tidxj);

    f[0] += amdDppUpdateShfl<float, /* row_shl:1 */ 0x101>(f[0]);
    f[1] += amdDppUpdateShfl<float, /* row_shr:1 */ 0x111>(f[1]);
    f[2] += amdDppUpdateShfl<float, /* row_shl:1 */ 0x101>(f[2]);

    if (tidxi & 1)
    {
        f[0] = f[1];
    }

    f[0] += amdDppUpdateShfl<float, /* row_shl:2 */ 0x102>(f[0]);
    f[2] += amdDppUpdateShfl<float, /* row_shr:2 */ 0x112>(f[2]);

    if (tidxi & 2)
    {
        f[0] = f[2];
    }

    f[0] += amdDppUpdateShfl<float, /* row_shl:4 */ 0x104>(f[0]);

    return f[0];
}

/*! \brief Lowest level i force reduction.
 *
 * Only works for array sizes that are power of 2.
 * Uses atomic operations instead of shuffles.
 */
template<PairlistType pairlistType, bool calculateShift>
__device__ static inline float3
reduceForceIAtomics(AmdPackedFloat3 input, float3* result, const int tidxj, const int aidx)
{

    static_assert(isPowerOfTwo(sc_gpuClusterPerSuperCluster(pairlistType)));
    constexpr int c_clSize                 = sc_gpuClusterSize(pairlistType);
    constexpr int c_parallelExecutionWidth = sc_gpuParallelExecutionWidth(pairlistType);

#pragma unroll
    for (int offset = c_parallelExecutionWidth >> 1; offset >= c_clSize; offset >>= 1)
    {
        input[0] += __shfl_down(input[0], offset);
        input[1] += __shfl_down(input[1], offset);
        input[2] += __shfl_down(input[2], offset);
    }

    float3 shiftForce = make_float3(0.0F);
    if (tidxj % (c_parallelExecutionWidth / c_clSize) == 0)
    {
        atomicAdd(&(result[aidx].x), input[0]);
        atomicAdd(&(result[aidx].y), input[1]);
        atomicAdd(&(result[aidx].z), input[2]);
        if constexpr (calculateShift)
        {
            shiftForce.x = input[0];
            shiftForce.y = input[1];
            shiftForce.z = input[2];
        }
    }
    return shiftForce;
}

/*! \brief Reduce i forces.
 *
 * Only works for array sizes that are power of 2.
 * Depending on architecture, reduce using DPP shuffles for main forces
 * or atomics.  Final accumulation is always done using atomics, while
 * shift forces are using DPP shuffles for non CDNA architectures.
 */
template<bool calculateShift, PairlistType pairlistType>
__device__ static inline void reduceForceI(AmdPackedFloat3* input,
                                           float3*          result,
                                           const int        tidxi,
                                           const int        tidxj,
                                           const int        tidx,
                                           const int        sci,
                                           float3*          fShift,
                                           const int        shiftBase)
{
    constexpr int c_clusterPerSuperCluster = sc_gpuClusterPerSuperCluster(pairlistType);
    constexpr int c_clSize                 = sc_gpuClusterSize(pairlistType);

    if constexpr (gmx::sc_gpuParallelExecutionWidth(pairlistType) == 64)
    {
        float shiftForceBuffer = 0.0F;
        float fci[c_clusterPerSuperCluster];
        for (int i = 0; i < c_clusterPerSuperCluster; i++)
        {
            fci[i] = reduceForceIWarpShuffle<pairlistType>(input[i], tidxi, tidxj);
            shiftForceBuffer += fci[i];
        }
        if (tidxi < 3)
        {
            for (int i = 0; i < c_clusterPerSuperCluster; i++)
            {
                const int ai = (sci * c_clusterPerSuperCluster + i) * c_clSize + tidxj;
                amdFastAtomicAddForce(result, ai, tidxi, fci[i]);
            }
        }

        if constexpr (calculateShift)
        {
            if (tidxi < 3)
            {
                amdFastAtomicAddForce(fShift, shiftBase, tidxi, shiftForceBuffer);
            }
        }
    }
    else
    {
        float3 shiftForceBuffer = make_float3(0.0F);
        for (int i = 0; i < c_clusterPerSuperCluster; i++)
        {
            const int ai = (sci * c_clusterPerSuperCluster + i) * c_clSize + tidxi;
            shiftForceBuffer +=
                    reduceForceIAtomics<pairlistType, calculateShift>(input[i], result, tidxj, ai);
        }

        if constexpr (calculateShift)
        {
            shiftForceBuffer.x += amdDppUpdateShfl<float, 0xb1>(shiftForceBuffer.x);
            shiftForceBuffer.y += amdDppUpdateShfl<float, 0xb1>(shiftForceBuffer.y);
            shiftForceBuffer.z += amdDppUpdateShfl<float, 0xb1>(shiftForceBuffer.z);

            shiftForceBuffer.x += amdDppUpdateShfl<float, 0x4e>(shiftForceBuffer.x);
            shiftForceBuffer.y += amdDppUpdateShfl<float, 0x4e>(shiftForceBuffer.y);
            shiftForceBuffer.z += amdDppUpdateShfl<float, 0x4e>(shiftForceBuffer.z);

            shiftForceBuffer.x += amdDppUpdateShfl<float, 0x114>(shiftForceBuffer.x);
            shiftForceBuffer.y += amdDppUpdateShfl<float, 0x114>(shiftForceBuffer.y);
            shiftForceBuffer.z += amdDppUpdateShfl<float, 0x114>(shiftForceBuffer.z);
            if (tidx == (c_clSize - 1) || tidx == (c_subWarp<pairlistType> + c_clSize - 1))
            {

                atomicAdd(&(fShift[shiftBase].x), shiftForceBuffer.x);
                atomicAdd(&(fShift[shiftBase].y), shiftForceBuffer.y);
                atomicAdd(&(fShift[shiftBase].z), shiftForceBuffer.z);
            }
        }
    }
}

/*! \brief Energy reduction kernel for regular nbnxm kernels.
 *
 * Only works for power of two array sizes.
 */
template<PairlistType pairlistType>
__device__ static inline void
reduceEnergyWarpShuffle(float localLJ, float localEl, float* gm_LJ, float* gm_El, int tidx)
{
    static_assert(isPowerOfTwo(sc_gpuClusterPerSuperCluster(pairlistType)));
    constexpr int c_parallelExecutionWidth = sc_gpuParallelExecutionWidth(pairlistType);
    localLJ += amdDppUpdateShfl<float, 0xb1>(localLJ);
    localEl += amdDppUpdateShfl<float, 0xb1>(localEl);

    localLJ += amdDppUpdateShfl<float, 0x4e>(localLJ);
    localEl += amdDppUpdateShfl<float, 0x4e>(localEl);

    // DPP: row_shr:4
    localLJ += amdDppUpdateShfl<float, 0x114>(localLJ);
    localEl += amdDppUpdateShfl<float, 0x114>(localEl);

    // DPP: row_shr:8
    localLJ += amdDppUpdateShfl<float, 0x118>(localLJ);
    localEl += amdDppUpdateShfl<float, 0x118>(localEl);

    // only CDNA (wave64) devices support broadcasts
    // so we can only use those dpp instructions on devices
    // with wave64 (aka CDNA)
    if constexpr (deviceWavefrontSize() == 64)
    {
        // DPP: row_bcast15 (Broadcast thread 15 of each row to next row)
        localLJ += amdDppUpdateShfl<float, 0x142>(localLJ);
        localEl += amdDppUpdateShfl<float, 0x142>(localEl);

        if constexpr (sc_gpuParallelExecutionWidth(pairlistType) == 64)
        {

            // DPP: row_bcast31 (Broadcast thread 31 to rows 2 and 3)
            localLJ += amdDppUpdateShfl<float, 0x143>(localLJ);
            localEl += amdDppUpdateShfl<float, 0x143>(localEl);
        }
    }
    else
    {
        localLJ += __shfl(localLJ, 15);
        localEl += __shfl(localEl, 15);
    }

    /* The last thread in the subWarp writes the reduced energies */
    if ((tidx & (c_parallelExecutionWidth - 1)) == (c_parallelExecutionWidth - 1))
    {
        atomicAdd(gm_LJ, localLJ);
        atomicAdd(gm_El, localEl);
    }
}

} // namespace gmx

#endif /* NBNXM_HIP_KERNEL_UTILS_HPP */
