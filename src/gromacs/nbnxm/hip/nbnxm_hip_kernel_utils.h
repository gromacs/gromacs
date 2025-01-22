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
 *  kernels are included (has to be preceded by nbnxn_hip_types.h).
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

#include "gromacs/gpu_utils/hip_kernel_utils.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/gpu_utils/vectype_ops_hip.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/nbnxm/nbnxm_enums.h"

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
    if constexpr (c_subWarp<pairlistType> == warpSize)
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

} // namespace gmx

#endif /* NBNXN_HIP_KERNEL_UTILS_HPP */
