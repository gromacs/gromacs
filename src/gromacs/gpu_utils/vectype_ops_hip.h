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

#ifndef GMX_GPU_UTILS_VECTYPE_OPS_HIP_H
#define GMX_GPU_UTILS_VECTYPE_OPS_HIP_H

#if !defined(__HIPCC__)
#    error Including header specific for HIP device code without compiling the file with the correct compiler
#endif

#include "gromacs/gpu_utils/vectype_ops_hip_base_math.h"

#include "hip_sycl_kernel_utils.h"
#include "vectype_ops_cuda_hip_shared_declarations.h"
#include "vectype_ops_cuda_hip_shared_trig_math.h"

namespace internal
{

/*!\brief Shuffle based reduction to avoid atomic operations during force accumulation.
 *
 * Uses warp based operations to find out which lanes are currently active and need to contribute
 * to the final reduction.
 */
__device__ __forceinline__ float3 hipHeadSegmentedSum(float3 input, const bool flag, const int localId)
{
    // get bitmask for which lane in the warp is active
    unsigned long long warpFlags = __ballot(flag);

    // shift bitmask to the right once and assign value back
    warpFlags >>= 1;
    const int threadLaneId = (localId & (warpSize - 1));


    // Create bitmask of which elements to accumulate based on which lane
    // we are in the current warp
    warpFlags &= static_cast<unsigned long long>(-1)
                 ^ ((static_cast<unsigned long long>(1) << threadLaneId) - 1U);
    warpFlags |= static_cast<unsigned long long>(1) << (warpSize - 1U);

    // get number of valid entries in our warp by seeking the position
    // of the least significant bit set to 1
    const int validWarpItems = __ffsll(warpFlags);

    float3 output = input;
#pragma unroll
    for (int offset = 1; offset < warpSize; offset *= 2)
    {
        __builtin_assume(offset > 0);
        float3 value = make_float3(__shfl_down(output.x, offset),
                                   __shfl_down(output.y, offset),
                                   __shfl_down(output.z, offset));
        if (threadLaneId + offset < validWarpItems)
        {
            output += value;
        }
    }
    return output;
}

} // namespace internal

__device__ __forceinline__ void staggeredAtomicAddForce(float3* __restrict__ targetPtr,
                                                        float3    input,
                                                        const int index,
                                                        const int localId)
{
    // Combine forces of consecutive lanes that write forces for the same atom
    // but do it only if all lanes are active (the last block may have fewer lanes active or some
    // lanes may not pass tolerance conditions etc.)
    if (static_cast<int>(__popcll(__ballot(1))) == warpSize)
    {
        const int previousLaneId = __shfl_up(index, 1);
        // Is this thread the first one in the lane
        const bool head = (localId & (warpSize - 1)) == 0 || index != previousLaneId;

        input = internal::hipHeadSegmentedSum(input, head, localId);
        if (head)
        {
            // Reduce the number of conflicts that are left after combining consecutive forces
            // by a factor of 3: different lanes write x, y and z in a different order
            int3 offset = { 0, 1, 2 };
            input       = (localId % 3 == 0) ? make_float3(input.y, input.z, input.x) : input;
            offset      = (localId % 3 == 0) ? make_int3(offset.y, offset.z, offset.x) : offset;
            input       = (localId % 3 <= 1) ? make_float3(input.y, input.z, input.x) : input;
            offset      = (localId % 3 <= 1) ? make_int3(offset.y, offset.z, offset.x) : offset;

            auto* baseAddress = &indexedAddress(targetPtr, index)->x;

            atomicAdd(indexedAddress(baseAddress, offset.x), input.x);
            atomicAdd(indexedAddress(baseAddress, offset.y), input.y);
            atomicAdd(indexedAddress(baseAddress, offset.z), input.z);
        }
    }
    else
    {
        atomicAdd(indexedAddress(targetPtr, index), input);
    }
}

#define atomicFetchAdd atomicAdd

__device__ __forceinline__ void atomicFetchAddLocal(float3* targetPtr, const int index, float3 input)
{
    atomicAdd(indexedAddress(targetPtr, index), input);
}

#endif /* VECTYPE_OPS_HPP */
