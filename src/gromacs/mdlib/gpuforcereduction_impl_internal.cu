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
 *
 * \brief Implements GPU Force Reduction using CUDA
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "gpuforcereduction_impl_internal.h"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/gpu_utils/vectype_ops.cuh"

#if GMX_NVSHMEM
#    include <nvshmemx.h>
#endif

namespace gmx
{

constexpr static int c_threadsPerBlock = 128;

template<bool addRvecForce, bool accumulateForce>
static __global__ void reduceKernel(const float3* __restrict__ gm_nbnxmForce,
                                    const float3* __restrict__ rvecForceToAdd,
                                    float3*    gm_fTotal,
                                    const int* gm_cell,
                                    const int  numAtoms,
                                    uint64_t* forcesReadyNvshmemFlags, // NOLINT(readability-non-const-parameter)
                                    const uint64_t forcesReadyNvshmemFlagsCounter)
{
    // map particle-level parallelism to 1D CUDA thread and block index
    const int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;

    float3* gm_fDest = &gm_fTotal[threadIndex];
    float3  temp;

    if (addRvecForce && forcesReadyNvshmemFlags)
    {
#if GMX_NVSHMEM
        if (threadIdx.x == 0)
        {
            nvshmem_signal_wait_until(
                    forcesReadyNvshmemFlags, NVSHMEM_CMP_EQ, forcesReadyNvshmemFlagsCounter);
        }
#else
        GMX_UNUSED_VALUE(forcesReadyNvshmemFlagsCounter);
#endif
    }

    // perform addition for each particle
    if (threadIndex < numAtoms)
    {
        // Accumulate or set nbnxm force
        if (accumulateForce)
        {
            temp = *gm_fDest;
            temp += gm_nbnxmForce[gm_cell[threadIndex]];
        }
        else
        {
            temp = gm_nbnxmForce[gm_cell[threadIndex]];
        }
    }

    if (addRvecForce && forcesReadyNvshmemFlags)
    {
#if GMX_NVSHMEM
        __syncthreads();
#endif
    }

    if (threadIndex < numAtoms)
    {
        if (addRvecForce)
        {
            temp += rvecForceToAdd[threadIndex];
        }

        *gm_fDest = temp;
    }
}

void launchForceReductionKernel(int                        numAtoms,
                                int                        atomStart,
                                bool                       addRvecForce,
                                bool                       accumulate,
                                const DeviceBuffer<Float3> d_nbnxmForceToAdd,
                                const DeviceBuffer<Float3> d_rvecForceToAdd,
                                DeviceBuffer<Float3>       d_baseForce,
                                DeviceBuffer<int>          d_cell,
                                const DeviceStream&        deviceStream,
                                DeviceBuffer<uint64_t>     d_forcesReadyNvshmemFlags,
                                const uint64_t             forcesReadyNvshmemFlagsCounter)
{
    float3* d_baseForcePtr      = &(asFloat3(d_baseForce)[atomStart]);
    float3* d_nbnxmForcePtr     = asFloat3(d_nbnxmForceToAdd);
    float3* d_rvecForceToAddPtr = &(asFloat3(d_rvecForceToAdd)[atomStart]);

    // Configure and launch kernel
    KernelLaunchConfig config;
    config.blockSize[0]     = c_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = divideRoundUp((numAtoms + 1), c_threadsPerBlock);
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    auto kernelFn = addRvecForce
                            ? (accumulate ? reduceKernel<true, true> : reduceKernel<true, false>)
                            : (accumulate ? reduceKernel<false, true> : reduceKernel<false, false>);

    const auto kernelArgs = prepareGpuKernelArguments(kernelFn,
                                                      config,
                                                      &d_nbnxmForcePtr,
                                                      &d_rvecForceToAddPtr,
                                                      &d_baseForcePtr,
                                                      &d_cell,
                                                      &numAtoms,
                                                      &d_forcesReadyNvshmemFlags,
                                                      &forcesReadyNvshmemFlagsCounter);

    launchGpuKernel(kernelFn, config, deviceStream, nullptr, "Force Reduction", kernelArgs);
}

} // namespace gmx
