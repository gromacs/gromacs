/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief CUDA kernels for virtual site (vsite) construction and force spreading.
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "vsite_gpu_internal.h"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/utility/functions.h"

#include "vsite_gpu_internal_shared.h"

namespace gmx
{

//! Threads per block for all vsite kernels.
constexpr static int sc_threadsPerBlock = 128;

//! \brief Constructs vsite \p idx of type \p vsiteType.
template<InteractionFunction vsiteType>
static __device__ inline void constructVSiteGpu(int idx,
                                                const int* __restrict__ gm_iatoms,
                                                const float* __restrict__ gm_params,
                                                float3* __restrict__ gm_x,
                                                const PbcAiuc& pbcAiuc)
{
    if constexpr (vsiteType == InteractionFunction::VirtualSite1)
    {
        constructVSite1Gpu(idx, gm_iatoms, gm_x);
    }
    else if constexpr (vsiteType == InteractionFunction::VirtualSite2)
    {
        constructVSite2Gpu(idx, gm_iatoms, gm_params, gm_x, pbcAiuc);
    }
    else if constexpr (vsiteType == InteractionFunction::VirtualSite3)
    {
        constructVSite3Gpu(idx, gm_iatoms, gm_params, gm_x, pbcAiuc);
    }
}

//! \brief Spreads the force on vsite \p idx of type \p vsiteType to its constructing atoms.
template<InteractionFunction vsiteType>
static __device__ inline void spreadVSiteForceGpu(int idx,
                                                  const int* __restrict__ gm_iatoms,
                                                  const float* __restrict__ gm_params,
                                                  float3* __restrict__ gm_f)
{
    if constexpr (vsiteType == InteractionFunction::VirtualSite1)
    {
        spreadVSite1ForceGpu(idx, gm_iatoms, gm_f);
    }
    else if constexpr (vsiteType == InteractionFunction::VirtualSite2)
    {
        spreadVSite2ForceGpu(idx, gm_iatoms, gm_params, gm_f);
    }
    else if constexpr (vsiteType == InteractionFunction::VirtualSite3)
    {
        spreadVSite3ForceGpu(idx, gm_iatoms, gm_params, gm_f);
    }
}

//! \brief VSite construction kernel for all supported vsite types. Each thread constructs one vsite.
template<InteractionFunction vsiteType>
__launch_bounds__(sc_threadsPerBlock) __global__
        void constructVSiteKernel(int numVsites,
                                  const int* __restrict__ gm_iatoms,
                                  const float* __restrict__ gm_params,
                                  float3* __restrict__ gm_x,
                                  const PbcAiuc pbcAiuc)
{
    int idx = blockIdx.x * sc_threadsPerBlock + threadIdx.x;
    if (idx >= numVsites)
    {
        return;
    }
    constructVSiteGpu<vsiteType>(idx, gm_iatoms, gm_params, gm_x, pbcAiuc);
}

//! \brief VSite force spreading kernel for all supported vsite types. Each thread spreads one vsite's force.
template<InteractionFunction vsiteType>
__launch_bounds__(sc_threadsPerBlock) __global__
        void spreadVSiteForceKernel(int numVsites,
                                    const int* __restrict__ gm_iatoms,
                                    const float* __restrict__ gm_params,
                                    float3* __restrict__ gm_f)
{
    int idx = blockIdx.x * sc_threadsPerBlock + threadIdx.x;
    if (idx >= numVsites)
    {
        return;
    }
    spreadVSiteForceGpu<vsiteType>(idx, gm_iatoms, gm_params, gm_f);
}

/*! \brief Fused vsite construction kernel: constructs the vsites of all types in one launch.
 *
 * Each thread constructs one vsite; the type of vsite a thread works on is determined by the
 * thread ranges in \p kernelParams. As all types are handled concurrently, this kernel can only be
 * used when no vsite is a constructing atom of another vsite.
 */
__launch_bounds__(sc_threadsPerBlock) __global__
        void constructVSitesFusedKernel(VsiteGpuKernelParameters kernelParams,
                                        VsiteGpuKernelBuffers    kernelBuffers,
                                        float3* __restrict__ gm_x)
{
    const int tid = blockIdx.x * sc_threadsPerBlock + threadIdx.x;
    // The thread ranges are unrolled, so that the buffers are indexed at compile time.
    // A run-time index would put the kernel argument structs in local memory.
#pragma unroll
    for (int j = 0; j < c_numVSiteTypesOnGpu; j++)
    {
        if (tid <= kernelParams.typeRangeEnd[j])
        {
            const int idx = tid - kernelParams.typeRangeStart[j];
            // Threads in the padding at the end of a range have no vsite to construct
            if (idx < kernelParams.numVsites[j])
            {
                // We need to switch-case by thread range index here instead of the more sightly
                // vsite type, because the latter is a run-time value and would put the kernel
                // argument structs silently into local memory, at a large cost.
                switch (j)
                {
                    case 0:
                        constructVSiteGpu<InteractionFunction::VirtualSite1>(idx,
                                                                             kernelBuffers.d_iAtoms[j],
                                                                             kernelBuffers.d_params[j],
                                                                             gm_x,
                                                                             kernelParams.pbcAiuc);
                        break;
                    case 1:
                        constructVSiteGpu<InteractionFunction::VirtualSite2>(idx,
                                                                             kernelBuffers.d_iAtoms[j],
                                                                             kernelBuffers.d_params[j],
                                                                             gm_x,
                                                                             kernelParams.pbcAiuc);
                        break;
                    case 2:
                        constructVSiteGpu<InteractionFunction::VirtualSite3>(idx,
                                                                             kernelBuffers.d_iAtoms[j],
                                                                             kernelBuffers.d_params[j],
                                                                             gm_x,
                                                                             kernelParams.pbcAiuc);
                        break;
                    default: break;
                }
            }
            break;
        }
    }
}

/*! \brief Fused vsite force spreading kernel: spreads the forces of all types in one launch.
 *
 * Each thread spreads the force of one vsite; the type of vsite a thread works on is determined
 * by the thread ranges in \p kernelParams. As all types are handled concurrently, this kernel can
 * only be used when no vsite is a constructing atom of another vsite.
 */
__launch_bounds__(sc_threadsPerBlock) __global__
        void spreadVSiteForcesFusedKernel(VsiteGpuKernelParameters kernelParams,
                                          VsiteGpuKernelBuffers    kernelBuffers,
                                          float3* __restrict__ gm_f)
{
    const int tid = blockIdx.x * sc_threadsPerBlock + threadIdx.x;
    // The thread ranges are unrolled, so that the buffers are indexed at compile time.
    // A run-time index would put the kernel argument structs in local memory.
#pragma unroll
    for (int j = 0; j < c_numVSiteTypesOnGpu; j++)
    {
        if (tid <= kernelParams.typeRangeEnd[j])
        {
            const int idx = tid - kernelParams.typeRangeStart[j];
            // Threads in the padding at the end of a range have no vsite force to spread
            if (idx < kernelParams.numVsites[j])
            {
                // NOTE: the case labels are thread range indices, so the vsite types have to be
                // listed in the order of c_vsiteTypesOnGpu. Dispatching on c_vsiteTypesOnGpu[j]
                // instead is a run-time value, which makes the compiler place the kernel argument
                // structs in local memory, at a large cost.
                switch (j)
                {
                    case 0:
                        spreadVSiteForceGpu<InteractionFunction::VirtualSite1>(
                                idx, kernelBuffers.d_iAtoms[j], kernelBuffers.d_params[j], gm_f);
                        break;
                    case 1:
                        spreadVSiteForceGpu<InteractionFunction::VirtualSite2>(
                                idx, kernelBuffers.d_iAtoms[j], kernelBuffers.d_params[j], gm_f);
                        break;
                    case 2:
                        spreadVSiteForceGpu<InteractionFunction::VirtualSite3>(
                                idx, kernelBuffers.d_iAtoms[j], kernelBuffers.d_params[j], gm_f);
                        break;
                    default: break;
                }
            }
            break;
        }
    }
}

// ---------------------------------------------------------------------------
// Host kernel launch functions
// ---------------------------------------------------------------------------

template<InteractionFunction vsiteType>
void launchConstructVSiteKernel(const VsiteGpuKernelParameters& kernelParams,
                                const VsiteGpuKernelBuffers&    kernelBuffers,
                                DeviceBuffer<Float3>            d_x,
                                const DeviceStream&             deviceStream)
{
    constexpr int typeIndex = vsiteTypeIndex(vsiteType);
    int           numVsites = kernelParams.numVsites[typeIndex];
    if (numVsites == 0)
    {
        return;
    }
    DeviceBuffer<t_iatom> d_iatoms = kernelBuffers.d_iAtoms[typeIndex];
    DeviceBuffer<float>   d_params = kernelBuffers.d_params[typeIndex];
    KernelLaunchConfig    config;
    config.blockSize[0]     = sc_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = divideRoundUp(numVsites, sc_threadsPerBlock);
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    auto kernelPtr  = constructVSiteKernel<vsiteType>;
    auto kernelArgs = prepareGpuKernelArguments(kernelPtr,
                                                config,
                                                &numVsites,
                                                &d_iatoms,
                                                &d_params,
                                                asFloat3Pointer(&d_x),
                                                &kernelParams.pbcAiuc);
    launchGpuKernel(kernelPtr, config, deviceStream, nullptr, "constructVSite", kernelArgs);
}

template<InteractionFunction vsiteType>
void launchSpreadVSiteForceKernel(const VsiteGpuKernelParameters& kernelParams,
                                  const VsiteGpuKernelBuffers&    kernelBuffers,
                                  DeviceBuffer<Float3>            d_f,
                                  const DeviceStream&             deviceStream)
{
    constexpr int typeIndex = vsiteTypeIndex(vsiteType);
    int           numVsites = kernelParams.numVsites[typeIndex];
    if (numVsites == 0)
    {
        return;
    }
    DeviceBuffer<t_iatom> d_iatoms = kernelBuffers.d_iAtoms[typeIndex];
    DeviceBuffer<float>   d_params = kernelBuffers.d_params[typeIndex];
    KernelLaunchConfig    config;
    config.blockSize[0]     = sc_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = divideRoundUp(numVsites, sc_threadsPerBlock);
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    auto kernelPtr  = spreadVSiteForceKernel<vsiteType>;
    auto kernelArgs = prepareGpuKernelArguments(
            kernelPtr, config, &numVsites, &d_iatoms, &d_params, asFloat3Pointer(&d_f));
    launchGpuKernel(kernelPtr, config, deviceStream, nullptr, "spreadVSiteForce", kernelArgs);
}

void launchConstructVSitesFusedKernel(const VsiteGpuKernelParameters& kernelParams,
                                      const VsiteGpuKernelBuffers&    kernelBuffers,
                                      DeviceBuffer<Float3>            d_x,
                                      const DeviceStream&             deviceStream)
{
    const int numThreads = kernelParams.typeRangeEnd[c_numVSiteTypesOnGpu - 1] + 1;
    if (numThreads == 0)
    {
        return;
    }
    KernelLaunchConfig config;
    config.blockSize[0]     = sc_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = divideRoundUp(numThreads, sc_threadsPerBlock);
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    auto kernelPtr  = constructVSitesFusedKernel;
    auto kernelArgs = prepareGpuKernelArguments(
            kernelPtr, config, &kernelParams, &kernelBuffers, asFloat3Pointer(&d_x));
    launchGpuKernel(kernelPtr, config, deviceStream, nullptr, "constructVSitesFused", kernelArgs);
}

void launchSpreadVSiteForcesFusedKernel(const VsiteGpuKernelParameters& kernelParams,
                                        const VsiteGpuKernelBuffers&    kernelBuffers,
                                        DeviceBuffer<Float3>            d_f,
                                        const DeviceStream&             deviceStream)
{
    const int numThreads = kernelParams.typeRangeEnd[c_numVSiteTypesOnGpu - 1] + 1;
    if (numThreads == 0)
    {
        return;
    }
    KernelLaunchConfig config;
    config.blockSize[0]     = sc_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = divideRoundUp(numThreads, sc_threadsPerBlock);
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    auto kernelPtr  = spreadVSiteForcesFusedKernel;
    auto kernelArgs = prepareGpuKernelArguments(
            kernelPtr, config, &kernelParams, &kernelBuffers, asFloat3Pointer(&d_f));
    launchGpuKernel(kernelPtr, config, deviceStream, nullptr, "spreadVSiteForcesFused", kernelArgs);
}

//! \brief Explicitly instantiate the launch functions for every vsite type supported on GPU.
#define INSTANTIATE_VSITE_LAUNCHERS(vsiteType)                                             \
    template void launchConstructVSiteKernel<vsiteType>(const VsiteGpuKernelParameters&,   \
                                                        const VsiteGpuKernelBuffers&,      \
                                                        DeviceBuffer<Float3>,              \
                                                        const DeviceStream&);              \
    template void launchSpreadVSiteForceKernel<vsiteType>(const VsiteGpuKernelParameters&, \
                                                          const VsiteGpuKernelBuffers&,    \
                                                          DeviceBuffer<Float3>,            \
                                                          const DeviceStream&);

INSTANTIATE_VSITE_LAUNCHERS(InteractionFunction::VirtualSite1)
INSTANTIATE_VSITE_LAUNCHERS(InteractionFunction::VirtualSite2)
INSTANTIATE_VSITE_LAUNCHERS(InteractionFunction::VirtualSite3)

#undef INSTANTIATE_VSITE_LAUNCHERS

} // namespace gmx
