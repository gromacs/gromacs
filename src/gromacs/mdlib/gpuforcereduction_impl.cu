/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
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

#include "gpuforcereduction_impl.cuh"

#include <stdio.h>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.cuh"
#include "gromacs/gpu_utils/typecasts.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/utility/gmxassert.h"

#include "gpuforcereduction.h"

namespace gmx
{

constexpr static int c_threadsPerBlock = 128;

typedef struct rvecDeviceForceData rvecDeviceForceData_t;


template<bool addRvecForce, bool accumulateForce>
static __global__ void reduceKernel(const float3* __restrict__ gm_nbnxmForce,
                                    const float3* __restrict__ rvecForceToAdd,
                                    float3*    gm_fTotal,
                                    const int* gm_cell,
                                    const int  numAtoms)
{

    // map particle-level parallelism to 1D CUDA thread and block index
    const int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;

    // perform addition for each particle
    if (threadIndex < numAtoms)
    {

        float3* gm_fDest = &gm_fTotal[threadIndex];
        float3  temp;

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

        if (addRvecForce)
        {
            temp += rvecForceToAdd[threadIndex];
        }

        *gm_fDest = temp;
    }
    return;
}

GpuForceReduction::Impl::Impl(const DeviceContext& deviceContext,
                              const DeviceStream&  deviceStream,
                              gmx_wallcycle*       wcycle) :
    deviceContext_(deviceContext),
    deviceStream_(deviceStream),
    wcycle_(wcycle){};

void GpuForceReduction::Impl::reinit(float3*               baseForcePtr,
                                     const int             numAtoms,
                                     ArrayRef<const int>   cell,
                                     const int             atomStart,
                                     const bool            accumulate,
                                     GpuEventSynchronizer* completionMarker)
{
    GMX_ASSERT((baseForcePtr != nullptr), "Input base force for reduction has no data");
    baseForce_        = &(baseForcePtr[atomStart]);
    numAtoms_         = numAtoms;
    atomStart_        = atomStart;
    accumulate_       = accumulate;
    completionMarker_ = completionMarker;
    cellInfo_.cell    = cell.data();

    wallcycle_start_nocount(wcycle_, ewcLAUNCH_GPU);
    reallocateDeviceBuffer(&cellInfo_.d_cell, numAtoms_, &cellInfo_.cellSize,
                           &cellInfo_.cellSizeAlloc, deviceContext_);
    copyToDeviceBuffer(&cellInfo_.d_cell, &(cellInfo_.cell[atomStart]), 0, numAtoms_, deviceStream_,
                       GpuApiCallBehavior::Async, nullptr);
    wallcycle_stop(wcycle_, ewcLAUNCH_GPU);

    dependencyList_.clear();
};

void GpuForceReduction::Impl::registerNbnxmForce(DeviceBuffer<RVec> forcePtr)
{
    GMX_ASSERT((forcePtr != nullptr), "Input force for reduction has no data");
    nbnxmForceToAdd_ = forcePtr;
};

void GpuForceReduction::Impl::registerRvecForce(DeviceBuffer<RVec> forcePtr)
{
    GMX_ASSERT((forcePtr != nullptr), "Input force for reduction has no data");
    rvecForceToAdd_ = forcePtr;
};

void GpuForceReduction::Impl::addDependency(GpuEventSynchronizer* const dependency)
{
    dependencyList_.push_back(dependency);
}

void GpuForceReduction::Impl::execute()
{
    wallcycle_start_nocount(wcycle_, ewcLAUNCH_GPU);
    wallcycle_sub_start(wcycle_, ewcsLAUNCH_GPU_NB_F_BUF_OPS);

    if (numAtoms_ == 0)
    {
        return;
    }

    GMX_ASSERT((nbnxmForceToAdd_ != nullptr), "Nbnxm force for reduction has no data");

    // Enqueue wait on all dependencies passed
    for (auto const synchronizer : dependencyList_)
    {
        synchronizer->enqueueWaitEvent(deviceStream_);
    }

    float3* d_nbnxmForce     = asFloat3(nbnxmForceToAdd_);
    float3* d_rvecForceToAdd = &(asFloat3(rvecForceToAdd_)[atomStart_]);

    // Configure and launch kernel
    KernelLaunchConfig config;
    config.blockSize[0]     = c_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = ((numAtoms_ + 1) + c_threadsPerBlock - 1) / c_threadsPerBlock;
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    auto kernelFn = (rvecForceToAdd_ != nullptr)
                            ? (accumulate_ ? reduceKernel<true, true> : reduceKernel<true, false>)
                            : (accumulate_ ? reduceKernel<false, true> : reduceKernel<false, false>);

    const auto kernelArgs = prepareGpuKernelArguments(kernelFn, config, &d_nbnxmForce, &d_rvecForceToAdd,
                                                      &baseForce_, &cellInfo_.d_cell, &numAtoms_);

    launchGpuKernel(kernelFn, config, deviceStream_, nullptr, "Force Reduction", kernelArgs);

    // Mark that kernel has been launched
    if (completionMarker_ != nullptr)
    {
        completionMarker_->markEvent(deviceStream_);
    }

    wallcycle_sub_stop(wcycle_, ewcsLAUNCH_GPU_NB_F_BUF_OPS);
    wallcycle_stop(wcycle_, ewcLAUNCH_GPU);
}

GpuForceReduction::Impl::~Impl(){};

GpuForceReduction::GpuForceReduction(const DeviceContext& deviceContext,
                                     const DeviceStream&  deviceStream,
                                     gmx_wallcycle*       wcycle) :
    impl_(new Impl(deviceContext, deviceStream, wcycle))
{
}

void GpuForceReduction::registerNbnxmForce(void* forcePtr)
{
    impl_->registerNbnxmForce(reinterpret_cast<DeviceBuffer<RVec>>(forcePtr));
}

void GpuForceReduction::registerRvecForce(void* forcePtr)
{
    impl_->registerRvecForce(reinterpret_cast<DeviceBuffer<RVec>>(forcePtr));
}

void GpuForceReduction::addDependency(GpuEventSynchronizer* const dependency)
{
    impl_->addDependency(dependency);
}

void GpuForceReduction::reinit(DeviceBuffer<RVec>    baseForcePtr,
                               const int             numAtoms,
                               ArrayRef<const int>   cell,
                               const int             atomStart,
                               const bool            accumulate,
                               GpuEventSynchronizer* completionMarker)
{
    impl_->reinit(asFloat3(baseForcePtr), numAtoms, cell, atomStart, accumulate, completionMarker);
}
void GpuForceReduction::execute()
{
    impl_->execute();
}

GpuForceReduction::~GpuForceReduction() = default;

} // namespace gmx
