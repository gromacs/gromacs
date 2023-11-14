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
 * \brief Implements backend-agnostic GPU Force Reduction functions
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "gpuforcereduction_impl.h"

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/mdlib/gpuforcereduction_impl_internal.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

GpuForceReduction::Impl::Impl(const DeviceContext& deviceContext,
                              const DeviceStream&  deviceStream,
                              gmx_wallcycle*       wcycle) :
    baseForce_(),
    deviceContext_(deviceContext),
    deviceStream_(deviceStream),
    nbnxmForceToAdd_(),
    rvecForceToAdd_(),
    forcesReadyNvshmemFlags(nullptr),
    forcesReadyNvshmemFlagsCounter(0),
    wcycle_(wcycle)
{
    cellInfo_.d_cell = nullptr;
}

void GpuForceReduction::Impl::reinit(DeviceBuffer<Float3>  baseForcePtr,
                                     const int             numAtoms,
                                     ArrayRef<const int>   cell,
                                     const int             atomStart,
                                     const bool            accumulate,
                                     GpuEventSynchronizer* completionMarker)
{
    GMX_ASSERT(baseForcePtr, "Input base force for reduction has no data");
    baseForce_        = baseForcePtr;
    numAtoms_         = numAtoms;
    atomStart_        = atomStart;
    accumulate_       = accumulate;
    completionMarker_ = completionMarker;
    cellInfo_.cell    = cell.data();

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    reallocateDeviceBuffer(
            &cellInfo_.d_cell, numAtoms_, &cellInfo_.cellSize, &cellInfo_.cellSizeAlloc, deviceContext_);
    copyToDeviceBuffer(&cellInfo_.d_cell,
                       &(cellInfo_.cell[atomStart]),
                       0,
                       numAtoms_,
                       deviceStream_,
                       GpuApiCallBehavior::Async,
                       nullptr);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);

    dependencyList_.clear();
};

void GpuForceReduction::Impl::registerNbnxmForce(DeviceBuffer<RVec> forcePtr)
{
    GMX_ASSERT(forcePtr, "Input force for reduction has no data");
    nbnxmForceToAdd_ = forcePtr;
};

void GpuForceReduction::Impl::registerRvecForce(DeviceBuffer<RVec> forcePtr)
{
    GMX_ASSERT(forcePtr, "Input force for reduction has no data");
    rvecForceToAdd_ = forcePtr;
};

// NOLINTNEXTLINE(readability-convert-member-functions-to-static, readability-non-const-parameter)
void GpuForceReduction::Impl::registerForcesReadyNvshmemFlags(DeviceBuffer<uint64_t> syncObj)
{
    GMX_ASSERT(syncObj, "Input force for reduction has no data");
    forcesReadyNvshmemFlags = syncObj;
}

void GpuForceReduction::Impl::addDependency(GpuEventSynchronizer* dependency)
{
    GMX_ASSERT(dependency != nullptr, "Force reduction dependency synchronizer should not be NULL");
    dependencyList_.push_back(dependency);
}

void GpuForceReduction::Impl::execute()
{
    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuNBFBufOps);

    if (numAtoms_ != 0)
    {
        GMX_ASSERT(nbnxmForceToAdd_, "Nbnxm force for reduction has no data");

        // Enqueue wait on all dependencies passed
        for (auto* synchronizer : dependencyList_)
        {
            synchronizer->enqueueWaitEvent(deviceStream_);
        }

        const bool addRvecForce = static_cast<bool>(rvecForceToAdd_); // True iff initialized

        if (addRvecForce && forcesReadyNvshmemFlags)
        {
            forcesReadyNvshmemFlagsCounter++;
        }

        launchForceReductionKernel(numAtoms_,
                                   atomStart_,
                                   addRvecForce,
                                   accumulate_,
                                   nbnxmForceToAdd_,
                                   rvecForceToAdd_,
                                   baseForce_,
                                   cellInfo_.d_cell,
                                   deviceStream_,
                                   forcesReadyNvshmemFlags,
                                   forcesReadyNvshmemFlagsCounter);
    }
    else
    {
        /* In case we have nothing to do, but still have dependencies, we need
         * to consume them and mark our own event.
         * Happens sometimes in MdrunVsitesTest.
         * Issue #3988, #4227. */
        for (auto* synchronizer : dependencyList_)
        {
            synchronizer->consume();
        }
    }

    /* Mark that kernel has been launched.
     * Even if we have no work to do and have not launched the kernel, we still mark the event
     * in order to ensure proper marking/consumption balance, see Issue #3988, #4227. */
    if (completionMarker_ != nullptr)
    {
        completionMarker_->markEvent(deviceStream_);
    }

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuNBFBufOps);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

GpuForceReduction::Impl::~Impl()
{
    freeDeviceBuffer(&cellInfo_.d_cell);
}

GpuForceReduction::GpuForceReduction(const DeviceContext& deviceContext,
                                     const DeviceStream&  deviceStream,
                                     gmx_wallcycle*       wcycle) :
    impl_(new Impl(deviceContext, deviceStream, wcycle))
{
}

void GpuForceReduction::registerNbnxmForce(DeviceBuffer<RVec> forcePtr)
{
    impl_->registerNbnxmForce(forcePtr);
}

void GpuForceReduction::registerRvecForce(DeviceBuffer<RVec> forcePtr)
{
    impl_->registerRvecForce(forcePtr);
}

void GpuForceReduction::addDependency(GpuEventSynchronizer* dependency)
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
    impl_->reinit(baseForcePtr, numAtoms, cell, atomStart, accumulate, completionMarker);
}

void GpuForceReduction::execute()
{
    impl_->execute();
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static, readability-non-const-parameter)
void GpuForceReduction::registerForcesReadyNvshmemFlags(DeviceBuffer<uint64_t> syncObj)
{
    GMX_ASSERT(syncObj, "force sync object is NULL");
    impl_->registerForcesReadyNvshmemFlags(syncObj);
}

GpuForceReduction::~GpuForceReduction() = default;

} // namespace gmx
