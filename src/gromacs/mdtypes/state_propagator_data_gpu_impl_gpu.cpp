/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief Definitions of interfaces for GPU state data propagator object.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include "config.h"

#if GMX_GPU

#    include "gromacs/gpu_utils/device_stream_manager.h"
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/gpueventsynchronizer.h"
#    include "gromacs/math/vectypes.h"
#    include "gromacs/mdtypes/state_propagator_data_gpu.h"
#    include "gromacs/timing/wallcycle.h"
#    include "gromacs/utility/classhelpers.h"

#    include "state_propagator_data_gpu_impl.h"


namespace gmx
{

StatePropagatorDataGpu::Impl::Impl(const DeviceStreamManager& deviceStreamManager,
                                   GpuApiCallBehavior         transferKind,
                                   int                        allocationBlockSizeDivisor,
                                   gmx_wallcycle*             wcycle) :
    deviceContext_(deviceStreamManager.context()),
    transferKind_(transferKind),
    allocationBlockSizeDivisor_(allocationBlockSizeDivisor),
    wcycle_(wcycle)
{
    static_assert(
            GMX_GPU,
            "GPU state propagator data object should only be constructed on the GPU code-paths.");

    // We need to keep local copies for re-initialization.
    pmeStream_      = deviceStreamManager.streamIsValid(DeviceStreamType::Pme)
                              ? &deviceStreamManager.stream(DeviceStreamType::Pme)
                              : nullptr;
    localStream_    = deviceStreamManager.streamIsValid(DeviceStreamType::NonBondedLocal)
                              ? &deviceStreamManager.stream(DeviceStreamType::NonBondedLocal)
                              : nullptr;
    nonLocalStream_ = deviceStreamManager.streamIsValid(DeviceStreamType::NonBondedNonLocal)
                              ? &deviceStreamManager.stream(DeviceStreamType::NonBondedNonLocal)
                              : nullptr;
    updateStream_   = deviceStreamManager.streamIsValid(DeviceStreamType::UpdateAndConstraints)
                              ? &deviceStreamManager.stream(DeviceStreamType::UpdateAndConstraints)
                              : nullptr;


    // Map the atom locality to the stream that will be used for coordinates,
    // velocities and forces transfers. Same streams are used for H2D and D2H copies.
    // Note, that nullptr stream is used here to indicate that the copy is not supported.
    xCopyStreams_[AtomLocality::Local]    = updateStream_;
    xCopyStreams_[AtomLocality::NonLocal] = nonLocalStream_;
    xCopyStreams_[AtomLocality::All]      = nullptr;

    vCopyStreams_[AtomLocality::Local]    = updateStream_;
    vCopyStreams_[AtomLocality::NonLocal] = nullptr;
    vCopyStreams_[AtomLocality::All]      = nullptr;

    fCopyStreams_[AtomLocality::Local]    = localStream_;
    fCopyStreams_[AtomLocality::NonLocal] = nonLocalStream_;
    fCopyStreams_[AtomLocality::All]      = updateStream_;

    copyInStream_ = std::make_unique<DeviceStream>(deviceContext_, DeviceStreamPriority::Normal, false);
    memsetStream_ = std::make_unique<DeviceStream>(deviceContext_, DeviceStreamPriority::Normal, false);
}

StatePropagatorDataGpu::Impl::Impl(const DeviceStream*  pmeStream,
                                   const DeviceContext& deviceContext,
                                   GpuApiCallBehavior   transferKind,
                                   int                  allocationBlockSizeDivisor,
                                   gmx_wallcycle*       wcycle) :
    deviceContext_(deviceContext),
    transferKind_(transferKind),
    allocationBlockSizeDivisor_(allocationBlockSizeDivisor),
    wcycle_(wcycle)
{
    static_assert(
            GMX_GPU,
            "GPU state propagator data object should only be constructed on the GPU code-paths.");

    GMX_ASSERT(pmeStream->isValid(), "GPU PME stream should be valid.");
    pmeStream_      = pmeStream;
    localStream_    = pmeStream; // For clearing the force buffer
    nonLocalStream_ = nullptr;
    updateStream_   = nullptr;

    // Only local/all coordinates are allowed to be copied in PME-only rank/ PME tests.
    // This it temporary measure to make it safe to use this class in those cases.
    xCopyStreams_[AtomLocality::Local]    = pmeStream_;
    xCopyStreams_[AtomLocality::NonLocal] = nullptr;
    xCopyStreams_[AtomLocality::All]      = nullptr;

    vCopyStreams_[AtomLocality::Local]    = nullptr;
    vCopyStreams_[AtomLocality::NonLocal] = nullptr;
    vCopyStreams_[AtomLocality::All]      = nullptr;

    fCopyStreams_[AtomLocality::Local]    = nullptr;
    fCopyStreams_[AtomLocality::NonLocal] = nullptr;
    fCopyStreams_[AtomLocality::All]      = nullptr;
}

StatePropagatorDataGpu::Impl::~Impl()
{
    // Flush all the streams before freeing memory. See #4519.
    const std::array<const DeviceStream*, 6> allStreams{ pmeStream_,          localStream_,
                                                         nonLocalStream_,     updateStream_,
                                                         copyInStream_.get(), memsetStream_.get() };
    for (const DeviceStream* stream : allStreams)
    {
        if (stream)
        {
            stream->synchronize();
        }
    }

    freeDeviceBuffer(&d_x_);
    freeDeviceBuffer(&d_v_);
    freeDeviceBuffer(&d_f_);
}

void StatePropagatorDataGpu::Impl::reinit(int numAtomsLocal, int numAtomsAll)
{
    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start_nocount(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);

    numAtomsLocal_ = numAtomsLocal;
    numAtomsAll_   = numAtomsAll;

    int numAtomsPadded;
    if (allocationBlockSizeDivisor_ > 0)
    {
        numAtomsPadded = ((numAtomsAll_ + allocationBlockSizeDivisor_ - 1) / allocationBlockSizeDivisor_)
                         * allocationBlockSizeDivisor_;
    }
    else
    {
        numAtomsPadded = numAtomsAll_;
    }

    reallocateDeviceBuffer(&d_x_, numAtomsPadded, &d_xSize_, &d_xCapacity_, deviceContext_);

    const size_t paddingAllocationSize = numAtomsPadded - numAtomsAll_;
    if (paddingAllocationSize > 0)
    {
        // The PME stream is used here because the padding region of d_x_ is only in the PME task.
        clearDeviceBufferAsync(&d_x_, numAtomsAll_, paddingAllocationSize, *pmeStream_);
        // Wait for clearing to complete since with PME-PP pipelining PME will use different streams.
        pmeStream_->synchronize();
    }

    reallocateDeviceBuffer(&d_v_, numAtomsAll_, &d_vSize_, &d_vCapacity_, deviceContext_);
    reallocateDeviceBuffer(&d_f_, numAtomsAll_, &d_fSize_, &d_fCapacity_, deviceContext_);

    // Clearing of the forces can be done in local stream since the nonlocal stream cannot reach
    // the force accumulation stage before syncing with the local stream. Only done in CUDA and
    // SYCL, since the force buffer ops are not implemented in OpenCL.
    static constexpr bool sc_haveGpuFBufferOps = ((GMX_GPU_CUDA != 0) || (GMX_GPU_SYCL != 0));
    if (sc_haveGpuFBufferOps)
    {
        clearDeviceBufferAsync(&d_f_, 0, d_fCapacity_, *localStream_);
        // We need to synchronize to avoid a data race with copyForcesToGpu(Local).
        localStream_->synchronize();
    }

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

std::tuple<int, int> StatePropagatorDataGpu::Impl::getAtomRangesFromAtomLocality(AtomLocality atomLocality) const
{
    int atomsStartAt   = 0;
    int numAtomsToCopy = 0;
    switch (atomLocality)
    {
        case AtomLocality::All:
            atomsStartAt   = 0;
            numAtomsToCopy = numAtomsAll_;
            break;
        case AtomLocality::Local:
            atomsStartAt   = 0;
            numAtomsToCopy = numAtomsLocal_;
            break;
        case AtomLocality::NonLocal:
            atomsStartAt   = numAtomsLocal_;
            numAtomsToCopy = numAtomsAll_ - numAtomsLocal_;
            break;
        default:
            GMX_RELEASE_ASSERT(false,
                               "Wrong range of atoms requested in GPU state data manager. Should "
                               "be All, Local or NonLocal.");
    }
    GMX_ASSERT(atomsStartAt >= 0,
               "The first element to copy has negative index. Probably, the GPU propagator state "
               "was not initialized.");
    GMX_ASSERT(numAtomsToCopy >= 0,
               "Number of atoms to copy is negative. Probably, the GPU propagator state was not "
               "initialized.");
    return std::make_tuple(atomsStartAt, numAtomsToCopy);
}

void StatePropagatorDataGpu::Impl::copyToDevice(DeviceBuffer<RVec>                   d_data,
                                                const gmx::ArrayRef<const gmx::RVec> h_data,
                                                int                                  dataSize,
                                                AtomLocality                         atomLocality,
                                                const DeviceStream&                  deviceStream)
{
    GMX_UNUSED_VALUE(dataSize);

    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");

    GMX_ASSERT(dataSize >= 0, "Trying to copy to device buffer before it was allocated.");

    GMX_ASSERT(deviceStream.isValid(), "No stream is valid for copying with given atom locality.");

    int atomsStartAt, numAtomsToCopy;
    std::tie(atomsStartAt, numAtomsToCopy) = getAtomRangesFromAtomLocality(atomLocality);

    if (numAtomsToCopy != 0)
    {
        GMX_ASSERT(atomsStartAt + numAtomsToCopy <= dataSize,
                   "The device allocation is smaller than requested copy range.");
        GMX_ASSERT(atomsStartAt + numAtomsToCopy <= h_data.ssize(),
                   "The host buffer is smaller than the requested copy range.");

        copyToDeviceBuffer(&d_data,
                           reinterpret_cast<const RVec*>(&h_data.data()[atomsStartAt]),
                           atomsStartAt,
                           numAtomsToCopy,
                           deviceStream,
                           transferKind_,
                           nullptr);
    }
}

void StatePropagatorDataGpu::Impl::copyFromDevice(gmx::ArrayRef<gmx::RVec> h_data,
                                                  DeviceBuffer<RVec>       d_data,
                                                  int                      dataSize,
                                                  AtomLocality             atomLocality,
                                                  const DeviceStream&      deviceStream)
{
    GMX_UNUSED_VALUE(dataSize);

    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");

    GMX_ASSERT(dataSize >= 0, "Trying to copy from device buffer before it was allocated.");

    GMX_ASSERT(deviceStream.isValid(), "No stream is valid for copying with given atom locality.");

    int atomsStartAt, numAtomsToCopy;
    std::tie(atomsStartAt, numAtomsToCopy) = getAtomRangesFromAtomLocality(atomLocality);

    if (numAtomsToCopy != 0)
    {
        GMX_ASSERT(atomsStartAt + numAtomsToCopy <= dataSize,
                   "The device allocation is smaller than requested copy range.");
        GMX_ASSERT(atomsStartAt + numAtomsToCopy <= h_data.ssize(),
                   "The host buffer is smaller than the requested copy range.");

        copyFromDeviceBuffer(reinterpret_cast<RVec*>(&h_data.data()[atomsStartAt]),
                             &d_data,
                             atomsStartAt,
                             numAtomsToCopy,
                             deviceStream,
                             transferKind_,
                             nullptr);
    }
}

void StatePropagatorDataGpu::Impl::clearOnDevice(DeviceBuffer<RVec>  d_data,
                                                 int                 dataSize,
                                                 AtomLocality        atomLocality,
                                                 const DeviceStream& deviceStream) const
{
    GMX_UNUSED_VALUE(dataSize);

    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");

    GMX_ASSERT(dataSize >= 0, "Trying to clear to device buffer before it was allocated.");

    GMX_ASSERT(deviceStream.isValid(), "No stream is valid for clearing with given atom locality.");

    int atomsStartAt, numAtomsToClear;
    std::tie(atomsStartAt, numAtomsToClear) = getAtomRangesFromAtomLocality(atomLocality);

    if (numAtomsToClear != 0)
    {
        GMX_ASSERT(atomsStartAt + numAtomsToClear <= dataSize,
                   "The device allocation is smaller than requested clear range.");

        clearDeviceBufferAsync(&d_data, atomsStartAt, numAtomsToClear, deviceStream);
    }
}

DeviceBuffer<RVec> StatePropagatorDataGpu::Impl::getCoordinates()
{
    return d_x_;
}

void StatePropagatorDataGpu::Impl::copyCoordinatesToGpu(const gmx::ArrayRef<const gmx::RVec> h_x,
                                                        AtomLocality atomLocality,
                                                        int          expectedConsumptionCount)
{
    GMX_ASSERT(atomLocality < AtomLocality::All,
               formatString("Wrong atom locality. Only Local and NonLocal are allowed for "
                            "coordinate transfers, passed value is \"%s\"",
                            enumValueToString(atomLocality))
                       .c_str());

    const DeviceStream* deviceStream = xCopyStreams_[atomLocality];
    GMX_ASSERT(deviceStream != nullptr,
               "No stream is valid for copying positions with given atom locality.");

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);

    copyToDevice(d_x_, h_x, d_xSize_, atomLocality, *deviceStream);

    if (expectedConsumptionCount > 0)
    {
        xReadyOnDevice_[atomLocality].markEvent(*deviceStream);
    }
    xReadyOnDevice_[atomLocality].setConsumptionLimits(expectedConsumptionCount, expectedConsumptionCount);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

GpuEventSynchronizer* StatePropagatorDataGpu::Impl::getCoordinatesReadyOnDeviceEvent(
        AtomLocality              atomLocality,
        const SimulationWorkload& simulationWork,
        const StepWorkload&       stepWork,
        GpuEventSynchronizer*     gpuCoordinateHaloLaunched)
{
    // The provider of the coordinates may be different for local atoms. If the update is offloaded
    // and this is not a neighbor search step, then the consumer needs to wait for the update
    // to complete. Otherwise, the coordinates are copied from the host and we need to wait for
    // the copy event. Non-local coordinates are provided by the GPU halo exchange (if active), otherwise by H2D copy.

    if (atomLocality == AtomLocality::NonLocal && stepWork.useGpuXHalo)
    {
        GMX_ASSERT(gpuCoordinateHaloLaunched != nullptr,
                   "GPU halo exchange is active but its completion event is null.");
        return gpuCoordinateHaloLaunched;
    }
    if (atomLocality == AtomLocality::Local && simulationWork.useGpuUpdate && !stepWork.doNeighborSearch)
    {
        GMX_ASSERT(xUpdatedOnDeviceEvent_ != nullptr, "The event synchronizer can not be nullptr.");
        return xUpdatedOnDeviceEvent_;
    }
    else
    {
        if (stepWork.doNeighborSearch && xUpdatedOnDeviceEvent_)
        {
            /* On search steps, we do not consume the result of the GPU update
             * but rather that of a H2D transfer. So, we reset the event triggered after
             * update to avoid leaving it unconsumed.
             * Unfortunately, we don't always have the event marked either (e.g., on the
             * first step) so we just reset it here.
             * See Issue #3988. */
            xUpdatedOnDeviceEvent_->reset();
        }
        return &xReadyOnDevice_[atomLocality];
    }
}

void StatePropagatorDataGpu::Impl::waitCoordinatesCopiedToDevice(AtomLocality atomLocality)
{
    wallcycle_start(wcycle_, WallCycleCounter::WaitGpuStatePropagatorData);
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    xReadyOnDevice_[atomLocality].waitForEvent();
    wallcycle_stop(wcycle_, WallCycleCounter::WaitGpuStatePropagatorData);
}

void StatePropagatorDataGpu::Impl::consumeCoordinatesCopiedToDeviceEvent(AtomLocality atomLocality)
{
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    xReadyOnDevice_[atomLocality].consume();
}

void StatePropagatorDataGpu::Impl::resetCoordinatesCopiedToDeviceEvent(AtomLocality atomLocality)
{
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    xReadyOnDevice_[atomLocality].reset();
}

void StatePropagatorDataGpu::Impl::setXUpdatedOnDeviceEvent(GpuEventSynchronizer* xUpdatedOnDeviceEvent)
{
    GMX_ASSERT(xUpdatedOnDeviceEvent != nullptr, "The event synchronizer can not be nullptr.");
    xUpdatedOnDeviceEvent_ = xUpdatedOnDeviceEvent;
}

void StatePropagatorDataGpu::Impl::setXUpdatedOnDeviceEventExpectedConsumptionCount(int expectedConsumptionCount)
{
    xUpdatedOnDeviceEvent_->setConsumptionLimits(expectedConsumptionCount, expectedConsumptionCount);
}

void StatePropagatorDataGpu::Impl::setFReadyOnDeviceEventExpectedConsumptionCount(AtomLocality atomLocality,
                                                                                  int expectedConsumptionCount)
{
    fReadyOnDevice_[atomLocality].setConsumptionLimits(expectedConsumptionCount, expectedConsumptionCount);
}

void StatePropagatorDataGpu::Impl::copyCoordinatesFromGpu(gmx::ArrayRef<gmx::RVec> h_x,
                                                          AtomLocality             atomLocality,
                                                          GpuEventSynchronizer*    dependency)
{
    GMX_ASSERT(atomLocality < AtomLocality::All,
               formatString("Wrong atom locality. Only Local and NonLocal are allowed for "
                            "coordinate transfers, passed value is \"%s\"",
                            enumValueToString(atomLocality))
                       .c_str());
    const DeviceStream* deviceStream = xCopyStreams_[atomLocality];
    GMX_ASSERT(deviceStream != nullptr,
               "No stream is valid for copying positions with given atom locality.");

    if (dependency != nullptr)
    {
        dependency->enqueueWaitEvent(*deviceStream);
    }

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);

    copyFromDevice(h_x, d_x_, d_xSize_, atomLocality, *deviceStream);
    // Note: unlike copyCoordinatesToGpu this is not used in OpenCL, and the conditional is not needed.
    xReadyOnHost_[atomLocality].markEvent(*deviceStream);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

void StatePropagatorDataGpu::Impl::waitCoordinatesReadyOnHost(AtomLocality atomLocality)
{
    wallcycle_start(wcycle_, WallCycleCounter::WaitGpuStatePropagatorData);
    xReadyOnHost_[atomLocality].waitForEvent();
    wallcycle_stop(wcycle_, WallCycleCounter::WaitGpuStatePropagatorData);
}


DeviceBuffer<RVec> StatePropagatorDataGpu::Impl::getVelocities()
{
    return d_v_;
}

void StatePropagatorDataGpu::Impl::copyVelocitiesToGpu(const gmx::ArrayRef<const gmx::RVec> h_v,
                                                       AtomLocality atomLocality)
{
    GMX_ASSERT(atomLocality == AtomLocality::Local,
               formatString("Wrong atom locality. Only Local is allowed for "
                            "velocity transfers, passed value is \"%s\"",
                            enumValueToString(atomLocality))
                       .c_str());
    const DeviceStream* deviceStream = vCopyStreams_[atomLocality];
    GMX_ASSERT(deviceStream != nullptr,
               "No stream is valid for copying velocities with given atom locality.");

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);

    copyToDevice(d_v_, h_v, d_vSize_, atomLocality, *deviceStream);
    /* Not marking the event, because it is not used anywhere.
     * Since we only use velocities on the device for update, and we launch the copy in
     * the "update" stream, that should be safe.
     */

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

void StatePropagatorDataGpu::Impl::copyVelocitiesFromGpu(gmx::ArrayRef<gmx::RVec> h_v, AtomLocality atomLocality)
{
    GMX_ASSERT(atomLocality == AtomLocality::Local,
               formatString("Wrong atom locality. Only Local is allowed for "
                            "velocity transfers, passed value is \"%s\"",
                            enumValueToString(atomLocality))
                       .c_str());
    const DeviceStream* deviceStream = vCopyStreams_[atomLocality];
    GMX_ASSERT(deviceStream != nullptr,
               "No stream is valid for copying velocities with given atom locality.");

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);

    copyFromDevice(h_v, d_v_, d_vSize_, atomLocality, *deviceStream);
    vReadyOnHost_[atomLocality].markEvent(*deviceStream);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

void StatePropagatorDataGpu::Impl::waitVelocitiesReadyOnHost(AtomLocality atomLocality)
{
    wallcycle_start(wcycle_, WallCycleCounter::WaitGpuStatePropagatorData);
    vReadyOnHost_[atomLocality].waitForEvent();
    wallcycle_stop(wcycle_, WallCycleCounter::WaitGpuStatePropagatorData);
}


DeviceBuffer<RVec> StatePropagatorDataGpu::Impl::getForces()
{
    return d_f_;
}

// Copy CPU forces to GPU using stream internal to this module to allow overlap
// with GPU force calculations.
void StatePropagatorDataGpu::Impl::copyForcesToGpu(const gmx::ArrayRef<const gmx::RVec> h_f,
                                                   AtomLocality atomLocality)
{
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    DeviceStream* deviceStream = copyInStream_.get();
    GMX_ASSERT(deviceStream != nullptr,
               "No stream is valid for copying forces with given atom locality.");

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);

    copyToDevice(d_f_, h_f, d_fSize_, atomLocality, *deviceStream);
    fReadyOnDevice_[atomLocality].markEvent(*deviceStream);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

void StatePropagatorDataGpu::Impl::clearForcesOnGpu(AtomLocality atomLocality, GpuEventSynchronizer* dependency)
{
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    DeviceStream* deviceStream = memsetStream_.get();

    if (dependency != nullptr)
    {
        dependency->enqueueWaitEvent(*deviceStream);
    }

    GMX_ASSERT(deviceStream != nullptr,
               "No stream is valid for clearing forces with given atom locality.");

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);

    clearOnDevice(d_f_, d_fSize_, atomLocality, *deviceStream);

    fReadyOnDevice_[atomLocality].markEvent(*deviceStream);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

GpuEventSynchronizer* StatePropagatorDataGpu::Impl::getLocalForcesReadyOnDeviceEvent(StepWorkload stepWork,
                                                                                     SimulationWorkload simulationWork)
{
    if (stepWork.useGpuFBufferOps && !simulationWork.useCpuPmePpCommunication)
    {
        return &fReducedOnDevice_[AtomLocality::Local];
    }
    else
    {
        return &fReadyOnDevice_[AtomLocality::Local];
    }
}

GpuEventSynchronizer* StatePropagatorDataGpu::Impl::fReducedOnDevice(AtomLocality atomLocality)
{
    return &fReducedOnDevice_[atomLocality];
}

void StatePropagatorDataGpu::Impl::consumeForcesReducedOnDeviceEvent(AtomLocality atomLocality)
{
    fReducedOnDevice_[atomLocality].consume();
}

GpuEventSynchronizer* StatePropagatorDataGpu::Impl::fReadyOnDevice(AtomLocality atomLocality)
{
    return &fReadyOnDevice_[atomLocality];
}

void StatePropagatorDataGpu::Impl::copyForcesFromGpu(gmx::ArrayRef<gmx::RVec> h_f, AtomLocality atomLocality)
{
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    const DeviceStream* deviceStream = fCopyStreams_[atomLocality];
    GMX_ASSERT(deviceStream != nullptr,
               "No stream is valid for copying forces with given atom locality.");

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);

    copyFromDevice(h_f, d_f_, d_fSize_, atomLocality, *deviceStream);
    fReadyOnHost_[atomLocality].markEvent(*deviceStream);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchStatePropagatorData);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

void StatePropagatorDataGpu::Impl::waitForcesReadyOnHost(AtomLocality atomLocality)
{
    wallcycle_start(wcycle_, WallCycleCounter::WaitGpuStatePropagatorData);
    fReadyOnHost_[atomLocality].waitForEvent();
    wallcycle_stop(wcycle_, WallCycleCounter::WaitGpuStatePropagatorData);
}

const DeviceStream* StatePropagatorDataGpu::Impl::getUpdateStream()
{
    return updateStream_;
}

int StatePropagatorDataGpu::Impl::numAtomsLocal() const
{
    return numAtomsLocal_;
}

int StatePropagatorDataGpu::Impl::numAtomsAll() const
{
    return numAtomsAll_;
}


StatePropagatorDataGpu::StatePropagatorDataGpu(const DeviceStreamManager& deviceStreamManager,
                                               GpuApiCallBehavior         transferKind,
                                               int            allocationBlockSizeDivisor,
                                               gmx_wallcycle* wcycle) :
    impl_(new Impl(deviceStreamManager, transferKind, allocationBlockSizeDivisor, wcycle))
{
}

StatePropagatorDataGpu::StatePropagatorDataGpu(const DeviceStream*  pmeStream,
                                               const DeviceContext& deviceContext,
                                               GpuApiCallBehavior   transferKind,
                                               int                  allocationBlockSizeDivisor,
                                               gmx_wallcycle*       wcycle) :
    impl_(new Impl(pmeStream, deviceContext, transferKind, allocationBlockSizeDivisor, wcycle))
{
}

StatePropagatorDataGpu::StatePropagatorDataGpu(StatePropagatorDataGpu&& /* other */) noexcept = default;

StatePropagatorDataGpu& StatePropagatorDataGpu::operator=(StatePropagatorDataGpu&& /* other */) noexcept = default;

StatePropagatorDataGpu::~StatePropagatorDataGpu() = default;


void StatePropagatorDataGpu::reinit(int numAtomsLocal, int numAtomsAll)
{
    return impl_->reinit(numAtomsLocal, numAtomsAll);
}

std::tuple<int, int> StatePropagatorDataGpu::getAtomRangesFromAtomLocality(AtomLocality atomLocality) const
{
    return impl_->getAtomRangesFromAtomLocality(atomLocality);
}


DeviceBuffer<RVec> StatePropagatorDataGpu::getCoordinates()
{
    return impl_->getCoordinates();
}

void StatePropagatorDataGpu::copyCoordinatesToGpu(const gmx::ArrayRef<const gmx::RVec> h_x,
                                                  AtomLocality                         atomLocality,
                                                  int expectedConsumptionCount)
{
    return impl_->copyCoordinatesToGpu(h_x, atomLocality, expectedConsumptionCount);
}

GpuEventSynchronizer*
StatePropagatorDataGpu::getCoordinatesReadyOnDeviceEvent(AtomLocality              atomLocality,
                                                         const SimulationWorkload& simulationWork,
                                                         const StepWorkload&       stepWork,
                                                         GpuEventSynchronizer* gpuCoordinateHaloLaunched)
{
    return impl_->getCoordinatesReadyOnDeviceEvent(
            atomLocality, simulationWork, stepWork, gpuCoordinateHaloLaunched);
}

void StatePropagatorDataGpu::waitCoordinatesCopiedToDevice(AtomLocality atomLocality)
{
    return impl_->waitCoordinatesCopiedToDevice(atomLocality);
}

void StatePropagatorDataGpu::consumeCoordinatesCopiedToDeviceEvent(AtomLocality atomLocality)
{
    return impl_->consumeCoordinatesCopiedToDeviceEvent(atomLocality);
}

void StatePropagatorDataGpu::resetCoordinatesCopiedToDeviceEvent(AtomLocality atomLocality)
{
    return impl_->resetCoordinatesCopiedToDeviceEvent(atomLocality);
}

void StatePropagatorDataGpu::setXUpdatedOnDeviceEvent(GpuEventSynchronizer* xUpdatedOnDeviceEvent)
{
    impl_->setXUpdatedOnDeviceEvent(xUpdatedOnDeviceEvent);
}

void StatePropagatorDataGpu::setXUpdatedOnDeviceEventExpectedConsumptionCount(int expectedConsumptionCount)
{
    impl_->setXUpdatedOnDeviceEventExpectedConsumptionCount(expectedConsumptionCount);
}

void StatePropagatorDataGpu::setFReadyOnDeviceEventExpectedConsumptionCount(AtomLocality atomLocality,
                                                                            int expectedConsumptionCount)
{
    impl_->setFReadyOnDeviceEventExpectedConsumptionCount(atomLocality, expectedConsumptionCount);
}

void StatePropagatorDataGpu::copyCoordinatesFromGpu(gmx::ArrayRef<RVec>   h_x,
                                                    AtomLocality          atomLocality,
                                                    GpuEventSynchronizer* dependency)
{
    return impl_->copyCoordinatesFromGpu(h_x, atomLocality, dependency);
}

void StatePropagatorDataGpu::waitCoordinatesReadyOnHost(AtomLocality atomLocality)
{
    return impl_->waitCoordinatesReadyOnHost(atomLocality);
}


DeviceBuffer<RVec> StatePropagatorDataGpu::getVelocities()
{
    return impl_->getVelocities();
}

void StatePropagatorDataGpu::copyVelocitiesToGpu(const gmx::ArrayRef<const gmx::RVec> h_v,
                                                 AtomLocality                         atomLocality)
{
    return impl_->copyVelocitiesToGpu(h_v, atomLocality);
}

void StatePropagatorDataGpu::copyVelocitiesFromGpu(gmx::ArrayRef<RVec> h_v, AtomLocality atomLocality)
{
    return impl_->copyVelocitiesFromGpu(h_v, atomLocality);
}

void StatePropagatorDataGpu::waitVelocitiesReadyOnHost(AtomLocality atomLocality)
{
    return impl_->waitVelocitiesReadyOnHost(atomLocality);
}


DeviceBuffer<RVec> StatePropagatorDataGpu::getForces()
{
    return impl_->getForces();
}

void StatePropagatorDataGpu::copyForcesToGpu(const gmx::ArrayRef<const gmx::RVec> h_f, AtomLocality atomLocality)
{
    return impl_->copyForcesToGpu(h_f, atomLocality);
}

void StatePropagatorDataGpu::clearForcesOnGpu(AtomLocality atomLocality, GpuEventSynchronizer* dependency)
{
    return impl_->clearForcesOnGpu(atomLocality, dependency);
}

GpuEventSynchronizer* StatePropagatorDataGpu::getLocalForcesReadyOnDeviceEvent(StepWorkload stepWork,
                                                                               SimulationWorkload simulationWork)
{
    return impl_->getLocalForcesReadyOnDeviceEvent(stepWork, simulationWork);
}

GpuEventSynchronizer* StatePropagatorDataGpu::fReducedOnDevice(AtomLocality atomLocality)
{
    return impl_->fReducedOnDevice(atomLocality);
}

void StatePropagatorDataGpu::consumeForcesReducedOnDeviceEvent(AtomLocality atomLocality)
{
    impl_->consumeForcesReducedOnDeviceEvent(atomLocality);
}

GpuEventSynchronizer* StatePropagatorDataGpu::fReadyOnDevice(AtomLocality atomLocality)
{
    return impl_->fReadyOnDevice(atomLocality);
}

void StatePropagatorDataGpu::copyForcesFromGpu(gmx::ArrayRef<RVec> h_f, AtomLocality atomLocality)
{
    return impl_->copyForcesFromGpu(h_f, atomLocality);
}

void StatePropagatorDataGpu::waitForcesReadyOnHost(AtomLocality atomLocality)
{
    return impl_->waitForcesReadyOnHost(atomLocality);
}


const DeviceStream* StatePropagatorDataGpu::getUpdateStream()
{
    return impl_->getUpdateStream();
}

int StatePropagatorDataGpu::numAtomsLocal() const
{
    return impl_->numAtomsLocal();
}

int StatePropagatorDataGpu::numAtomsAll() const
{
    return impl_->numAtomsAll();
}

void StatePropagatorDataGpu::waitCoordinatesUpdatedOnDevice()
{
    impl_->waitCoordinatesUpdatedOnDevice();
}

} // namespace gmx

#endif // GMX_GPU
