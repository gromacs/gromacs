/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Definitions of interfaces for GPU state data propagator object.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include "config.h"

#if GMX_GPU != GMX_GPU_NONE

#    if GMX_GPU == GMX_GPU_CUDA
#        include "gromacs/gpu_utils/cudautils.cuh"
#    endif
#    include "gromacs/gpu_utils/devicebuffer.h"
#    if GMX_GPU == GMX_GPU_OPENCL
#        include "gromacs/gpu_utils/oclutils.h"
#    endif
#    include "gromacs/math/vectypes.h"
#    include "gromacs/mdtypes/state_propagator_data_gpu.h"
#    include "gromacs/timing/wallcycle.h"
#    include "gromacs/utility/classhelpers.h"

#    include "state_propagator_data_gpu_impl.h"


namespace gmx
{

StatePropagatorDataGpu::Impl::Impl(const void*        pmeStream,
                                   const void*        localStream,
                                   const void*        nonLocalStream,
                                   const void*        deviceContext,
                                   GpuApiCallBehavior transferKind,
                                   int                paddingSize,
                                   gmx_wallcycle*     wcycle) :
    transferKind_(transferKind),
    paddingSize_(paddingSize),
    wcycle_(wcycle)
{
    static_assert(GMX_GPU != GMX_GPU_NONE,
                  "This object should only be constructed on the GPU code-paths.");

    // TODO: Refactor when the StreamManager is introduced.
    if (GMX_GPU == GMX_GPU_OPENCL)
    {
        GMX_ASSERT(deviceContext != nullptr, "GPU context should be set in OpenCL builds.");
        GMX_ASSERT(pmeStream != nullptr, "GPU PME stream should be set in OpenCL builds.");

        // The update stream is set to the PME stream in OpenCL, since PME stream is the only stream created in the PME context.
        pmeStream_     = *static_cast<const CommandStream*>(pmeStream);
        updateStream_  = *static_cast<const CommandStream*>(pmeStream);
        deviceContext_ = *static_cast<const DeviceContext*>(deviceContext);
        GMX_UNUSED_VALUE(localStream);
        GMX_UNUSED_VALUE(nonLocalStream);
    }

    if (GMX_GPU == GMX_GPU_CUDA)
    {
        if (pmeStream != nullptr)
        {
            pmeStream_ = *static_cast<const CommandStream*>(pmeStream);
        }
        if (localStream != nullptr)
        {
            localStream_ = *static_cast<const CommandStream*>(localStream);
        }
        if (nonLocalStream != nullptr)
        {
            nonLocalStream_ = *static_cast<const CommandStream*>(nonLocalStream);
        }

        // TODO: The update stream should be created only when it is needed.
#    if (GMX_GPU == GMX_GPU_CUDA)
        cudaError_t stat;
        stat = cudaStreamCreate(&updateStream_);
        CU_RET_ERR(stat, "CUDA stream creation failed in StatePropagatorDataGpu");
#    endif
        GMX_UNUSED_VALUE(deviceContext);
    }

    // Map the atom locality to the stream that will be used for coordinates,
    // velocities and forces transfers. Same streams are used for H2D and D2H copies.
    // Note, that nullptr stream is used here to indicate that the copy is not supported.
    xCopyStreams_[AtomLocality::Local]    = updateStream_;
    xCopyStreams_[AtomLocality::NonLocal] = nonLocalStream_;
    xCopyStreams_[AtomLocality::All]      = updateStream_;

    vCopyStreams_[AtomLocality::Local]    = updateStream_;
    vCopyStreams_[AtomLocality::NonLocal] = nullptr;
    vCopyStreams_[AtomLocality::All]      = updateStream_;

    fCopyStreams_[AtomLocality::Local]    = localStream_;
    fCopyStreams_[AtomLocality::NonLocal] = nonLocalStream_;
    fCopyStreams_[AtomLocality::All]      = updateStream_;
}

StatePropagatorDataGpu::Impl::Impl(const void*        pmeStream,
                                   const void*        deviceContext,
                                   GpuApiCallBehavior transferKind,
                                   int                paddingSize,
                                   gmx_wallcycle*     wcycle) :
    transferKind_(transferKind),
    paddingSize_(paddingSize),
    wcycle_(wcycle)
{
    static_assert(GMX_GPU != GMX_GPU_NONE,
                  "This object should only be constructed on the GPU code-paths.");

    if (GMX_GPU == GMX_GPU_OPENCL)
    {
        GMX_ASSERT(deviceContext != nullptr, "GPU context should be set in OpenCL builds.");
        deviceContext_ = *static_cast<const DeviceContext*>(deviceContext);
    }

    GMX_ASSERT(pmeStream != nullptr, "GPU PME stream should be set.");
    pmeStream_ = *static_cast<const CommandStream*>(pmeStream);

    localStream_    = nullptr;
    nonLocalStream_ = nullptr;
    updateStream_   = nullptr;


    // Only local/all coordinates are allowed to be copied in PME-only rank/ PME tests.
    // This it temporary measure to make it safe to use this class in those cases.
    xCopyStreams_[AtomLocality::Local]    = pmeStream_;
    xCopyStreams_[AtomLocality::NonLocal] = nullptr;
    xCopyStreams_[AtomLocality::All]      = pmeStream_;

    vCopyStreams_[AtomLocality::Local]    = nullptr;
    vCopyStreams_[AtomLocality::NonLocal] = nullptr;
    vCopyStreams_[AtomLocality::All]      = nullptr;

    fCopyStreams_[AtomLocality::Local]    = nullptr;
    fCopyStreams_[AtomLocality::NonLocal] = nullptr;
    fCopyStreams_[AtomLocality::All]      = nullptr;
}

StatePropagatorDataGpu::Impl::~Impl() {}

void StatePropagatorDataGpu::Impl::reinit(int numAtomsLocal, int numAtomsAll)
{
    wallcycle_start_nocount(wcycle_, ewcLAUNCH_GPU);
    wallcycle_sub_start_nocount(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);

    numAtomsLocal_ = numAtomsLocal;
    numAtomsAll_   = numAtomsAll;

    int numAtomsPadded;
    if (paddingSize_ > 0)
    {
        numAtomsPadded = ((numAtomsAll_ + paddingSize_ - 1) / paddingSize_) * paddingSize_;
    }
    else
    {
        numAtomsPadded = numAtomsAll_;
    }

    reallocateDeviceBuffer(&d_x_, DIM * numAtomsPadded, &d_xSize_, &d_xCapacity_, deviceContext_);

    const size_t paddingAllocationSize = numAtomsPadded - numAtomsAll_;
    if (paddingAllocationSize > 0)
    {
        // The PME stream is used here because the padding region of d_x_ is only in the PME task.
        clearDeviceBufferAsync(&d_x_, DIM * numAtomsAll_, DIM * paddingAllocationSize, pmeStream_);
    }

    reallocateDeviceBuffer(&d_v_, DIM * numAtomsAll_, &d_vSize_, &d_vCapacity_, deviceContext_);
    const int d_fOldCapacity = d_fCapacity_;
    reallocateDeviceBuffer(&d_f_, DIM * numAtomsAll_, &d_fSize_, &d_fCapacity_, deviceContext_);
    // Clearing of the forces can be done in local stream since the nonlocal stream cannot reach
    // the force accumulation stage before syncing with the local stream. Only done in CUDA,
    // since the force buffer ops are not implemented in OpenCL.
    if (GMX_GPU == GMX_GPU_CUDA && d_fCapacity_ != d_fOldCapacity)
    {
        clearDeviceBufferAsync(&d_f_, 0, d_fCapacity_, localStream_);
    }

    wallcycle_sub_stop(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);
    wallcycle_stop(wcycle_, ewcLAUNCH_GPU);
}

std::tuple<int, int> StatePropagatorDataGpu::Impl::getAtomRangesFromAtomLocality(AtomLocality atomLocality)
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
               "The first elemtnt to copy has negative index. Probably, the GPU propagator state "
               "was not initialized.");
    GMX_ASSERT(numAtomsToCopy >= 0,
               "Number of atoms to copy is negative. Probably, the GPU propagator state was not "
               "initialized.");
    return std::make_tuple(atomsStartAt, numAtomsToCopy);
}

void StatePropagatorDataGpu::Impl::copyToDevice(DeviceBuffer<float>                  d_data,
                                                const gmx::ArrayRef<const gmx::RVec> h_data,
                                                int                                  dataSize,
                                                AtomLocality                         atomLocality,
                                                CommandStream                        commandStream)
{
    GMX_UNUSED_VALUE(dataSize);

    GMX_ASSERT(dataSize >= 0, "Trying to copy to device buffer before it was allocated.");

    wallcycle_start_nocount(wcycle_, ewcLAUNCH_GPU);
    wallcycle_sub_start(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);

    int atomsStartAt, numAtomsToCopy;
    std::tie(atomsStartAt, numAtomsToCopy) = getAtomRangesFromAtomLocality(atomLocality);

    int elementsStartAt   = atomsStartAt * DIM;
    int numElementsToCopy = numAtomsToCopy * DIM;

    if (numAtomsToCopy != 0)
    {
        GMX_ASSERT(elementsStartAt + numElementsToCopy <= dataSize,
                   "The device allocation is smaller than requested copy range.");
        GMX_ASSERT(atomsStartAt + numAtomsToCopy <= h_data.ssize(),
                   "The host buffer is smaller than the requested copy range.");

        copyToDeviceBuffer(&d_data, reinterpret_cast<const float*>(&h_data.data()[atomsStartAt]),
                           elementsStartAt, numElementsToCopy, commandStream, transferKind_, nullptr);
    }

    wallcycle_sub_stop(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);
    wallcycle_stop(wcycle_, ewcLAUNCH_GPU);
}

void StatePropagatorDataGpu::Impl::copyFromDevice(gmx::ArrayRef<gmx::RVec> h_data,
                                                  DeviceBuffer<float>      d_data,
                                                  int                      dataSize,
                                                  AtomLocality             atomLocality,
                                                  CommandStream            commandStream)
{
    GMX_UNUSED_VALUE(dataSize);

    GMX_ASSERT(dataSize >= 0, "Trying to copy from device buffer before it was allocated.");

    wallcycle_start_nocount(wcycle_, ewcLAUNCH_GPU);
    wallcycle_sub_start(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);

    int atomsStartAt, numAtomsToCopy;
    std::tie(atomsStartAt, numAtomsToCopy) = getAtomRangesFromAtomLocality(atomLocality);

    int elementsStartAt   = atomsStartAt * DIM;
    int numElementsToCopy = numAtomsToCopy * DIM;

    if (numAtomsToCopy != 0)
    {
        GMX_ASSERT(elementsStartAt + numElementsToCopy <= dataSize,
                   "The device allocation is smaller than requested copy range.");
        GMX_ASSERT(atomsStartAt + numAtomsToCopy <= h_data.ssize(),
                   "The host buffer is smaller than the requested copy range.");

        copyFromDeviceBuffer(reinterpret_cast<float*>(&h_data.data()[atomsStartAt]), &d_data,
                             elementsStartAt, numElementsToCopy, commandStream, transferKind_, nullptr);
    }

    wallcycle_sub_stop(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);
    wallcycle_stop(wcycle_, ewcLAUNCH_GPU);
}

DeviceBuffer<float> StatePropagatorDataGpu::Impl::getCoordinates()
{
    return d_x_;
}

void StatePropagatorDataGpu::Impl::copyCoordinatesToGpu(const gmx::ArrayRef<const gmx::RVec> h_x,
                                                        AtomLocality atomLocality)
{
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    CommandStream commandStream = xCopyStreams_[atomLocality];
    GMX_ASSERT(commandStream != nullptr,
               "No stream is valid for copying positions with given atom locality.");

    wallcycle_start_nocount(wcycle_, ewcLAUNCH_GPU);
    wallcycle_sub_start(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);

    copyToDevice(d_x_, h_x, d_xSize_, atomLocality, commandStream);

    // markEvent is skipped in OpenCL as:
    //   - it's not needed, copy is done in the same stream as the only consumer task (PME)
    //   - we don't consume the events in OpenCL which is not allowed by GpuEventSynchronizer (would leak memory).
    // TODO: remove this by adding an event-mark free flavor of this function
    if (GMX_GPU == GMX_GPU_CUDA)
    {
        xReadyOnDevice_[atomLocality].markEvent(commandStream);
    }

    wallcycle_sub_stop(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);
    wallcycle_stop(wcycle_, ewcLAUNCH_GPU);
}

GpuEventSynchronizer*
StatePropagatorDataGpu::Impl::getCoordinatesReadyOnDeviceEvent(AtomLocality atomLocality,
                                                               const SimulationWorkload& simulationWork,
                                                               const StepWorkload&       stepWork)
{
    // The provider of the coordinates may be different for local atoms. If the update is offloaded
    // and this is not a neighbor search step, then the consumer needs to wait for the update
    // to complete. Otherwise, the coordinates are copied from the host and we need to wait for
    // the copy event. Non-local coordinates are always provided by the H2D copy.
    //
    // TODO: This should be reconsidered to support the halo exchange.
    //
    // In OpenCL no events are used as coordinate sync is not necessary
    if (GMX_GPU == GMX_GPU_OPENCL)
    {
        return nullptr;
    }
    if (atomLocality == AtomLocality::Local && simulationWork.useGpuUpdate && !stepWork.doNeighborSearch)
    {
        return &xUpdatedOnDevice_;
    }
    else
    {
        return &xReadyOnDevice_[atomLocality];
    }
}

void StatePropagatorDataGpu::Impl::waitCoordinatesCopiedToDevice(AtomLocality atomLocality)
{
    wallcycle_start(wcycle_, ewcWAIT_GPU_STATE_PROPAGATOR_DATA);
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    xReadyOnDevice_[atomLocality].waitForEvent();
    wallcycle_stop(wcycle_, ewcWAIT_GPU_STATE_PROPAGATOR_DATA);
}

GpuEventSynchronizer* StatePropagatorDataGpu::Impl::xUpdatedOnDevice()
{
    return &xUpdatedOnDevice_;
}

void StatePropagatorDataGpu::Impl::copyCoordinatesFromGpu(gmx::ArrayRef<gmx::RVec> h_x, AtomLocality atomLocality)
{
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    CommandStream commandStream = xCopyStreams_[atomLocality];
    GMX_ASSERT(commandStream != nullptr,
               "No stream is valid for copying positions with given atom locality.");

    wallcycle_start_nocount(wcycle_, ewcLAUNCH_GPU);
    wallcycle_sub_start(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);

    copyFromDevice(h_x, d_x_, d_xSize_, atomLocality, commandStream);
    // Note: unlike copyCoordinatesToGpu this is not used in OpenCL, and the conditional is not needed.
    xReadyOnHost_[atomLocality].markEvent(commandStream);

    wallcycle_sub_stop(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);
    wallcycle_stop(wcycle_, ewcLAUNCH_GPU);
}

void StatePropagatorDataGpu::Impl::waitCoordinatesReadyOnHost(AtomLocality atomLocality)
{
    wallcycle_start(wcycle_, ewcWAIT_GPU_STATE_PROPAGATOR_DATA);
    xReadyOnHost_[atomLocality].waitForEvent();
    wallcycle_stop(wcycle_, ewcWAIT_GPU_STATE_PROPAGATOR_DATA);
}


DeviceBuffer<float> StatePropagatorDataGpu::Impl::getVelocities()
{
    return d_v_;
}

void StatePropagatorDataGpu::Impl::copyVelocitiesToGpu(const gmx::ArrayRef<const gmx::RVec> h_v,
                                                       AtomLocality atomLocality)
{
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    CommandStream commandStream = vCopyStreams_[atomLocality];
    GMX_ASSERT(commandStream != nullptr,
               "No stream is valid for copying velocities with given atom locality.");

    wallcycle_start_nocount(wcycle_, ewcLAUNCH_GPU);
    wallcycle_sub_start(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);

    copyToDevice(d_v_, h_v, d_vSize_, atomLocality, commandStream);
    vReadyOnDevice_[atomLocality].markEvent(commandStream);

    wallcycle_sub_stop(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);
    wallcycle_stop(wcycle_, ewcLAUNCH_GPU);
}

GpuEventSynchronizer* StatePropagatorDataGpu::Impl::getVelocitiesReadyOnDeviceEvent(AtomLocality atomLocality)
{
    return &vReadyOnDevice_[atomLocality];
}


void StatePropagatorDataGpu::Impl::copyVelocitiesFromGpu(gmx::ArrayRef<gmx::RVec> h_v, AtomLocality atomLocality)
{
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    CommandStream commandStream = vCopyStreams_[atomLocality];
    GMX_ASSERT(commandStream != nullptr,
               "No stream is valid for copying velocities with given atom locality.");

    wallcycle_start_nocount(wcycle_, ewcLAUNCH_GPU);
    wallcycle_sub_start(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);

    copyFromDevice(h_v, d_v_, d_vSize_, atomLocality, commandStream);
    vReadyOnHost_[atomLocality].markEvent(commandStream);

    wallcycle_sub_stop(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);
    wallcycle_stop(wcycle_, ewcLAUNCH_GPU);
}

void StatePropagatorDataGpu::Impl::waitVelocitiesReadyOnHost(AtomLocality atomLocality)
{
    wallcycle_start(wcycle_, ewcWAIT_GPU_STATE_PROPAGATOR_DATA);
    vReadyOnHost_[atomLocality].waitForEvent();
    wallcycle_stop(wcycle_, ewcWAIT_GPU_STATE_PROPAGATOR_DATA);
}


DeviceBuffer<float> StatePropagatorDataGpu::Impl::getForces()
{
    return d_f_;
}

void StatePropagatorDataGpu::Impl::copyForcesToGpu(const gmx::ArrayRef<const gmx::RVec> h_f,
                                                   AtomLocality atomLocality)
{
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    CommandStream commandStream = fCopyStreams_[atomLocality];
    GMX_ASSERT(commandStream != nullptr,
               "No stream is valid for copying forces with given atom locality.");

    wallcycle_start_nocount(wcycle_, ewcLAUNCH_GPU);
    wallcycle_sub_start(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);

    copyToDevice(d_f_, h_f, d_fSize_, atomLocality, commandStream);
    fReadyOnDevice_[atomLocality].markEvent(commandStream);

    wallcycle_sub_stop(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);
    wallcycle_stop(wcycle_, ewcLAUNCH_GPU);
}

GpuEventSynchronizer* StatePropagatorDataGpu::Impl::getForcesReadyOnDeviceEvent(AtomLocality atomLocality,
                                                                                bool useGpuFBufferOps)
{
    if ((atomLocality == AtomLocality::Local || atomLocality == AtomLocality::NonLocal) && useGpuFBufferOps)
    {
        return &fReducedOnDevice_;
    }
    else
    {
        return &fReadyOnDevice_[atomLocality];
    }
}

GpuEventSynchronizer* StatePropagatorDataGpu::Impl::fReducedOnDevice()
{
    return &fReducedOnDevice_;
}

void StatePropagatorDataGpu::Impl::copyForcesFromGpu(gmx::ArrayRef<gmx::RVec> h_f, AtomLocality atomLocality)
{
    GMX_ASSERT(atomLocality < AtomLocality::Count, "Wrong atom locality.");
    CommandStream commandStream = fCopyStreams_[atomLocality];
    GMX_ASSERT(commandStream != nullptr,
               "No stream is valid for copying forces with given atom locality.");

    wallcycle_start_nocount(wcycle_, ewcLAUNCH_GPU);
    wallcycle_sub_start(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);

    copyFromDevice(h_f, d_f_, d_fSize_, atomLocality, commandStream);
    fReadyOnHost_[atomLocality].markEvent(commandStream);

    wallcycle_sub_stop(wcycle_, ewcsLAUNCH_STATE_PROPAGATOR_DATA);
    wallcycle_stop(wcycle_, ewcLAUNCH_GPU);
}

void StatePropagatorDataGpu::Impl::waitForcesReadyOnHost(AtomLocality atomLocality)
{
    wallcycle_start(wcycle_, ewcWAIT_GPU_STATE_PROPAGATOR_DATA);
    fReadyOnHost_[atomLocality].waitForEvent();
    wallcycle_stop(wcycle_, ewcWAIT_GPU_STATE_PROPAGATOR_DATA);
}

void* StatePropagatorDataGpu::Impl::getUpdateStream()
{
    return &updateStream_;
}

int StatePropagatorDataGpu::Impl::numAtomsLocal()
{
    return numAtomsLocal_;
}

int StatePropagatorDataGpu::Impl::numAtomsAll()
{
    return numAtomsAll_;
}


StatePropagatorDataGpu::StatePropagatorDataGpu(const void*        pmeStream,
                                               const void*        localStream,
                                               const void*        nonLocalStream,
                                               const void*        deviceContext,
                                               GpuApiCallBehavior transferKind,
                                               int                paddingSize,
                                               gmx_wallcycle*     wcycle) :
    impl_(new Impl(pmeStream, localStream, nonLocalStream, deviceContext, transferKind, paddingSize, wcycle))
{
}

StatePropagatorDataGpu::StatePropagatorDataGpu(const void*        pmeStream,
                                               const void*        deviceContext,
                                               GpuApiCallBehavior transferKind,
                                               int                paddingSize,
                                               gmx_wallcycle*     wcycle) :
    impl_(new Impl(pmeStream, deviceContext, transferKind, paddingSize, wcycle))
{
}

StatePropagatorDataGpu::StatePropagatorDataGpu(StatePropagatorDataGpu&& /* other */) noexcept = default;

StatePropagatorDataGpu& StatePropagatorDataGpu::operator=(StatePropagatorDataGpu&& /* other */) noexcept = default;

StatePropagatorDataGpu::~StatePropagatorDataGpu() = default;


void StatePropagatorDataGpu::reinit(int numAtomsLocal, int numAtomsAll)
{
    return impl_->reinit(numAtomsLocal, numAtomsAll);
}

std::tuple<int, int> StatePropagatorDataGpu::getAtomRangesFromAtomLocality(AtomLocality atomLocality)
{
    return impl_->getAtomRangesFromAtomLocality(atomLocality);
}


DeviceBuffer<float> StatePropagatorDataGpu::getCoordinates()
{
    return impl_->getCoordinates();
}

void StatePropagatorDataGpu::copyCoordinatesToGpu(const gmx::ArrayRef<const gmx::RVec> h_x,
                                                  AtomLocality                         atomLocality)
{
    return impl_->copyCoordinatesToGpu(h_x, atomLocality);
}

GpuEventSynchronizer*
StatePropagatorDataGpu::getCoordinatesReadyOnDeviceEvent(AtomLocality              atomLocality,
                                                         const SimulationWorkload& simulationWork,
                                                         const StepWorkload&       stepWork)
{
    return impl_->getCoordinatesReadyOnDeviceEvent(atomLocality, simulationWork, stepWork);
}

void StatePropagatorDataGpu::waitCoordinatesCopiedToDevice(AtomLocality atomLocality)
{
    return impl_->waitCoordinatesCopiedToDevice(atomLocality);
}

GpuEventSynchronizer* StatePropagatorDataGpu::xUpdatedOnDevice()
{
    return impl_->xUpdatedOnDevice();
}

void StatePropagatorDataGpu::copyCoordinatesFromGpu(gmx::ArrayRef<RVec> h_x, AtomLocality atomLocality)
{
    return impl_->copyCoordinatesFromGpu(h_x, atomLocality);
}

void StatePropagatorDataGpu::waitCoordinatesReadyOnHost(AtomLocality atomLocality)
{
    return impl_->waitCoordinatesReadyOnHost(atomLocality);
}


DeviceBuffer<float> StatePropagatorDataGpu::getVelocities()
{
    return impl_->getVelocities();
}

void StatePropagatorDataGpu::copyVelocitiesToGpu(const gmx::ArrayRef<const gmx::RVec> h_v,
                                                 AtomLocality                         atomLocality)
{
    return impl_->copyVelocitiesToGpu(h_v, atomLocality);
}

GpuEventSynchronizer* StatePropagatorDataGpu::getVelocitiesReadyOnDeviceEvent(AtomLocality atomLocality)
{
    return impl_->getVelocitiesReadyOnDeviceEvent(atomLocality);
}

void StatePropagatorDataGpu::copyVelocitiesFromGpu(gmx::ArrayRef<RVec> h_v, AtomLocality atomLocality)
{
    return impl_->copyVelocitiesFromGpu(h_v, atomLocality);
}

void StatePropagatorDataGpu::waitVelocitiesReadyOnHost(AtomLocality atomLocality)
{
    return impl_->waitVelocitiesReadyOnHost(atomLocality);
}


DeviceBuffer<float> StatePropagatorDataGpu::getForces()
{
    return impl_->getForces();
}

void StatePropagatorDataGpu::copyForcesToGpu(const gmx::ArrayRef<const gmx::RVec> h_f, AtomLocality atomLocality)
{
    return impl_->copyForcesToGpu(h_f, atomLocality);
}

GpuEventSynchronizer* StatePropagatorDataGpu::getForcesReadyOnDeviceEvent(AtomLocality atomLocality,
                                                                          bool useGpuFBufferOps)
{
    return impl_->getForcesReadyOnDeviceEvent(atomLocality, useGpuFBufferOps);
}

GpuEventSynchronizer* StatePropagatorDataGpu::fReducedOnDevice()
{
    return impl_->fReducedOnDevice();
}

void StatePropagatorDataGpu::copyForcesFromGpu(gmx::ArrayRef<RVec> h_f, AtomLocality atomLocality)
{
    return impl_->copyForcesFromGpu(h_f, atomLocality);
}

void StatePropagatorDataGpu::waitForcesReadyOnHost(AtomLocality atomLocality)
{
    return impl_->waitForcesReadyOnHost(atomLocality);
}


void* StatePropagatorDataGpu::getUpdateStream()
{
    return impl_->getUpdateStream();
}

int StatePropagatorDataGpu::numAtomsLocal()
{
    return impl_->numAtomsLocal();
}

int StatePropagatorDataGpu::numAtomsAll()
{
    return impl_->numAtomsAll();
}

} // namespace gmx

#endif // GMX_GPU == GMX_GPU_NONE
