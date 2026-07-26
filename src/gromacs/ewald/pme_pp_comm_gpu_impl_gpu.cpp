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
 * \brief Implements PME-PP communication using GPU direct communication
 *
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "config.h"

#include <cstdlib>

#include <array>
#include <atomic>

#include "gromacs/ewald/pme_pp_communication.h"
#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"

#include "pme_pp_comm_gpu_impl.h"
#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/cudautils.cuh"
#    include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#endif
#if GMX_GPU_HIP
#    include "gromacs/gpu_utils/hiputils.h"
#    include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#endif
#if GMX_GPU_SYCL
#    include "gromacs/gpu_utils/gmxsycl.h"
#endif
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{

PmePpCommGpu::Impl::Impl(MPI_Comm             comm,
                         int                  pmeRank,
                         const DeviceContext& deviceContext,
                         const DeviceStream&  deviceStream,
                         bool                 useNvshmem) :
    deviceContext_(deviceContext),
    pmePpCommStream_(deviceStream),
    comm_(comm),
    pmeRank_(pmeRank),
    d_pmeForces_(nullptr),
    forcesReadyNvshmemFlags(nullptr),
    useNvshmem_(useNvshmem)
{
    stageLibMpiGpuCpuComm_ = (std::getenv("GMX_DISABLE_STAGED_GPU_TO_CPU_PMEPP_COMM") == nullptr);
}

PmePpCommGpu::Impl::~Impl()
{
    try
    {
        freeDeviceBuffer(&d_pmeForces_);
    }
    catch (gmx::InternalError& e)
    {
        fprintf(stderr, "Internal error in destructor of PmePpCommGpu: %s\n", e.what());
    }
}

void PmePpCommGpu::Impl::reinit(ArrayRef<RVec> pmeCpuForceReceiveBuffer)
{
    pmeCpuForceReceiveBuffer_ = pmeCpuForceReceiveBuffer;
    int newSize               = pmeCpuForceReceiveBuffer_.size();
    if (useNvshmem_)
    {
        /* Find the max nAtomsAlloc among all ranks for symmetric force-buffer allocation
         * and the max number of PP peers per PME rank for symmetric signal allocation.
         * PME ranks participate in the matching collective from pme_gpu_reinit_atoms. */
        std::array<int, 2> buf{ newSize, 0 };
#if GMX_MPI
        MPI_Allreduce(MPI_IN_PLACE, buf.data(), buf.size(), MPI_INT, MPI_MAX, comm_);
#endif
        newSize              = buf[0];
        const int numPpRanks = buf[1];
        // symmetric buffer used for synchronization purpose 1 to be used to signal PME to PP rank
        // of put, and numPpRanks is intended to be used for each PP rank buffer consumption
        // completion signal to PME to allow to produce it again. this a collective call.
        reallocateDeviceBuffer(&forcesReadyNvshmemFlags,
                               1 + numPpRanks,
                               &forcesReadyNvshmemFlagsSize_,
                               &forcesReadyNvshmemFlagsSizeAlloc_,
                               deviceContext_,
                               true);
    }

    // Reallocate device buffer used for staging PME force.
    // if useNvshmem_ is true, this a collective call, resulting in symmetric allocation involving PME + PP ranks.
    reallocateDeviceBuffer(
            &d_pmeForces_, newSize, &d_pmeForcesSize_, &d_pmeForcesSizeAlloc_, deviceContext_, useNvshmem_);

    // This rank will access PME rank memory directly, so needs to receive the remote PME buffer addresses.
#if GMX_MPI
    if (GMX_THREAD_MPI)
    {
        // receive device coordinate buffer address from PME rank
        MPI_Recv(&remotePmeXBuffer_,
                 sizeof(Float3*),
                 MPI_BYTE,
                 pmeRank_,
                 eCommType_COORD_GPU_REMOTE_GPU_PTR,
                 comm_,
                 MPI_STATUS_IGNORE);
        // send host and device force buffer addresses to PME rank
        Float3* pmeForcePtr = asRawDevicePointer(d_pmeForces_);
        MPI_Send(&pmeForcePtr, sizeof(Float3*), MPI_BYTE, pmeRank_, eCommType_FORCES_GPU_REMOTE_GPU_PTR, comm_);
        RVec* pmeCpuForceReceiveBufferData = pmeCpuForceReceiveBuffer_.data();
        MPI_Send(&pmeCpuForceReceiveBufferData,
                 sizeof(RVec*),
                 MPI_BYTE,
                 pmeRank_,
                 eCommType_FORCES_GPU_REMOTE_CPU_PTR,
                 comm_);
        // Receive address of event and associated flag from PME rank, to allow sync to local stream after force transfer
        // NOLINTNEXTLINE(bugprone-sizeof-expression)
        MPI_Recv(&remotePmeForceSendEvent_,
                 sizeof(GpuEventSynchronizer*),
                 MPI_BYTE,
                 pmeRank_,
                 eCommType_FORCES_GPU_SYNCHRONIZER,
                 comm_,
                 MPI_STATUS_IGNORE);
        MPI_Recv(&remotePmeForceSendEventRecorded_,
                 sizeof(std::atomic<bool>*),
                 MPI_BYTE,
                 pmeRank_,
                 eCommType_FORCES_GPU_EVENT_RECORDED,
                 comm_,
                 MPI_STATUS_IGNORE);
    }
#endif
}

void PmePpCommGpu::Impl::receiveForceFromPmeDirectGpu(bool receivePmeForceToGpu)
{
    // Wait until remote PME task has pushed data, and then enqueue remote event to local stream.

    // Spin until PME rank sets flag
    while (!(remotePmeForceSendEventRecorded_->load(std::memory_order_acquire))) {}

    // Enqueue remote event
    remotePmeForceSendEvent_->enqueueWaitEvent(pmePpCommStream_);

    // Reset the flag
    remotePmeForceSendEventRecorded_->store(false, std::memory_order_release);

    if (receivePmeForceToGpu)
    {
        // Peer-to-peer copy
        //
        // Record event to be enqueued in the GPU local buffer operations, to
        // satisfy dependency on receiving the PME force data before
        // reducing it with the other force contributions.
        forcesReadySynchronizer_.markEvent(pmePpCommStream_);
    }
    else
    {
        // Ensure CPU waits for PME forces to be copied before reducing
        // them with other forces on the CPU
        pmePpCommStream_.synchronize();
    }
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void PmePpCommGpu::Impl::receiveForceFromPmeGpuAwareMpi(const bool receivePmeForceToGpu)
{
#if GMX_LIB_MPI
    GMX_RELEASE_ASSERT(coordinateSendRequestIsActive_,
                       "Coordinates must have been sent for a domain (even if empty)");
    // Wait on previous non-blocking coordinate send. This already
    // must have completed for PME forces to be ready, but the wait is
    // required by the MPI standard to assure completion.
    MPI_Wait(&coordinateSendRequest_, MPI_STATUS_IGNORE);
    coordinateSendRequestIsActive_ = false;

    // The PME rank always sends forces, even when the domain is
    // empty, so each PP rank must always post a receive.
    Float3*   pmeForcePtr = receivePmeForceToGpu ? asRawDevicePointer(d_pmeForces_)
                                                 : pmeCpuForceReceiveBuffer_.data();
    const int recvSize = receivePmeForceToGpu ? d_pmeForcesSize_ : pmeCpuForceReceiveBuffer_.size();
    // Ensure that the two PME force receive buffers were resized
    // consistently, when relevant.
    GMX_ASSERT(receivePmeForceToGpu || (d_pmeForcesSize_ == gmx::ssize(pmeCpuForceReceiveBuffer_))
                       || (useNvshmem_ && d_pmeForcesSize_ >= gmx::ssize(pmeCpuForceReceiveBuffer_)),
               "Mismatch between sizes of GPU and CPU PME force receive buffers");
    if (forceRecvRequestIsActive_)
    {
        MPI_Wait(&forceRecvRequest_, MPI_STATUS_IGNORE);
        forceRecvRequestIsActive_ = false;
    }
    else if (!stageLibMpiGpuCpuComm_ && !(useNvshmem_ && receivePmeForceToGpu))
    {
        MPI_Recv(pmeForcePtr, recvSize * DIM, MPI_FLOAT, pmeRank_, eCommType_FORCES_GPU, comm_, MPI_STATUS_IGNORE);
    }

    if (stageLibMpiGpuCpuComm_
        && !receivePmeForceToGpu) // destination is CPU memory, so finalize transfer with local D2H
    {
        copyFromDeviceBuffer(reinterpret_cast<RVec*>(pmeForcePtr),
                             &d_pmeForces_,
                             0,
                             recvSize,
                             pmePpCommStream_,
                             GpuApiCallBehavior::Sync,
                             nullptr);
    }

#else
    GMX_UNUSED_VALUE(receivePmeForceToGpu);
#endif
}

void PmePpCommGpu::Impl::receiveForceFromPme(const bool receivePmeForceToGpu)
{
    if (GMX_THREAD_MPI)
    {
        receiveForceFromPmeDirectGpu(receivePmeForceToGpu);
    }
    else
    {
        receiveForceFromPmeGpuAwareMpi(receivePmeForceToGpu);
    }
}

void PmePpCommGpu::Impl::sendCoordinatesToPmeDirectGpu(const Float3* sendPtr,
                                                       const int     sendSize,
                                                       const bool    sendCoordinatesFromGpu,
                                                       GpuEventSynchronizer* coordinatesReadyOnDeviceEvent)
{
    GMX_ASSERT(GpuConfigurationCapabilities::ThreadMpiDirectComm,
               "Direct GPU communications only supported with CUDA and HIP.");
    // ensure stream waits until coordinate data is available on device
    if (coordinatesReadyOnDeviceEvent)
    {
        coordinatesReadyOnDeviceEvent->enqueueWaitEvent(pmePpCommStream_);
    }

    // Destination is always a remote device
    if (sendCoordinatesFromGpu)
    {
        // Peer-to-peer copy
        copyBetweenDeviceBuffers(
                remotePmeXBuffer_, sendPtr, 0, sendSize, pmePpCommStream_, GpuApiCallBehavior::Async, nullptr);
    }
    else
    {
        copyToDeviceBuffer(
                &remotePmeXBuffer_, sendPtr, 0, sendSize, pmePpCommStream_, GpuApiCallBehavior::Async, nullptr);
    }

#if GMX_MPI
    // Record and send event to allow PME task to sync to above transfer before commencing force calculations
    pmeCoordinatesSynchronizer_.markEvent(pmePpCommStream_);
    GpuEventSynchronizer* pmeSync = &pmeCoordinatesSynchronizer_;
    // NOLINTNEXTLINE(bugprone-sizeof-expression)
    MPI_Send(&pmeSync, sizeof(GpuEventSynchronizer*), MPI_BYTE, pmeRank_, eCommType_COORD_GPU_SYNCHRONIZER, comm_);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void PmePpCommGpu::Impl::sendCoordinatesToPmeGpuAwareMpi(const Float3* sendPtr,
                                                         int           sendSize,
                                                         GpuEventSynchronizer* coordinatesReadyOnDeviceEvent,
                                                         bool receivePmeForceToGpu)
{
    // Post the non-blocking send first, as it is definitely on the
    // critical path.

    // ensure coordinate data is available on device before we start transfer
    if (coordinatesReadyOnDeviceEvent)
    {
        coordinatesReadyOnDeviceEvent->waitForEvent();
    }

#if GMX_LIB_MPI
    // The corresponding wait for the below non-blocking coordinate send is in receiveForceFromPmeGpuAwareMpi.
    MPI_Isend(sendPtr, sendSize * DIM, MPI_FLOAT, pmeRank_, eCommType_COORD_GPU, comm_, &coordinateSendRequest_);
    coordinateSendRequestIsActive_ = true;
#else
    GMX_UNUSED_VALUE(sendPtr);
    GMX_UNUSED_VALUE(sendSize);
#endif

    // The PME rank always sends forces, even when the domain is
    // empty, so each PP rank must always post a receive. This is
    // posted after the possible coordinate-send, so the latter is not
    // delayed.
    if (stageLibMpiGpuCpuComm_ && !(useNvshmem_ && receivePmeForceToGpu))
    {
#if GMX_LIB_MPI
        // A non-blocking receive is used so that the data transfer
        // neither blocks nor is blocked by other PP MPI activities.
        MPI_Irecv(asRawDevicePointer(d_pmeForces_),
                  sendSize * DIM,
                  MPI_FLOAT,
                  pmeRank_,
                  eCommType_FORCES_GPU,
                  comm_,
                  &forceRecvRequest_);
        forceRecvRequestIsActive_ = true;
#endif
    }
}

std::optional<DeviceBuffer<Float3>> PmePpCommGpu::Impl::getGpuForceStagingPtr()
{
    if (d_pmeForcesSize_ > 0)
    {
        return d_pmeForces_;
    }
    else
    {
        return std::nullopt;
    }
}

std::optional<GpuEventSynchronizer*> PmePpCommGpu::Impl::getForcesReadySynchronizer()
{
    if (GMX_THREAD_MPI && d_pmeForcesSize_ > 0)
    {
        return &forcesReadySynchronizer_;
    }
    else
    {
        return std::nullopt;
    }
}

DeviceBuffer<uint64_t> PmePpCommGpu::Impl::getGpuForcesSyncObj()
{
    return forcesReadyNvshmemFlags;
}

PmePpCommGpu::PmePpCommGpu(MPI_Comm             comm,
                           int                  pmeRank,
                           const DeviceContext& deviceContext,
                           const DeviceStream&  deviceStream,
                           bool                 useNvshmem) :
    impl_(new Impl(comm, pmeRank, deviceContext, deviceStream, useNvshmem))
{
}

PmePpCommGpu::~PmePpCommGpu() = default;

void PmePpCommGpu::reinit(ArrayRef<RVec> pmeCpuForceReceiveBuffer)
{
    impl_->reinit(pmeCpuForceReceiveBuffer);
}

void PmePpCommGpu::receiveForceFromPme(const bool receivePmeForceToGpu)
{
    impl_->receiveForceFromPme(receivePmeForceToGpu);
}

void PmePpCommGpu::sendCoordinatesToPmeFromGpu(DeviceBuffer<RVec>    sendPtr,
                                               int                   sendSize,
                                               GpuEventSynchronizer* coordinatesReadyOnDeviceEvent,
                                               bool                  receiveForcesToGpu)
{
    if (GMX_THREAD_MPI)
    {
        const bool sendCoordinatesFromGpu = true;
        impl_->sendCoordinatesToPmeDirectGpu(
                asRawDevicePointer(sendPtr), sendSize, sendCoordinatesFromGpu, coordinatesReadyOnDeviceEvent);
    }
    else
    {
        impl_->sendCoordinatesToPmeGpuAwareMpi(
                asRawDevicePointer(sendPtr), sendSize, coordinatesReadyOnDeviceEvent, receiveForcesToGpu);
    }
}

void PmePpCommGpu::sendCoordinatesToPmeFromCpu(const RVec* sendPtr, int sendSize, bool receiveForcesToGpu)
{
    GpuEventSynchronizer* coordinatesReadyOnDeviceEvent = nullptr;
    if (GMX_THREAD_MPI)
    {
        const bool sendCoordinatesFromGpu = false;
        impl_->sendCoordinatesToPmeDirectGpu(
                sendPtr, sendSize, sendCoordinatesFromGpu, coordinatesReadyOnDeviceEvent);
    }
    else
    {
        impl_->sendCoordinatesToPmeGpuAwareMpi(
                sendPtr, sendSize, coordinatesReadyOnDeviceEvent, receiveForcesToGpu);
    }
}

std::optional<DeviceBuffer<Float3>> PmePpCommGpu::getGpuForceStagingPtr()
{
    return impl_->getGpuForceStagingPtr();
}

std::optional<GpuEventSynchronizer*> PmePpCommGpu::getForcesReadySynchronizer()
{
    return impl_->getForcesReadySynchronizer();
}

DeviceBuffer<uint64_t> PmePpCommGpu::getGpuForcesSyncObj()
{
    return impl_->getGpuForcesSyncObj();
}

} // namespace gmx
