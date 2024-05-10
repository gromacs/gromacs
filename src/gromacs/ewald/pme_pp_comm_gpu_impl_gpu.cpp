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

#include "gromacs/ewald/pme_pp_communication.h"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"

#include "pme_pp_comm_gpu_impl.h"
#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/cudautils.cuh"
#    include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#endif
#if GMX_GPU_SYCL
#    include "gromacs/gpu_utils/gmxsycl.h"
#endif
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{

PmePpCommGpu::Impl::Impl(MPI_Comm                    comm,
                         int                         pmeRank,
                         gmx::HostVector<gmx::RVec>* pmeCpuForceBuffer,
                         const DeviceContext&        deviceContext,
                         const DeviceStream&         deviceStream,
                         bool                        useNvshmem) :
    deviceContext_(deviceContext),
    pmePpCommStream_(deviceStream),
    comm_(comm),
    pmeRank_(pmeRank),
    pmeCpuForceBuffer_(pmeCpuForceBuffer),
    d_pmeForces_(nullptr),
    forcesReadyNvshmemFlags(nullptr),
    useNvshmem_(useNvshmem)
{
    stageLibMpiGpuCpuComm_ = (getenv("GMX_DISABLE_STAGED_GPU_TO_CPU_PMEPP_COMM") == nullptr);
}

PmePpCommGpu::Impl::~Impl()
{
    freeDeviceBuffer(&d_pmeForces_);
}

void PmePpCommGpu::Impl::reinit(int size)
{
    int newSize = size;
    if (useNvshmem_)
    {
        MPI_Allreduce(&size, &newSize, 1, MPI_INT, MPI_MAX, comm_);

        int numPpRanks = 0;
        MPI_Bcast(&numPpRanks, 1, MPI_INT, pmeRank_, comm_);
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
#if GMX_THREAD_MPI
    if (GMX_THREAD_MPI)
    {
        // receive device coordinate buffer address from PME rank
        MPI_Recv(&remotePmeXBuffer_, sizeof(Float3*), MPI_BYTE, pmeRank_, 0, comm_, MPI_STATUS_IGNORE);
        // send host and device force buffer addresses to PME rank
        MPI_Send(&d_pmeForces_, sizeof(Float3*), MPI_BYTE, pmeRank_, 0, comm_);
        RVec* pmeCpuForceBufferData = pmeCpuForceBuffer_->data();
        MPI_Send(&pmeCpuForceBufferData, sizeof(RVec*), MPI_BYTE, pmeRank_, 0, comm_);
        // Receive address of event and associated flag from PME rank, to allow sync to local stream after force transfer
        // NOLINTNEXTLINE(bugprone-sizeof-expression)
        MPI_Recv(&remotePmeForceSendEvent_, sizeof(GpuEventSynchronizer*), MPI_BYTE, pmeRank_, 0, comm_, MPI_STATUS_IGNORE);
        MPI_Recv(&remotePmeForceSendEventRecorded_, sizeof(std::atomic<bool>*), MPI_BYTE, pmeRank_, 0, comm_, MPI_STATUS_IGNORE);
    }

#endif
}

void PmePpCommGpu::Impl::receiveForceFromPmePeerToPeer(bool receivePmeForceToGpu)
{
#if GMX_MPI
    // Wait until remote PME task has pushed data, and then enqueue remote event to local stream.

    if (d_pmeForcesSize_ <= 0)
    {
        return;
    }

    // Spin until PME rank sets flag
    while (!(remotePmeForceSendEventRecorded_->load(std::memory_order_acquire))) {}

    // Enqueue remote event
    remotePmeForceSendEvent_->enqueueWaitEvent(pmePpCommStream_);

    // Reset the flag
    remotePmeForceSendEventRecorded_->store(false, std::memory_order_release);

    if (receivePmeForceToGpu)
    {
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
#endif
}

// NOLINTNEXTLINE readability-convert-member-functions-to-static
void PmePpCommGpu::Impl::receiveForceFromPmeGpuAwareMpi(Float3* pmeForcePtr, int recvSize)
{
#if GMX_LIB_MPI

    // Wait on previous non-blocking coordinate send. This already must have completed for PME
    // forces to be ready, but the wait is necessary to avoid issues with certain MPI libraries.
    GMX_ASSERT(coordinateSendRequestIsActive_,
               "A coordinate send request should be active before force is recieved");
    MPI_Wait(&coordinateSendRequest_, MPI_STATUS_IGNORE);
    coordinateSendRequestIsActive_ = false;

    if (!stageLibMpiGpuCpuComm_)
    {
        MPI_Recv(pmeForcePtr, recvSize * DIM, MPI_FLOAT, pmeRank_, 0, comm_, MPI_STATUS_IGNORE);
    }
    else
    {
        if (useNvshmem_)
        {
            // destination is CPU memory, so finalize transfer with local D2H
            if (pmeForcePtr != asMpiPointer(d_pmeForces_))
            {
                // Receive data from remote GPU in memory of local GPU
                MPI_Recv(asMpiPointer(d_pmeForces_), recvSize * DIM, MPI_FLOAT, pmeRank_, 0, comm_, MPI_STATUS_IGNORE);
            }
        }
        else
        {
            // Receive data from remote GPU in memory of local GPU
            MPI_Recv(asMpiPointer(d_pmeForces_), recvSize * DIM, MPI_FLOAT, pmeRank_, 0, comm_, MPI_STATUS_IGNORE);
        }

        if (pmeForcePtr != asMpiPointer(d_pmeForces_)) // destination is CPU memory, so finalize transfer with local D2H
        {
            copyFromDeviceBuffer(reinterpret_cast<RVec*>(pmeForcePtr),
                                 &d_pmeForces_,
                                 0,
                                 recvSize,
                                 pmePpCommStream_,
                                 GpuApiCallBehavior::Sync,
                                 nullptr);
        }
    }
#else
    GMX_UNUSED_VALUE(pmeForcePtr);
    GMX_UNUSED_VALUE(recvSize);
#endif
}

void PmePpCommGpu::Impl::receiveForceFromPme(Float3* recvPtr, int recvSize, bool receivePmeForceToGpu)
{
    Float3* pmeForcePtr = receivePmeForceToGpu ? asMpiPointer(d_pmeForces_) : recvPtr;
    if (GMX_THREAD_MPI)
    {
        receiveForceFromPmePeerToPeer(receivePmeForceToGpu);
    }
    else
    {
        receiveForceFromPmeGpuAwareMpi(pmeForcePtr, recvSize);
    }
}

// NOLINTNEXTLINE readability-convert-member-functions-to-static
void PmePpCommGpu::Impl::sendCoordinatesToPmeGpuAwareMpi(Float3*               sendPtr,
                                                         int                   sendSize,
                                                         GpuEventSynchronizer* coordinatesReadyOnDeviceEvent)
{
    // ensure coordinate data is available on device before we start transfer
    if (coordinatesReadyOnDeviceEvent)
    {
        coordinatesReadyOnDeviceEvent->waitForEvent();
    }

#if GMX_LIB_MPI
    Float3* sendptr_x = sendPtr;
    // The corresponding wait for the below non-blocking coordinate send is in receiveForceFromPmeGpuAwareMpi.
    // Strictly, a wait is not necessary since the recieve must complete before PME forces are calculated,
    // but it is required to avoid issues in certain MPI libraries.
    MPI_Isend(sendptr_x, sendSize * DIM, MPI_FLOAT, pmeRank_, eCommType_COORD_GPU, comm_, &coordinateSendRequest_);
    coordinateSendRequestIsActive_ = true;
#else
    GMX_UNUSED_VALUE(sendPtr);
    GMX_UNUSED_VALUE(sendSize);
#endif
}

void PmePpCommGpu::Impl::sendCoordinatesToPme(Float3*               sendPtr,
                                              int                   sendSize,
                                              GpuEventSynchronizer* coordinatesReadyOnDeviceEvent)
{
    if (GMX_THREAD_MPI)
    {
        sendCoordinatesToPmePeerToPeer(sendPtr, sendSize, coordinatesReadyOnDeviceEvent);
    }
    else
    {
        sendCoordinatesToPmeGpuAwareMpi(sendPtr, sendSize, coordinatesReadyOnDeviceEvent);
    }
}
DeviceBuffer<Float3> PmePpCommGpu::Impl::getGpuForceStagingPtr()
{
    return d_pmeForces_;
}

GpuEventSynchronizer* PmePpCommGpu::Impl::getForcesReadySynchronizer()
{
    if (GMX_THREAD_MPI)
    {
        return &forcesReadySynchronizer_;
    }
    else
    {
        return nullptr;
    }
}

DeviceBuffer<uint64_t> PmePpCommGpu::Impl::getGpuForcesSyncObj()
{
    return forcesReadyNvshmemFlags;
}

PmePpCommGpu::PmePpCommGpu(MPI_Comm                    comm,
                           int                         pmeRank,
                           gmx::HostVector<gmx::RVec>* pmeCpuForceBuffer,
                           const DeviceContext&        deviceContext,
                           const DeviceStream&         deviceStream,
                           bool                        useNvshmem) :
    impl_(new Impl(comm, pmeRank, pmeCpuForceBuffer, deviceContext, deviceStream, useNvshmem))
{
}

PmePpCommGpu::~PmePpCommGpu() = default;

void PmePpCommGpu::reinit(int size)
{
    impl_->reinit(size);
}

void PmePpCommGpu::receiveForceFromPme(RVec* recvPtr, int recvSize, bool receivePmeForceToGpu)
{
    impl_->receiveForceFromPme(recvPtr, recvSize, receivePmeForceToGpu);
}

void PmePpCommGpu::sendCoordinatesToPmeFromGpu(DeviceBuffer<RVec>    sendPtr,
                                               int                   sendSize,
                                               GpuEventSynchronizer* coordinatesReadyOnDeviceEvent)
{
    impl_->sendCoordinatesToPme(asMpiPointer(sendPtr), sendSize, coordinatesReadyOnDeviceEvent);
}

void PmePpCommGpu::sendCoordinatesToPmeFromCpu(RVec* sendPtr, int sendSize)
{
    impl_->sendCoordinatesToPme(sendPtr, sendSize, nullptr);
}

DeviceBuffer<Float3> PmePpCommGpu::getGpuForceStagingPtr()
{
    return impl_->getGpuForceStagingPtr();
}

GpuEventSynchronizer* PmePpCommGpu::getForcesReadySynchronizer()
{
    return impl_->getForcesReadySynchronizer();
}

DeviceBuffer<uint64_t> PmePpCommGpu::getGpuForcesSyncObj()
{
    return impl_->getGpuForcesSyncObj();
}

} // namespace gmx
