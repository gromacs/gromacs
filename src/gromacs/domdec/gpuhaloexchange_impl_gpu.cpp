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
 * \brief Implements shared code for GPU halo exchange.
 *
 *
 * \author Alan Gray <alang@nvidia.com>
 * \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \ingroup module_domdec
 */
#include "gmxpre.h"

#include "gpuhaloexchange_impl_gpu.h"

#include "config.h"

#include "gromacs/gpu_utils/gpu_utils.h"

#include "fused_gpuhaloexchange.h"

#if GMX_NVSHMEM
#    include <nvshmem.h>
#    include <nvshmemx.h>
#endif

#include <cstdint>

#include <numeric>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/vectypes.h"

#include "domdec_internal.h"
#include "halo_plan_utils.h"

// NOLINTNEXTLINE(misc-redundant-expression)
constexpr bool supportedLibMpiBuild =
        ((GMX_LIB_MPI != 0) && gmx::GpuConfigurationCapabilities::HaloExchangeDirectComm);
// NOLINTNEXTLINE(misc-redundant-expression)
constexpr bool supportedThreadMpiBuild =
        ((GMX_THREAD_MPI != 0) && gmx::GpuConfigurationCapabilities::HaloExchangeDirectComm);

namespace gmx
{

HaloPlan computeHaloPlan(const gmx_domdec_comm_t& comm,
                         int                      dimIndex,
                         int                      pulse,
                         MPI_Comm                 mpiCommMySim,
                         int                      sendRankX,
                         int                      recvRankX)
{
    HaloPlan plan;

    const gmx_domdec_comm_dim_t& cd = comm.cd[dimIndex];
    plan.ind                        = &cd.ind[pulse];
    plan.receiveInPlace             = cd.receiveInPlace;

    int numZoneTemp   = 1;
    int numAtomsTotal = comm.atomRanges.numHomeAtoms();
    for (int i = 0; i <= dimIndex; i++)
    {
        GMX_ASSERT(numZoneTemp <= sc_maxNumIZones, "Too many domain decomposition zones");
        int pulseMax = (i == dimIndex) ? pulse : (comm.cd[i].numPulses() - 1);
        for (int p = 0; p <= pulseMax; p++)
        {
            plan.atomOffset                 = numAtomsTotal;
            const gmx_domdec_ind_t& indTemp = comm.cd[i].ind[p];
            numAtomsTotal += indTemp.nrecv[numZoneTemp + 1];
        }
        plan.numZone = numZoneTemp;
        numZoneTemp += numZoneTemp;
    }

    plan.xSendSize = plan.ind->nsend[plan.numZone + 1];

#if GMX_MPI
    // exchange send/recv sizes for X
    int xRecv = 0;
    MPI_Sendrecv(&plan.xSendSize,
                 sizeof(int),
                 MPI_BYTE,
                 sendRankX,
                 0,
                 &xRecv,
                 sizeof(int),
                 MPI_BYTE,
                 recvRankX,
                 0,
                 mpiCommMySim,
                 MPI_STATUS_IGNORE);
    plan.xRecvSize = xRecv;
#else
    plan.xRecvSize = plan.xSendSize;
#endif

    return plan;
}

void GpuHaloExchange::Impl::reinitHalo(DeviceBuffer<Float3> d_coordinatesBuffer,
                                       DeviceBuffer<Float3> d_forcesBuffer)
{
    GMX_RELEASE_ASSERT(supportedLibMpiBuild || supportedThreadMpiBuild,
                       "Gpu Halo Exchange not supported in this build");

    wallcycle_start_nocount(wcycle_, WallCycleCounter::Domdec);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::DDGpu);

    d_x_ = d_coordinatesBuffer;
    d_f_ = d_forcesBuffer;

    const gmx_domdec_comm_t& comm = *dd_->comm;

    // Common plan for this dim/pulse
    const auto plan = computeHaloPlan(comm, dimIndex_, pulse_, mpi_comm_mysim_, sendRankX_, recvRankX_);
    ind_            = plan.ind;
    receiveInPlace_ = plan.receiveInPlace;
    numHomeAtoms_   = comm.atomRanges.numHomeAtoms();
    atomOffset_     = plan.atomOffset;
    numZone_        = plan.numZone;

    int newSize = plan.xSendSize;

    // reallocates only if needed
    h_indexMap_.resize(newSize);
    // reallocate on device only if needed
    if (newSize > maxPackedBufferSize_)
    {
        reallocateDeviceBuffer(&d_indexMap_, newSize, &indexMapSize_, &indexMapSizeAlloc_, deviceContext_);
        reallocateDeviceBuffer(&d_sendBuf_, newSize, &sendBufSize_, &sendBufSizeAlloc_, deviceContext_);
        // number of values/elems is same for indexMapSizeAlloc_/sendBufSizeAlloc_ so we can use either.
        maxPackedBufferSize_ = sendBufSizeAlloc_;
    }

    xSendSize_ = newSize;
    reallocateDeviceBuffer(&d_recvBuf_, newSize, &recvBufSize_, &recvBufSizeAlloc_, deviceContext_);

    xRecvSize_ = plan.xRecvSize;
    fSendSize_ = xRecvSize_;
    fRecvSize_ = xSendSize_;

    if (!receiveInPlace_)
    {
        // Same buffers will be used for both coordinates and forces
        changePinningPolicy(&h_outOfPlaceRecvBuffer_, PinningPolicy::PinnedIfSupported);
        changePinningPolicy(&h_outOfPlaceSendBuffer_, PinningPolicy::PinnedIfSupported);
        h_outOfPlaceSendBuffer_.resize(std::max(xSendSize_, fSendSize_));
        h_outOfPlaceRecvBuffer_.resize(std::max(xRecvSize_, fRecvSize_));
    }

    if (newSize > 0)
    {
        GMX_ASSERT(ind_->index.size() == h_indexMap_.size(),
                   "Size mismatch between domain decomposition communication index array and GPU "
                   "halo exchange index mapping array");
        std::copy(ind_->index.begin(), ind_->index.end(), h_indexMap_.begin());

        copyToDeviceBuffer(
                &d_indexMap_, h_indexMap_.data(), 0, newSize, *haloStream_, GpuApiCallBehavior::Async, nullptr);
    }

#if GMX_MPI && GMX_THREAD_MPI
    // Exchange of remote addresses from neighboring ranks is needed only with peer-to-peer copies as cudaMemcpy needs both src/dst pointer
    // MPI calls such as MPI_send doesn't worry about receiving address, that is taken care by MPI_recv call in neighboring rank
    if constexpr (supportedThreadMpiBuild)
    {
        // This rank will push data to its neighbor, so needs to know the remote receive address
        // and similarly send its receive address to other neighbour. We can do this here since
        // the pointers will not change until the next NS step.

        // Coordinates buffer:
        Float3* recvPtr = &asMpiPointer(d_x_)[atomOffset_];
        MPI_Sendrecv(&recvPtr,
                     sizeof(void*),
                     MPI_BYTE,
                     recvRankX_,
                     0,
                     &remoteXPtr_,
                     sizeof(void*),
                     MPI_BYTE,
                     sendRankX_,
                     0,
                     mpi_comm_mysim_,
                     MPI_STATUS_IGNORE);

        // Force buffer:
        recvPtr = asMpiPointer(d_recvBuf_);
        MPI_Sendrecv(&recvPtr,
                     sizeof(void*),
                     MPI_BYTE,
                     recvRankF_,
                     0,
                     &remoteFPtr_,
                     sizeof(void*),
                     MPI_BYTE,
                     sendRankF_,
                     0,
                     mpi_comm_mysim_,
                     MPI_STATUS_IGNORE);
    }
#endif

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::DDGpu);
    wallcycle_stop(wcycle_, WallCycleCounter::Domdec);
}

void GpuHaloExchange::Impl::enqueueWaitRemoteCoordinatesReadyEvent(GpuEventSynchronizer* coordinatesReadyOnDeviceEvent)
{
#if GMX_MPI
    GMX_ASSERT(coordinatesReadyOnDeviceEvent != nullptr,
               "Co-ordinate Halo exchange requires valid co-ordinate ready event");

    // Wait for event from receiving task that remote coordinates are ready, and enqueue that event to stream used
    // for subsequent data push. This avoids a race condition with the remote data being written in the previous timestep.
    // Similarly send event to task that will push data to this task.
    GpuEventSynchronizer* remoteCoordinatesReadyOnDeviceEvent;
    MPI_Sendrecv(&coordinatesReadyOnDeviceEvent,
                 sizeof(GpuEventSynchronizer*), //NOLINT(bugprone-sizeof-expression)
                 MPI_BYTE,
                 recvRankX_,
                 0,
                 &remoteCoordinatesReadyOnDeviceEvent,
                 sizeof(GpuEventSynchronizer*), //NOLINT(bugprone-sizeof-expression)
                 MPI_BYTE,
                 sendRankX_,
                 0,
                 mpi_comm_mysim_,
                 MPI_STATUS_IGNORE);
    remoteCoordinatesReadyOnDeviceEvent->enqueueWaitEvent(*haloStream_);
#else
    GMX_UNUSED_VALUE(coordinatesReadyOnDeviceEvent);
#endif
}

GpuEventSynchronizer* GpuHaloExchange::Impl::communicateHaloCoordinates(const matrix box,
                                                                        GpuEventSynchronizer* dependencyEvent)
{
    wallcycle_start(wcycle_, WallCycleCounter::LaunchGpuPp);

    // ensure stream waits until dependency has been satisfied
    dependencyEvent->enqueueWaitEvent(*haloStream_);

    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuMoveX);
    launchPackXKernel(box);
    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuMoveX);

    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);

    // Consider time spent in communicateHaloData as Comm.X counter
    // TODO: We need further refinement here as communicateHaloData includes launch time for async mem copy
    wallcycle_start(wcycle_, WallCycleCounter::MoveX);


    // wait for remote co-ordinates is implicit with process-MPI as non-local stream is synchronized before MPI calls
    // and MPI_Waitall call makes sure both neighboring ranks' non-local stream is synchronized before data transfer is initiated
    // For multidimensional halo exchanges, this needs to be done for every dimIndex_, since the remote ranks will be different
    // for each. But different pulses within a dimension will communicate with the same remote ranks, so we can restrict to the first pulse.
    if (GMX_THREAD_MPI && pulse_ == 0)
    {
        enqueueWaitRemoteCoordinatesReadyEvent(dependencyEvent);
    }

    if (receiveInPlace_)
    {
        communicateHaloData(d_sendBuf_,
                            0,
                            xSendSize_,
                            sendRankX_,
                            GMX_THREAD_MPI ? remoteXPtr_ : d_x_,
                            GMX_THREAD_MPI ? 0 : atomOffset_,
                            xRecvSize_,
                            recvRankX_,
                            HaloType::Coordinates);
    }
    else
    {
        communicateHaloCoordinatesOutOfPlace(d_sendBuf_, xSendSize_, sendRankX_, xRecvSize_, recvRankX_);
    }

    coordinateHaloLaunched_.markEvent(*haloStream_);

    wallcycle_stop(wcycle_, WallCycleCounter::MoveX);

    return &coordinateHaloLaunched_;
}

// The following method should be called after non-local buffer operations,
// and before the local buffer operations.
void GpuHaloExchange::Impl::communicateHaloForces(bool accumulateForces,
                                                  FixedCapacityVector<GpuEventSynchronizer*, 2>* dependencyEvents)
{

    // Consider time spent in communicateHaloData as Comm.F counter
    // TODO: We need further refinement here as communicateHaloData includes launch time for async mem copy
    wallcycle_start(wcycle_, WallCycleCounter::MoveF);

    while (!dependencyEvents->empty())
    {
        auto* dependency = dependencyEvents->back();
        dependency->enqueueWaitEvent(*haloStream_);
        dependencyEvents->pop_back();
    }

    // Communicate halo data
    if (receiveInPlace_)
    {
        communicateHaloData(d_f_,
                            atomOffset_,
                            fSendSize_,
                            sendRankF_,
                            GMX_THREAD_MPI ? remoteFPtr_ : d_recvBuf_,
                            0,
                            fRecvSize_,
                            recvRankF_,
                            HaloType::Forces);
    }
    else
    {
        communicateHaloForcesOutOfPlace(d_f_, fSendSize_, sendRankF_, fRecvSize_, recvRankF_);
    }

    wallcycle_stop(wcycle_, WallCycleCounter::MoveF);

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuMoveF);

    // Unpack halo buffer into force array
    if (pulse_ > 0 || dd_->ndim > 1)
    {
        // We need to accumulate rather than set, since it is possible
        // that, in this pulse/dim, a value could be written to a location
        // corresponding to the halo region of a following pulse/dim.
        accumulateForces = true;
    }

    launchUnpackFKernel(accumulateForces);

    fReadyOnDevice_.markEvent(*haloStream_);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuMoveF);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

void GpuHaloExchange::Impl::communicateHaloData(DeviceBuffer<Float3> sendPtr,
                                                int                  sendOffset,
                                                int                  sendSize,
                                                int                  sendRank,
                                                DeviceBuffer<Float3> recvPtr,
                                                int                  recvOffset,
                                                int                  recvSize,
                                                int                  recvRank,
                                                HaloType             haloType)
{
    if constexpr (supportedThreadMpiBuild)
    {
        // no need to explicitly sync with GMX_THREAD_MPI as all operations are
        // anyway launched in correct stream
        communicateHaloDataPeerToPeer(&sendPtr, sendOffset, sendSize, sendRank, &recvPtr, recvRank, haloType);
        GMX_RELEASE_ASSERT(recvOffset == 0,
                           "communicateHaloDataPeerToPeer does not support receiveOffset");
    }
    else if constexpr (supportedLibMpiBuild)
    {
        communicateHaloDataGpuAwareMpi(
                asMpiPointer(sendPtr), sendOffset, sendSize, sendRank, asMpiPointer(recvPtr), recvOffset, recvSize, recvRank);
    }
    else
    {
        GMX_RELEASE_ASSERT(false, "Calling GPU haloexchange without MPI");
    }
}

void GpuHaloExchange::Impl::communicateHaloDataGpuAwareMpi(Float3* sendPtr,
                                                           int     sendOffset,
                                                           int     sendSize,
                                                           int     sendRank,
                                                           Float3* recvPtr,
                                                           int     recvOffset,
                                                           int     recvSize,
                                                           int     recvRank)
{
    GMX_RELEASE_ASSERT(supportedLibMpiBuild, "Build does not support GPU-aware MPI");
    // no need to wait for haloDataReadyOnDevice event if this rank is not sending any data
    if (sendSize > 0)
    {
        // wait for halo stream to complete all outstanding
        // activities, to ensure that buffer is up-to-date in GPU memory
        // before transferring to remote rank

        // TODO: Replace stream synchronize with event synchronize
        haloStream_->synchronize();
    }

    // perform halo exchange directly in device buffers
#if GMX_MPI
    MPI_Request request;

    // recv remote data into halo region
    MPI_Irecv(recvPtr + recvOffset, recvSize * DIM, MPI_FLOAT, recvRank, 0, mpi_comm_mysim_, &request);

    // send data to remote halo region
    MPI_Send(sendPtr + sendOffset, sendSize * DIM, MPI_FLOAT, sendRank, 0, mpi_comm_mysim_);

    MPI_Wait(&request, MPI_STATUS_IGNORE);
#else
    GMX_UNUSED_VALUE(sendPtr);
    GMX_UNUSED_VALUE(sendRank);
    GMX_UNUSED_VALUE(recvPtr);
    GMX_UNUSED_VALUE(recvSize);
    GMX_UNUSED_VALUE(recvRank);
#endif
}

// Note that this coordinate and the below force methods involve
// staging data through host memory, which allows a simple
// implementation with unified codepaths for thread-MPI and
// lib-MPI. While it is possible to optimize, this path is only
// expected to be triggered in etremely over-decomposed corner cases
// which have very poor performance for other reasons.
void GpuHaloExchange::Impl::communicateHaloCoordinatesOutOfPlace(DeviceBuffer<Float3> d_sendPtr,
                                                                 int                  sendSize,
                                                                 int                  sendRank,
                                                                 int                  recvSize,
                                                                 int                  recvRank)
{
#if GMX_MPI
    MPI_Request request;
    // copy entire halo buffer to staging send buffer in host memory
    copyFromDeviceBuffer(
            h_outOfPlaceSendBuffer_.data(), &d_sendPtr, 0, sendSize, *haloStream_, GpuApiCallBehavior::Async, nullptr);
    haloStream_->synchronize();
    // exchange host staging buffers with MPI
    MPI_Irecv(h_outOfPlaceRecvBuffer_.data(), recvSize * DIM, MPI_FLOAT, recvRank, 0, mpi_comm_mysim_, &request);
    MPI_Send(h_outOfPlaceSendBuffer_.data(), sendSize * DIM, MPI_FLOAT, sendRank, 0, mpi_comm_mysim_);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    // copy sections of staging receive buffer, in turn, to final locations in device memory
    int stageBufIndex = 0;
    for (int zone = 0; zone < numZone_; zone++)
    {
        int numElements = ind_->cell2at1[zone] - ind_->cell2at0[zone];
        copyToDeviceBuffer(&d_x_,
                           &h_outOfPlaceRecvBuffer_[stageBufIndex],
                           ind_->cell2at0[zone],
                           numElements,
                           *haloStream_,
                           GpuApiCallBehavior::Async,
                           nullptr);
        stageBufIndex += numElements;
    }
    haloStream_->synchronize();
#else
    GMX_UNUSED_VALUE(d_sendPtr);
    GMX_UNUSED_VALUE(sendSize);
    GMX_UNUSED_VALUE(sendRank);
    GMX_UNUSED_VALUE(recvSize);
    GMX_UNUSED_VALUE(recvRank);
#endif
}

void GpuHaloExchange::Impl::communicateHaloForcesOutOfPlace(DeviceBuffer<Float3> d_sendPtr,
                                                            int                  sendSize,
                                                            int                  sendRank,
                                                            int                  recvSize,
                                                            int                  recvRank)
{
#if GMX_MPI
    MPI_Request request;
    // copy sections of force buffer in device memory, in turn, to host staging buffer
    int stageBufIndex = 0;
    for (int zone = 0; zone < numZone_; zone++)
    {
        int numElements = ind_->cell2at1[zone] - ind_->cell2at0[zone];
        copyFromDeviceBuffer(&h_outOfPlaceSendBuffer_[stageBufIndex],
                             &d_sendPtr,
                             ind_->cell2at0[zone],
                             numElements,
                             *haloStream_,
                             GpuApiCallBehavior::Async,
                             nullptr);
        stageBufIndex += numElements;
    }
    haloStream_->synchronize();
    // exchange host staging buffers with MPI
    MPI_Irecv(h_outOfPlaceRecvBuffer_.data(), recvSize * DIM, MPI_FLOAT, recvRank, 0, mpi_comm_mysim_, &request);
    MPI_Send(h_outOfPlaceSendBuffer_.data(), sendSize * DIM, MPI_FLOAT, sendRank, 0, mpi_comm_mysim_);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    // copy entire host staging receive buffer to device memory receive buffer
    copyToDeviceBuffer(
            &d_recvBuf_, h_outOfPlaceRecvBuffer_.data(), 0, recvSize, *haloStream_, GpuApiCallBehavior::Async, nullptr);
#else
    GMX_UNUSED_VALUE(d_sendPtr);
    GMX_UNUSED_VALUE(sendSize);
    GMX_UNUSED_VALUE(sendRank);
    GMX_UNUSED_VALUE(recvSize);
    GMX_UNUSED_VALUE(recvRank);
#endif
}

void GpuHaloExchange::Impl::communicateHaloDataPeerToPeer(DeviceBuffer<Float3>* sendPtr,
                                                          int                   sendOffset,
                                                          int                   sendSize,
                                                          int                   sendRank,
                                                          DeviceBuffer<Float3>* remotePtr,
                                                          int                   recvRank,
                                                          HaloType              haloType)
{
    GMX_RELEASE_ASSERT(supportedThreadMpiBuild, "Build does not support peer-to-peer communication");
#if GMX_GPU_CUDA || GMX_GPU_HIP
    DeviceBuffer<Float3> offsetSendPtr = *sendPtr + sendOffset;
    copyBetweenDeviceBuffers(
            remotePtr, &offsetSendPtr, sendSize, *haloStream_, GpuApiCallBehavior::Async, nullptr);

#    if GMX_THREAD_MPI
    // ensure pushed data has arrived before remote rank progresses
    // This rank records an event and sends it to the remote rank which has just been pushed data.
    // This rank receives event from remote rank which has pushed data here, and enqueues that event
    // to its stream.
    GpuEventSynchronizer* haloDataTransferRemote;
    GpuEventSynchronizer* haloDataTransferLaunched = (haloType == HaloType::Coordinates)
                                                             ? haloXDataTransferLaunched_.get()
                                                             : haloFDataTransferLaunched_.get();
    GMX_ASSERT(haloDataTransferLaunched != nullptr,
               "Halo exchange requires valid event to synchronize data transfer initiated in "
               "remote rank");
    haloDataTransferLaunched->markEvent(*haloStream_);

    MPI_Sendrecv(&haloDataTransferLaunched,
                 sizeof(GpuEventSynchronizer*), //NOLINT(bugprone-sizeof-expression)
                 MPI_BYTE,
                 sendRank,
                 0,
                 &haloDataTransferRemote,
                 sizeof(GpuEventSynchronizer*), //NOLINT(bugprone-sizeof-expression)
                 MPI_BYTE,
                 recvRank,
                 0,
                 mpi_comm_mysim_,
                 MPI_STATUS_IGNORE);

    haloDataTransferRemote->enqueueWaitEvent(*haloStream_);
#    else
    GMX_UNUSED_VALUE(sendRank);
    GMX_UNUSED_VALUE(recvRank);
    GMX_UNUSED_VALUE(haloType);
#    endif
#else
    GMX_UNUSED_VALUE(sendPtr);
    GMX_UNUSED_VALUE(sendOffset);
    GMX_UNUSED_VALUE(sendSize);
    GMX_UNUSED_VALUE(sendRank);
    GMX_UNUSED_VALUE(remotePtr);
    GMX_UNUSED_VALUE(recvRank);
    GMX_UNUSED_VALUE(haloType);
#endif
}

GpuEventSynchronizer* GpuHaloExchange::Impl::getForcesReadyOnDeviceEvent()
{
    return &fReadyOnDevice_;
}

/*! \brief Create GpuHaloExchange object */
GpuHaloExchange::Impl::Impl(gmx_domdec_t*        dd,
                            int                  dimIndex,
                            MPI_Comm             mpi_comm_mysim,
                            MPI_Comm             mpi_comm_mysim_world,
                            const DeviceContext& deviceContext,
                            int                  pulse,
                            gmx_wallcycle*       wcycle) :
    dd_(dd),
    sendRankX_(dd->neighbor[dimIndex][1]),
    recvRankX_(dd->neighbor[dimIndex][0]),
    sendRankF_(dd->neighbor[dimIndex][0]),
    recvRankF_(dd->neighbor[dimIndex][1]),
    usePBC_(dd->ci[dd->dim[dimIndex]] == 0),
    haloXDataTransferLaunched_(GMX_THREAD_MPI ? new GpuEventSynchronizer() : nullptr),
    haloFDataTransferLaunched_(GMX_THREAD_MPI ? new GpuEventSynchronizer() : nullptr),
    mpi_comm_mysim_(mpi_comm_mysim),
    mpi_comm_mysim_world_(mpi_comm_mysim_world),
    deviceContext_(deviceContext),
    haloStream_(new DeviceStream(deviceContext, DeviceStreamPriority::High, false)),
    dimIndex_(dimIndex),
    pulse_(pulse),
    wcycle_(wcycle)
{
    if (usePBC_ && dd->unitCellInfo.haveScrewPBC)
    {
        gmx_fatal(FARGS, "Error: screw is not yet supported in GPU halo exchange\n");
    }

    changePinningPolicy(&h_indexMap_, gmx::PinningPolicy::PinnedIfSupported);

    allocateDeviceBuffer(&d_fShift_, 1, deviceContext_);
}

GpuHaloExchange::Impl::~Impl()
{
    try
    {
        freeDeviceBuffer(&d_recvBuf_);
        freeDeviceBuffer(&d_indexMap_);
        freeDeviceBuffer(&d_sendBuf_);
        freeDeviceBuffer(&d_fShift_);
    }
    catch (gmx::InternalError& e)
    {
        fprintf(stderr, "Internal error in destructor of GpuHaloExchange: %s\n", e.what());
    }
}

GpuHaloExchange::GpuHaloExchange(gmx_domdec_t*        dd,
                                 int                  dimIndex,
                                 MPI_Comm             mpi_comm_mysim,
                                 MPI_Comm             mpi_comm_mysim_world_,
                                 const DeviceContext& deviceContext,
                                 int                  pulse,
                                 gmx_wallcycle*       wcycle) :
    impl_(new Impl(dd, dimIndex, mpi_comm_mysim, mpi_comm_mysim_world_, deviceContext, pulse, wcycle))
{
}

GpuHaloExchange::GpuHaloExchange(GpuHaloExchange&&) noexcept = default;

GpuHaloExchange& GpuHaloExchange::operator=(GpuHaloExchange&& other) noexcept
{
    std::swap(impl_, other.impl_);
    return *this;
}

GpuHaloExchange::~GpuHaloExchange() = default;

void GpuHaloExchange::reinitHalo(DeviceBuffer<RVec> d_coordinatesBuffer, DeviceBuffer<RVec> d_forcesBuffer)
{
    impl_->reinitHalo(d_coordinatesBuffer, d_forcesBuffer);
}

GpuEventSynchronizer* GpuHaloExchange::communicateHaloCoordinates(const matrix box,
                                                                  GpuEventSynchronizer* dependencyEvent)
{
    return impl_->communicateHaloCoordinates(box, dependencyEvent);
}

void GpuHaloExchange::communicateHaloForces(bool accumulateForces,
                                            FixedCapacityVector<GpuEventSynchronizer*, 2>* dependencyEvents)
{
    impl_->communicateHaloForces(accumulateForces, dependencyEvents);
}

GpuEventSynchronizer* GpuHaloExchange::getForcesReadyOnDeviceEvent()
{
    return impl_->getForcesReadyOnDeviceEvent();
}

GpuHaloExchangeNvshmemHelper::GpuHaloExchangeNvshmemHelper(const gmx_domdec_t&       dd,
                                                           const DeviceContext&      context,
                                                           const DeviceStream&       stream,
                                                           const std::optional<int>& peerRank,
                                                           gmx_wallcycle*            wcycle,
                                                           MPI_Comm mpi_comm_mygroup,
                                                           MPI_Comm mpi_comm_mysim_world) :
    dd_(dd), stream_(stream), peerRank_(peerRank), context_(context), wcycle_(wcycle)
{
#if GMX_NVSHMEM
    fusedPpHaloExchange_ = std::make_unique<gmx::FusedGpuHaloExchange>(
            context_, wcycle, mpi_comm_mygroup, mpi_comm_mysim_world);
#else
    GMX_UNUSED_VALUE(mpi_comm_mygroup);
    GMX_UNUSED_VALUE(mpi_comm_mysim_world);
#endif
}

GpuHaloExchangeNvshmemHelper::~GpuHaloExchangeNvshmemHelper()
{
    if (ppHaloExSyncBufCapacity_ > 0)
    {
        freeDeviceBuffer(&d_ppHaloExSyncBase_);
    }
}
// Fused pass-through API definitions to avoid incomplete-type usage in header
void GpuHaloExchangeNvshmemHelper::reinitAllHaloExchanges(const t_commrec&   cr,
                                                          DeviceBuffer<RVec> d_coordinatesBuffer,
                                                          DeviceBuffer<RVec> d_forcesBuffer)
{
#if GMX_NVSHMEM
    fusedPpHaloExchange_->reinitAllHaloExchanges(
            cr, d_coordinatesBuffer, d_forcesBuffer, d_ppHaloExSyncBase_, totalNumPulses_);
#else
    GMX_UNUSED_VALUE(cr);
    GMX_UNUSED_VALUE(d_coordinatesBuffer);
    GMX_UNUSED_VALUE(d_forcesBuffer);
#endif
}

GpuEventSynchronizer* GpuHaloExchangeNvshmemHelper::launchAllCoordinateExchanges(const matrix box,
                                                                                 GpuEventSynchronizer* dependencyEvent)
{
#if GMX_NVSHMEM
    return fusedPpHaloExchange_->launchAllCoordinateExchanges(box, dependencyEvent);
#else
    GMX_UNUSED_VALUE(box);
    return dependencyEvent;
#endif
}

GpuEventSynchronizer* GpuHaloExchangeNvshmemHelper::launchAllForceExchanges(
        bool                                           accumulateForces,
        FixedCapacityVector<GpuEventSynchronizer*, 2>* dependencyEvents)
{
#if GMX_NVSHMEM
    return fusedPpHaloExchange_->launchAllForceExchanges(accumulateForces, dependencyEvents);
#else
    GMX_UNUSED_VALUE(accumulateForces);
    GMX_UNUSED_VALUE(dependencyEvents);
    return nullptr;
#endif
}

void GpuHaloExchangeNvshmemHelper::destroyAllHaloExchangeBuffers()
{
#if GMX_NVSHMEM
    fusedPpHaloExchange_->destroyAllHaloExchangeBuffers();
#endif
}

GpuEventSynchronizer* GpuHaloExchangeNvshmemHelper::getForcesReadyOnDeviceEvent()
{
#if GMX_NVSHMEM
    return fusedPpHaloExchange_->getForcesReadyOnDeviceEvent();
#else
    return nullptr;
#endif
}

DeviceBuffer<uint64_t> GpuHaloExchangeNvshmemHelper::getSyncBuffer() const
{
    return d_ppHaloExSyncBase_;
}

int GpuHaloExchangeNvshmemHelper::totalPulsesAndDims() const
{
    return totalNumPulses_;
}

void GpuHaloExchangeNvshmemHelper::allocateAndInitSignalBufs(int totalNumPulses)
{
    int totalSyncBufSize = totalNumPulses * numOfPpHaloExSyncBufs;
    reallocateDeviceBuffer(
            &d_ppHaloExSyncBase_, totalSyncBufSize, &ppHaloExSyncBufSize_, &ppHaloExSyncBufCapacity_, context_, true);
    // If the num of dims/pulses have changed we initialize the signalling
    // buffer to max val.
    if (ppHaloExPerSyncBufSize_ < totalNumPulses)
    {
        // Initialize the signalling buffer with max value.
        ppHaloExPerSyncBufSize_ = totalNumPulses;

        // TODO host allocation should be minimized by making this a
        // member variable
        gmx::HostVector<uint64_t> hostBuffer = {
            {}, gmx::HostAllocationPolicy(gmx::PinningPolicy::PinnedIfSupported)
        };
        hostBuffer.resize(totalSyncBufSize, ~0);
        // TODO replace this D2H with cudaMemsetAsync
        copyToDeviceBuffer<uint64_t>(&d_ppHaloExSyncBase_,
                                     hostBuffer.data(),
                                     0,
                                     static_cast<size_t>(totalSyncBufSize),
                                     stream_,
                                     GpuApiCallBehavior::Async,
                                     nullptr);
#if GMX_NVSHMEM
        nvshmemx_sync_all_on_stream(stream_.stream());
#endif
    }
}

void GpuHaloExchangeNvshmemHelper::reinit()
{
    const bool isPmeRank      = peerRank_.has_value();
    int        totalNumPulses = 0;
    if (isPmeRank)
    {
#if GMX_MPI
        // recv remote data of num of pulses in each dim.
        MPI_Recv(&totalNumPulses, 1, MPI_INT, peerRank_.value(), 0, dd_.mpiCommMySim().comm(), MPI_STATUS_IGNORE);
#endif
    }
    else
    {
        for (int d = 0; d < dd_.ndim; d++)
        {
            totalNumPulses += dd_.comm->cd[d].numPulses();
        }
#if GMX_MPI
        // we use PP rank which receives virial and energy from PME rank
        // to send the number of pulses data to PME rank.
        if (dd_.pme_receive_vir_ener)
        {
            MPI_Send(&totalNumPulses, 1, MPI_INT, dd_.pme_nodeid, 0, dd_.mpiCommMySim().comm());
        }
#endif
    }

    if (totalNumPulses > 0)
    {
        // These allocations are required only for PP Haloexchange
        // which is enabled if totalNumPulses > 0
        if (totalNumPulses_ < totalNumPulses)
        {
            totalNumPulses_ = totalNumPulses;
            allocateAndInitSignalBufs(totalNumPulses);
        }

        if (isPmeRank)
        {
            int64_t newSize              = 1;
            int64_t totalSendRecvBufSize = 1;

            MPI_Allreduce(
                    &newSize, &totalSendRecvBufSize, 1, MPI_INT64_T, MPI_MAX, dd_.mpiCommMySim().comm());
            size_t totalSendSize = static_cast<size_t>(totalSendRecvBufSize);
            size_t totalRecvSize = static_cast<size_t>(totalSendRecvBufSize);

            reallocateDeviceBuffer(&d_unifiedSendBufForPpFromPme_,
                                   totalSendSize,
                                   &d_unifiedSendBufSizeForPpFromPme_,
                                   &d_unifiedSendBufCapacityForPpFromPme_,
                                   context_,
                                   true);
            reallocateDeviceBuffer(&d_unifiedRecvBufForPpFromPme_,
                                   totalRecvSize,
                                   &d_unifiedRecvBufSizeForPpFromPme_,
                                   &d_unifiedRecvBufCapacityForPpFromPme_,
                                   context_,
                                   true);
        }
    }
}

} // namespace gmx
