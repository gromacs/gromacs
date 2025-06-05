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

#if GMX_NVSHMEM
#    include <nvshmem.h>
#    include <nvshmemx.h>
#endif

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/gpu_utils/nvshmem_utils.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxmpi.h"

#include "domdec_internal.h"

// NOLINTNEXTLINE(misc-redundant-expression)
constexpr bool supportedLibMpiBuild =
        ((GMX_LIB_MPI != 0) && ((GMX_GPU_CUDA != 0) || (GMX_GPU_SYCL != 0)));
// NOLINTNEXTLINE(misc-redundant-expression)
constexpr bool supportedThreadMpiBuild = ((GMX_THREAD_MPI != 0) && (GMX_GPU_CUDA != 0));

namespace gmx
{

void GpuHaloExchange::Impl::reinitHalo(DeviceBuffer<Float3> d_coordinatesBuffer,
                                       DeviceBuffer<Float3> d_forcesBuffer)
{
    GMX_RELEASE_ASSERT(supportedLibMpiBuild || supportedThreadMpiBuild,
                       "Gpu Halo Exchange not supported in this build");

    wallcycle_start_nocount(wcycle_, WallCycleCounter::Domdec);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::DDGpu);

    d_x_ = d_coordinatesBuffer;
    d_f_ = d_forcesBuffer;

    const gmx_domdec_comm_t&     comm = *dd_->comm;
    const gmx_domdec_comm_dim_t& cd   = comm.cd[dimIndex_];
    ind_                              = &cd.ind[pulse_];
    receiveInPlace_                   = cd.receiveInPlace;

    numHomeAtoms_ = comm.atomRanges.numHomeAtoms(); // offset for data received by this rank

    // Determine receive offset for the dimension index and pulse of this halo exchange object
    int numZoneTemp   = 1;
    int numAtomsTotal = numHomeAtoms_;
    for (int i = 0; i <= dimIndex_; i++)
    {
        GMX_ASSERT(numZoneTemp <= sc_maxNumIZones, "Too many domain decomposition zones");
        int pulseMax = (i == dimIndex_) ? pulse_ : (comm.cd[i].numPulses() - 1);
        for (int p = 0; p <= pulseMax; p++)
        {
            atomOffset_                     = numAtomsTotal;
            const gmx_domdec_ind_t& indTemp = comm.cd[i].ind[p];
            numAtomsTotal += indTemp.nrecv[numZoneTemp + 1];
        }
        numZone_ = numZoneTemp;
        numZoneTemp += numZoneTemp;
    }

    int newSize = ind_->nsend[numZone_ + 1];

    // reallocates only if needed
    h_indexMap_.resize(newSize);
    // reallocate on device only if needed
    if (newSize > maxPackedBufferSize_)
    {
        reallocateDeviceBuffer(&d_indexMap_, newSize, &indexMapSize_, &indexMapSizeAlloc_, deviceContext_);

        if (useNvshmem_ && (maxPackedBufferSize_ > 0) && (newSize > sendBufSizeAlloc_))
        {
            GMX_RELEASE_ASSERT(d_sendBuf_ != nullptr,
                               "Halo exchange requires d_sendBuf_ buffer to be valid");
#if GMX_NVSHMEM
            // unregister only when previous d_sendBuf_ was registered previously
            // via nvshmemx_buffer_register
            GMX_RELEASE_ASSERT(nvshmemx_buffer_unregister(d_sendBuf_) == 0,
                               "NVSHMEM d_sendBuf_ Buffer unregistration failed");
#endif
        }
        reallocateDeviceBuffer(&d_sendBuf_, newSize, &sendBufSize_, &sendBufSizeAlloc_, deviceContext_);

        if (useNvshmem_ && (maxPackedBufferSize_ < sendBufSizeAlloc_))
        {
            // As d_sendBuf_ is a source buffer in the PP Halo exchange nvshmem_put
            // we do not need to do a symmetric allocation for it, registering it via
            // nvshmemx_buffer_register is sufficient.
#if GMX_NVSHMEM
            std::size_t bufLen = sendBufSizeAlloc_ * sizeof(float3);
            GMX_RELEASE_ASSERT(nvshmemx_buffer_register(d_sendBuf_, bufLen) == 0,
                               "NVSHMEM d_sendBuf_ Buffer registration failed");
#endif
        }
        // number of values/elems is same for indexMapSizeAlloc_/sendBufSizeAlloc_ so we can use either.
        maxPackedBufferSize_ = sendBufSizeAlloc_;
    }

    xSendSize_         = newSize;
    int recvBufNewSize = newSize;
    if (useNvshmem_)
    {
        reinitXGridSizeAndDevBarrier();
#if GMX_MPI
        MPI_Allreduce(&newSize, &recvBufNewSize, 1, MPI_INT, MPI_MAX, mpi_comm_mysim_world_);
        // remote PE atomOffset to nvshmem put halo coordinates
        MPI_Sendrecv(&atomOffset_,
                     sizeof(int),
                     MPI_BYTE,
                     recvRankX_,
                     0,
                     &nvshmemHaloExchange_.putAtomOffsetInReceiverRankXBuf_,
                     sizeof(int),
                     MPI_BYTE,
                     sendRankX_,
                     0,
                     mpi_comm_mysim_,
                     MPI_STATUS_IGNORE);
#endif
    }

    reallocateDeviceBuffer(
            &d_recvBuf_, recvBufNewSize, &recvBufSize_, &recvBufSizeAlloc_, deviceContext_, useNvshmem_);

#if GMX_MPI
    MPI_Sendrecv(&xSendSize_,
                 sizeof(int),
                 MPI_BYTE,
                 sendRankX_,
                 0,
                 &xRecvSize_,
                 sizeof(int),
                 MPI_BYTE,
                 recvRankX_,
                 0,
                 mpi_comm_mysim_,
                 MPI_STATUS_IGNORE);
#endif
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
    if (supportedThreadMpiBuild)
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

void GpuHaloExchange::Impl::reinitNvshmemSignal(const t_commrec& cr, int signalObjOffset)
{
    if (useNvshmem_)
    {
        GMX_RELEASE_ASSERT(cr.nvshmemHandlePtr->d_ppHaloExSyncBase_ != nullptr,
                           "NVSHMEM Coordinate Halo exchange requires valid signal buffer");
        nvshmemHaloExchange_.signalObjOffset_     = signalObjOffset;
        nvshmemHaloExchange_.d_signalSenderRankX_ = cr.nvshmemHandlePtr->d_ppHaloExSyncBase_;
        // As only CUDA DeviceBuffer<> supports pointer updates from host side
        // we guard these pointer update code by GMX_GPU_CUDA
#if GMX_GPU_CUDA
        int totalPulsesAndDims = cr.nvshmemHandlePtr->ppHaloExPerSyncBufSize_;
        nvshmemHaloExchange_.d_signalReceiverRankX_ =
                cr.nvshmemHandlePtr->d_ppHaloExSyncBase_ + totalPulsesAndDims;
        nvshmemHaloExchange_.d_signalReceiverRankF_ =
                cr.nvshmemHandlePtr->d_ppHaloExSyncBase_ + 2 * totalPulsesAndDims;
#endif
    }
}

void GpuHaloExchange::Impl::reinitXGridSizeAndDevBarrier()
{
    if (useNvshmem_)
    {
        // Add 1 to size to take care of case when xSendSize_ == 0 in such case
        // we still launch the kernel to signal the remote rank recvRankX_
        // from which this rank receives packed data as (xRecvSize_ > 0)
        nvshmemHaloExchange_.gridDimX_ =
                gmx::divideRoundUp(xSendSize_ + 1, nvshmemHaloExchange_.c_nvshmemThreadsPerBlock);

        if ((nvshmemHaloExchange_.arriveWaitBarrierVal_ != nvshmemHaloExchange_.gridDimX_)
            && (xSendSize_ > 0))
        {
            nvshmemHaloExchange_.arriveWaitBarrierVal_ = nvshmemHaloExchange_.gridDimX_;
#if GMX_NVSHMEM
            GMX_RELEASE_ASSERT(
                    nvshmemHaloExchange_.d_arriveWaitBarrier_ != nullptr,
                    "NVSHMEM Coordinate Halo exchange requires valid arrive wait barrier buffer");
            cuda::barrier<cuda::thread_scope_device> temp_bar(nvshmemHaloExchange_.arriveWaitBarrierVal_);
            cudaMemcpyAsync(nvshmemHaloExchange_.d_arriveWaitBarrier_,
                            &temp_bar,
                            sizeof(cuda::barrier<cuda::thread_scope_device>),
                            cudaMemcpyDefault,
                            haloStream_->stream());
#endif
        }
    }
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

    if (!useNvshmem_)
    {
        // wait for remote co-ordinates is implicit with process-MPI as non-local stream is synchronized before MPI calls
        // and MPI_Waitall call makes sure both neighboring ranks' non-local stream is synchronized before data transfer is initiated
        // For multidimensional halo exchanges, this needs to be done for every dimIndex_, since the remote ranks will be different
        // for each. But different pulses within a dimension will communicate with the same remote ranks, so we can restrict to the first pulse.
        if (GMX_THREAD_MPI && pulse_ == 0)
        {
            enqueueWaitRemoteCoordinatesReadyEvent(dependencyEvent);
        }

        Float3* recvPtr = GMX_THREAD_MPI ? asMpiPointer(remoteXPtr_) : (asMpiPointer(d_x_) + atomOffset_);

        if (receiveInPlace_)
        {
            communicateHaloData(
                    asMpiPointer(d_sendBuf_), xSendSize_, sendRankX_, recvPtr, xRecvSize_, recvRankX_, HaloType::Coordinates);
        }
        else
        {
            communicateHaloCoordinatesOutOfPlace(d_sendBuf_, xSendSize_, sendRankX_, xRecvSize_, recvRankX_);
        }
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

    Float3* recvPtr = asMpiPointer(GMX_THREAD_MPI ? remoteFPtr_ : d_recvBuf_);

    if (!useNvshmem_ || (useNvshmem_ && !receiveInPlace_))
    {
        // Communicate halo data
        if (receiveInPlace_)
        {
            communicateHaloData((asMpiPointer(d_f_) + atomOffset_),
                                fSendSize_,
                                sendRankF_,
                                recvPtr,
                                fRecvSize_,
                                recvRankF_,
                                HaloType::Forces);
        }
        else
        {
            communicateHaloForcesOutOfPlace(d_f_, fSendSize_, sendRankF_, fRecvSize_, recvRankF_);
        }
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

void GpuHaloExchange::Impl::communicateHaloData(Float3*  sendPtr,
                                                int      sendSize,
                                                int      sendRank,
                                                Float3*  recvPtr,
                                                int      recvSize,
                                                int      recvRank,
                                                HaloType haloType)
{
    if (supportedThreadMpiBuild)
    {
        // no need to explicitly sync with GMX_THREAD_MPI as all operations are
        // anyway launched in correct stream
        communicateHaloDataPeerToPeer(sendPtr, sendSize, sendRank, recvPtr, recvRank, haloType);
    }
    else
    {
        communicateHaloDataGpuAwareMpi(sendPtr, sendSize, sendRank, recvPtr, recvSize, recvRank);
    }
}

void GpuHaloExchange::Impl::communicateHaloDataGpuAwareMpi(Float3* sendPtr,
                                                           int     sendSize,
                                                           int     sendRank,
                                                           Float3* recvPtr,
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
    MPI_Irecv(recvPtr, recvSize * DIM, MPI_FLOAT, recvRank, 0, mpi_comm_mysim_, &request);

    // send data to remote halo region
    MPI_Send(sendPtr, sendSize * DIM, MPI_FLOAT, sendRank, 0, mpi_comm_mysim_);

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

void GpuHaloExchange::Impl::communicateHaloDataPeerToPeer(Float3*  sendPtr,
                                                          int      sendSize,
                                                          int      sendRank,
                                                          Float3*  remotePtr,
                                                          int      recvRank,
                                                          HaloType haloType)
{
    GMX_RELEASE_ASSERT(supportedThreadMpiBuild, "Build does not support peer-to-peer communication");
    // The code below can be made backend-agnostic once we add device-to-device copy functionality to DeviceBuffer
#if GMX_GPU_CUDA
    cudaError_t stat;

    // We asynchronously push data to remote rank. The remote
    // destination pointer has already been set in the init fn.  We
    // don't need to worry about overwriting data the remote ranks
    // still needs since the halo exchange is just done once per
    // timestep, for each of X and F.

    // send data to neighbor, if any data exists to send
    if (sendSize > 0)
    {
        stat = cudaMemcpyAsync(remotePtr,
                               sendPtr,
                               sendSize * DIM * sizeof(float),
                               cudaMemcpyDeviceToDevice,
                               haloStream_->stream());

        CU_RET_ERR(stat, "cudaMemcpyAsync on GPU Domdec CUDA direct data transfer failed");
    }

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
                            bool                 useNvshmem,
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
    wcycle_(wcycle),
    useNvshmem_(useNvshmem)
{
    if (usePBC_ && dd->unitCellInfo.haveScrewPBC)
    {
        gmx_fatal(FARGS, "Error: screw is not yet supported in GPU halo exchange\n");
    }

    changePinningPolicy(&h_indexMap_, gmx::PinningPolicy::PinnedIfSupported);

    allocateDeviceBuffer(&d_fShift_, 1, deviceContext_);

    if (useNvshmem_)
    {
#if GMX_NVSHMEM
        allocateDeviceBuffer(&nvshmemHaloExchange_.d_arriveWaitBarrier_, 1, deviceContext_);
#endif
    }
}

GpuHaloExchange::Impl::~Impl()
{
    if (useNvshmem_)
    {
#if GMX_NVSHMEM
        // As d_sendBuf_ is a source buffer in the PP Halo exchange nvshmem_put
        // we had registered it via nvshmemx_buffer_register, such registered buffer
        // need to be first unregistered via nvshmemx_buffer_unregister before freeing.
        if (d_sendBuf_)
        {
            GMX_RELEASE_ASSERT(nvshmemx_buffer_unregister(d_sendBuf_) == 0,
                               "NVSHMEM d_sendBuf_ Buffer unregistration failed");
        }
#endif
    }
    else
    {
        // For the NVSHMEM path the freeing of d_recvBuf_
        // happens in destroyGpuHaloExchangeNvshmemBuf() due to it
        // been a collective call calling it at this point is not appropriate
        freeDeviceBuffer(&d_recvBuf_);
    }

    freeDeviceBuffer(&d_indexMap_);
    freeDeviceBuffer(&d_sendBuf_);
    freeDeviceBuffer(&d_fShift_);
}

void GpuHaloExchange::Impl::destroyGpuHaloExchangeNvshmemBuf()
{
    // freeing the NVSHMEM symmetric buffer
    freeDeviceBuffer(&d_recvBuf_);
}

GpuHaloExchange::GpuHaloExchange(gmx_domdec_t*        dd,
                                 int                  dimIndex,
                                 MPI_Comm             mpi_comm_mysim,
                                 MPI_Comm             mpi_comm_mysim_world_,
                                 const DeviceContext& deviceContext,
                                 int                  pulse,
                                 bool                 useNvshmem,
                                 gmx_wallcycle*       wcycle) :
    impl_(new Impl(dd, dimIndex, mpi_comm_mysim, mpi_comm_mysim_world_, deviceContext, pulse, useNvshmem, wcycle))
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

void GpuHaloExchange::reinitNvshmemSignal(const t_commrec& cr, int signalObjOffset)
{
    impl_->reinitNvshmemSignal(cr, signalObjOffset);
}

void GpuHaloExchange::destroyGpuHaloExchangeNvshmemBuf()
{
    impl_->destroyGpuHaloExchangeNvshmemBuf();
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
} // namespace gmx
