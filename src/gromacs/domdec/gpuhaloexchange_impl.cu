/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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
 * \brief Implements GPU halo exchange using CUDA.
 *
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_domdec
 */
#include "gmxpre.h"

#include "gpuhaloexchange_impl.cuh"

#include "config.h"

#include <assert.h>
#include <stdio.h>

#include <utility>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/gpu_utils/typecasts.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxmpi.h"

#include "domdec_internal.h"

namespace gmx
{

//! Number of CUDA threads in a block
// TODO Optimize this through experimentation
constexpr static int c_threadsPerBlock = 256;

template<bool usePBC>
__global__ void packSendBufKernel(float3* __restrict__ dataPacked,
                                  const float3* __restrict__ data,
                                  const int* __restrict__ map,
                                  const int    mapSize,
                                  const float3 coordinateShift)
{
    int           threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    float3*       gm_dataDest = &dataPacked[threadIndex];
    const float3* gm_dataSrc  = &data[map[threadIndex]];

    if (threadIndex < mapSize)
    {
        if (usePBC)
        {
            *gm_dataDest = *gm_dataSrc + coordinateShift;
        }
        else
        {
            *gm_dataDest = *gm_dataSrc;
        }
    }
}

/*! \brief unpack non-local force data buffer on the GPU using pre-populated "map" containing index
 * information
 * \param[out] data        full array of force values
 * \param[in]  dataPacked  packed array of force values to be transferred
 * \param[in]  map         array of indices defining mapping from full to packed array
 * \param[in]  mapSize     number of elements in map array
 */
template<bool accumulate>
__global__ void unpackRecvBufKernel(float3* __restrict__ data,
                                    const float3* __restrict__ dataPacked,
                                    const int* __restrict__ map,
                                    const int mapSize)
{

    int           threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    const float3* gm_dataSrc  = &dataPacked[threadIndex];
    float3*       gm_dataDest = &data[map[threadIndex]];

    if (threadIndex < mapSize)
    {
        if (accumulate)
        {
            *gm_dataDest += *gm_dataSrc;
        }
        else
        {
            *gm_dataDest = *gm_dataSrc;
        }
    }
}

void GpuHaloExchange::Impl::reinitHalo(float3* d_coordinatesBuffer, float3* d_forcesBuffer)
{
    wallcycle_start(wcycle_, WallCycleCounter::Domdec);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::DDGpu);

    d_x_ = d_coordinatesBuffer;
    d_f_ = d_forcesBuffer;

    const gmx_domdec_comm_t&     comm = *dd_->comm;
    const gmx_domdec_comm_dim_t& cd   = comm.cd[dimIndex_];
    const gmx_domdec_ind_t&      ind  = cd.ind[pulse_];

    numHomeAtoms_ = comm.atomRanges.numHomeAtoms(); // offset for data received by this rank

    // Determine receive offset for the dimension index and pulse of this halo exchange object
    int numZoneTemp   = 1;
    int numZone       = 0;
    int numAtomsTotal = numHomeAtoms_;
    for (int i = 0; i <= dimIndex_; i++)
    {
        int pulseMax = (i == dimIndex_) ? pulse_ : (comm.cd[i].numPulses() - 1);
        for (int p = 0; p <= pulseMax; p++)
        {
            atomOffset_                     = numAtomsTotal;
            const gmx_domdec_ind_t& indTemp = comm.cd[i].ind[p];
            numAtomsTotal += indTemp.nrecv[numZoneTemp + 1];
        }
        numZone = numZoneTemp;
        numZoneTemp += numZoneTemp;
    }

    int newSize = ind.nsend[numZone + 1];

    GMX_ASSERT(cd.receiveInPlace, "Out-of-place receive is not yet supported in GPU halo exchange");

    // reallocates only if needed
    h_indexMap_.resize(newSize);
    // reallocate on device only if needed
    if (newSize > maxPackedBufferSize_)
    {
        reallocateDeviceBuffer(&d_indexMap_, newSize, &indexMapSize_, &indexMapSizeAlloc_, deviceContext_);
        reallocateDeviceBuffer(&d_sendBuf_, newSize, &sendBufSize_, &sendBufSizeAlloc_, deviceContext_);
        reallocateDeviceBuffer(&d_recvBuf_, newSize, &recvBufSize_, &recvBufSizeAlloc_, deviceContext_);
        maxPackedBufferSize_ = newSize;
    }

    xSendSize_ = newSize;
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

    if (newSize > 0)
    {
        GMX_ASSERT(ind.index.size() == h_indexMap_.size(),
                   "Size mismatch between domain decomposition communication index array and GPU "
                   "halo exchange index mapping array");
        std::copy(ind.index.begin(), ind.index.end(), h_indexMap_.begin());

        copyToDeviceBuffer(
                &d_indexMap_, h_indexMap_.data(), 0, newSize, *haloStream_, GpuApiCallBehavior::Async, nullptr);
    }

#if GMX_MPI
    // Exchange of remote addresses from neighboring ranks is needed only with CUDA-direct as cudamemcpy needs both src/dst pointer
    // MPI calls such as MPI_send doesn't worry about receiving address, that is taken care by MPI_recv call in neighboring rank
    if (GMX_THREAD_MPI)
    {
        // This rank will push data to its neighbor, so needs to know
        // the remote receive address and similarly send its receive
        // address to other neighbour. We can do this here in reinit fn
        // since the pointers will not change until the next NS step.

        // Coordinates buffer:
        float3* recvPtr = &d_x_[atomOffset_];
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
        recvPtr = d_recvBuf_;
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
}

GpuEventSynchronizer* GpuHaloExchange::Impl::communicateHaloCoordinates(const matrix box,
                                                                        GpuEventSynchronizer* dependencyEvent)
{
    wallcycle_start(wcycle_, WallCycleCounter::LaunchGpu);

    // ensure stream waits until dependency has been satisfied
    dependencyEvent->enqueueWaitEvent(*haloStream_);

    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuMoveX);

    // launch kernel to pack send buffer
    KernelLaunchConfig config;
    config.blockSize[0]     = c_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = (xSendSize_ + c_threadsPerBlock - 1) / c_threadsPerBlock;
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    const float3* sendBuf  = d_sendBuf_;
    const float3* d_x      = d_x_;
    const int*    indexMap = d_indexMap_;
    const int     size     = xSendSize_;
    // The coordinateShift changes between steps when we have
    // performed a DD partition, or have updated the box e.g. when
    // performing pressure coupling. So, for simplicity, the box
    // is used every step to pass the shift vector as an argument of
    // the packing kernel.
    const int    boxDimensionIndex = dd_->dim[dimIndex_];
    const float3 coordinateShift{ box[boxDimensionIndex][XX],
                                  box[boxDimensionIndex][YY],
                                  box[boxDimensionIndex][ZZ] };

    // Avoid launching kernel when there is no work to do
    if (size > 0)
    {
        auto kernelFn = usePBC_ ? packSendBufKernel<true> : packSendBufKernel<false>;

        const auto kernelArgs = prepareGpuKernelArguments(
                kernelFn, config, &sendBuf, &d_x, &indexMap, &size, &coordinateShift);

        launchGpuKernel(kernelFn, config, *haloStream_, nullptr, "Domdec GPU Apply X Halo Exchange", kernelArgs);
    }

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuMoveX);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpu);

    // Consider time spent in communicateHaloData as Comm.X counter
    // ToDo: We need further refinement here as communicateHaloData includes launch time for cudamemcpyasync
    wallcycle_start(wcycle_, WallCycleCounter::MoveX);

    // wait for remote co-ordinates is implicit with process-MPI as non-local stream is synchronized before MPI calls
    // and MPI_Waitall call makes sure both neighboring ranks' non-local stream is synchronized before data transfer is initiated
    if (GMX_THREAD_MPI && dimIndex_ == 0 && pulse_ == 0)
    {
        enqueueWaitRemoteCoordinatesReadyEvent(dependencyEvent);
    }

    float3* recvPtr = GMX_THREAD_MPI ? remoteXPtr_ : &d_x_[atomOffset_];
    communicateHaloData(d_sendBuf_, xSendSize_, sendRankX_, recvPtr, xRecvSize_, recvRankX_);

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
    // ToDo: We need further refinement here as communicateHaloData includes launch time for cudamemcpyasync
    wallcycle_start(wcycle_, WallCycleCounter::MoveF);

    while (!dependencyEvents->empty())
    {
        auto* dependency = dependencyEvents->back();
        dependency->enqueueWaitEvent(*haloStream_);
        dependencyEvents->pop_back();
    }

    float3* recvPtr = GMX_THREAD_MPI ? remoteFPtr_ : d_recvBuf_;

    // Communicate halo data
    communicateHaloData(&(d_f_[atomOffset_]), fSendSize_, sendRankF_, recvPtr, fRecvSize_, recvRankF_);

    wallcycle_stop(wcycle_, WallCycleCounter::MoveF);

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpu);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuMoveF);

    float3* d_f = d_f_;

    // Unpack halo buffer into force array

    KernelLaunchConfig config;
    config.blockSize[0]     = c_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = (fRecvSize_ + c_threadsPerBlock - 1) / c_threadsPerBlock;
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    const float3* recvBuf  = d_recvBuf_;
    const int*    indexMap = d_indexMap_;
    const int     size     = fRecvSize_;

    if (pulse_ > 0 || dd_->ndim > 1)
    {
        // We need to accumulate rather than set, since it is possible
        // that, in this pulse/dim, a value could be written to a location
        // corresponding to the halo region of a following pulse/dim.
        accumulateForces = true;
    }

    if (size > 0)
    {
        auto kernelFn = accumulateForces ? unpackRecvBufKernel<true> : unpackRecvBufKernel<false>;

        const auto kernelArgs =
                prepareGpuKernelArguments(kernelFn, config, &d_f, &recvBuf, &indexMap, &size);

        launchGpuKernel(kernelFn, config, *haloStream_, nullptr, "Domdec GPU Apply F Halo Exchange", kernelArgs);
    }

    fReadyOnDevice_.markEvent(*haloStream_);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuMoveF);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpu);
}

void GpuHaloExchange::Impl::communicateHaloData(float3* sendPtr,
                                                int     sendSize,
                                                int     sendRank,
                                                float3* recvPtr,
                                                int     recvSize,
                                                int     recvRank)
{
    if (GMX_THREAD_MPI)
    {
        // no need to explicitly sync with GMX_THREAD_MPI as all operations are
        // anyway launched in correct stream
        communicateHaloDataWithCudaDirect(sendPtr, sendSize, sendRank, recvPtr, recvRank);
    }
    else
    {
        communicateHaloDataWithCudaMPI(sendPtr, sendSize, sendRank, recvPtr, recvSize, recvRank);
    }
}

void GpuHaloExchange::Impl::communicateHaloDataWithCudaMPI(float3* sendPtr,
                                                           int     sendSize,
                                                           int     sendRank,
                                                           float3* recvPtr,
                                                           int     recvSize,
                                                           int     recvRank)
{
    // no need to wait for haloDataReadyOnDevice event if this rank is not sending any data
    if (sendSize > 0)
    {
        // wait for halo stream to complete all outstanding
        // activities, to ensure that buffer is up-to-date in GPU memory
        // before transferring to remote rank

        // ToDo: Replace stream synchronize with event synchronize
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
#endif
}

void GpuHaloExchange::Impl::communicateHaloDataWithCudaDirect(float3* sendPtr,
                                                              int     sendSize,
                                                              int     sendRank,
                                                              float3* remotePtr,
                                                              int     recvRank)
{

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

#if GMX_MPI
    // ensure pushed data has arrived before remote rank progresses
    // This rank records an event and sends it to the remote rank which has just been pushed data.
    // This rank recieves event from remote rank which has pushed data here, and enqueues that event
    // to its stream.
    GpuEventSynchronizer* haloDataTransferRemote;

    GMX_ASSERT(haloDataTransferLaunched_ != nullptr,
               "Halo exchange requires valid event to synchronize data transfer initiated in "
               "remote rank");
    haloDataTransferLaunched_->markEvent(*haloStream_);

    MPI_Sendrecv(&haloDataTransferLaunched_,
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
#else
    GMX_UNUSED_VALUE(sendRank);
    GMX_UNUSED_VALUE(recvRank);
#endif
}

GpuEventSynchronizer* GpuHaloExchange::Impl::getForcesReadyOnDeviceEvent()
{
    return &fReadyOnDevice_;
}

/*! \brief Create Domdec GPU object */
GpuHaloExchange::Impl::Impl(gmx_domdec_t*        dd,
                            int                  dimIndex,
                            MPI_Comm             mpi_comm_mysim,
                            const DeviceContext& deviceContext,
                            int                  pulse,
                            gmx_wallcycle*       wcycle) :
    dd_(dd),
    sendRankX_(dd->neighbor[dimIndex][1]),
    recvRankX_(dd->neighbor[dimIndex][0]),
    sendRankF_(dd->neighbor[dimIndex][0]),
    recvRankF_(dd->neighbor[dimIndex][1]),
    usePBC_(dd->ci[dd->dim[dimIndex]] == 0),
    haloDataTransferLaunched_(GMX_THREAD_MPI ? new GpuEventSynchronizer() : nullptr),
    mpi_comm_mysim_(mpi_comm_mysim),
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
    freeDeviceBuffer(&d_indexMap_);
    freeDeviceBuffer(&d_sendBuf_);
    freeDeviceBuffer(&d_recvBuf_);
    freeDeviceBuffer(&d_fShift_);
    delete haloDataTransferLaunched_;
}

GpuHaloExchange::GpuHaloExchange(gmx_domdec_t*        dd,
                                 int                  dimIndex,
                                 MPI_Comm             mpi_comm_mysim,
                                 const DeviceContext& deviceContext,
                                 int                  pulse,
                                 gmx_wallcycle*       wcycle) :
    impl_(new Impl(dd, dimIndex, mpi_comm_mysim, deviceContext, pulse, wcycle))
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
    impl_->reinitHalo(asFloat3(d_coordinatesBuffer), asFloat3(d_forcesBuffer));
}

GpuEventSynchronizer* GpuHaloExchange::communicateHaloCoordinates(const matrix          box,
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
