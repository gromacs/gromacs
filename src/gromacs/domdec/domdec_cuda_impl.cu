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
 * \brief Implements Domdec halo exchange using CUDA
 *
 *
 * \author Alan Gray <alang@nvidia.com.com>
 *
 * \ingroup module_domdec
 */
#include "gmxpre.h"

#include "domdec_cuda_impl.h"

#include "config.h"

#include <assert.h>
#include <stdio.h>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_cuda.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/pbcutil/ishift.h"

#include "domdec_internal.h"

namespace gmx
{

//! Number of CUDA threads in a block
//TODO Optimize this through experimentation
constexpr static int c_threadsPerBlock = 256;

/*! \brief pack non-local coordinate data buffer on the GPU using pre-populated "map" containing index information
 * \param[out] dataPacked  packed array of position values to be transferred
 * \param[in] data         full array of position values
 * \param[in] map          array of indices defining mapping from full to packed array
 * \param[in] mapSize      number of elements in map array
 * \param[in] xShift       shift buffer
 */
template <bool usePBC>
__global__ void packSendBufKernel(rvec                      *dataPacked,
                                  const rvec * __restrict__  data,
                                  const int * __restrict__   map,
                                  const int                  mapSize,
                                  const float3               xShift);

template <bool usePBC>
__global__ void packSendBufKernel(rvec                      *dataPacked,
                                  const rvec * __restrict__  data,
                                  const int * __restrict__   map,
                                  const int                  mapSize,
                                  const float3               xShift)
{
    int     threadIndex       = blockIdx.x*blockDim.x+threadIdx.x;
    float3 *dataDest          = (float3*) dataPacked[threadIndex];
    float3 *dataSrc           = (float3*) data[map[threadIndex]];

    if (threadIndex < mapSize)
    {
        if (usePBC)
        {
            *dataDest = *dataSrc + xShift;
        }
        else
        {
            *dataDest = *dataSrc;
        }

    }

    return;
}


/*! \brief unpack non-local force data buffer on the GPU using pre-populated "map" containing index information
 * \param[out] data        full array of force values
 * \param[out] fshift      force shift buffer
 * \param[in]  dataPacked  packed array of force values to be transferred
 * \param[in]  map         array of indices defining mapping from full to packed array
 * \param[in]  mapSize     number of elements in map array
 */
template <bool accumulateShiftForces>
__global__ void unpackRecvBufKernel(rvec                *data,
                                    float3              *fshift,
                                    const rvec          *dataPacked,
                                    const int           *map,
                                    const int            mapSize);

template <bool accumulateShiftForces>
__global__ void unpackRecvBufKernel(rvec                *data,
                                    float3              *fshift,
                                    const rvec          *dataPacked,
                                    const int           *map,
                                    const int            mapSize)
{

    int     threadIndex        = blockIdx.x*blockDim.x+threadIdx.x;
    float3 *dataSrc            = (float3*) dataPacked[threadIndex];
    float3 *dataDest           = (float3*) data[map[threadIndex]];

    if (threadIndex < mapSize)
    {
        *dataDest += *dataSrc;
    }

    if (accumulateShiftForces)
    {
        //Atomically add source data float3 element to force shift buffer.
        //Stage through shared memory for performance
        __shared__ float3 tmpFshift;
        tmpFshift = make_float3(0.0f);
        __syncthreads();
        //each thread in this block updates shared memory float3
        atomicAdd(&tmpFshift, *dataSrc);
        __syncthreads();
        //first thread in block adds shared memory float3 to global float3
        if (threadIdx.x == 0)
        {
            atomicAdd(fshift, tmpFshift);
        }
    }

    return;
}

/*! \brief
 * This function is based on move_x, but instead of actually sending and recieving data
 * a map structure is created containing information on indexes to be sent/recieved.  This
 * map will later be passed to the GPU and used to pack a buffer, which will be sent
 * directly between GPUs.
 *
 */
void DomdecCuda::Impl::initHaloExchange(gmx_domdec_t *dd, matrix box, rvec *d_x_ptr, void *stream_ptr)
{

    cudaStream_t                  stream            = (cudaStream_t) stream_ptr;
    int                           nzone             = 1;
    gmx_domdec_comm_t            *comm              = dd->comm;;
    gmx_domdec_comm_dim_t        *cd                = &comm->cd[0];
    const gmx_domdec_ind_t       &ind               = cd->ind[0];
    int                           newSize           = ind.nsend[nzone+1];

    if (newSize > xSendSize_) //(re-)allocate required structures
    {
        h_indexMap_.resize(newSize);
        h_sendBuf_.resize(newSize);
        h_recvBuf_.resize(newSize);
        reallocateDeviceBuffer(&d_indexMap_, newSize, &indexMap_size_, &indexMap_size_alloc_, nullptr);
        reallocateDeviceBuffer(&d_sendBuf_, newSize, &sendBuf_size_, &sendBuf_size_alloc_, nullptr);
        reallocateDeviceBuffer(&d_recvBuf_, newSize, &recvBuf_size_, &recvBuf_size_alloc_, nullptr);
    }

    xSendSize_ = newSize;
    MPI_Sendrecv(&xSendSize_, sizeof(int), MPI_BYTE, sendRankX_, 0,
                 &xRecvSize_, sizeof(int), MPI_BYTE, recvRankX_, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    fSendSize_ = xRecvSize_;
    fRecvSize_ = xSendSize_;

    localOffset_ = comm->atomRanges.numHomeAtoms();  //offset for data recieved by this rank

    int n = 0;
    for (int j : ind.index)
    {
        h_indexMap_[n] = j; //add index to index map
        n++;
    }
    copyToDeviceBuffer(&d_indexMap_, h_indexMap_.data(), 0, newSize, stream, GpuApiCallBehavior::Async, nullptr);

    //three components of shift for PBC
    xShift_ = make_float3(box[dd->dim[0]][0],
                          box[dd->dim[0]][1],
                          box[dd->dim[0]][2]);

    bShiftForcesNeedPbc_   = (dd->ci[0] == 0);

    if (GMX_THREAD_MPI)
    {
        // This rank will push data to its neighbor, so needs to know
        // the remote receive address and similarly send its receive
        // address to other neighbour. We can do this here in init fn
        // since the pointers will not change untul the next NS step.

        //Position buffer:
        void* recvPtr  = static_cast<void*> (&d_x_ptr[localOffset_][0]);
        MPI_Sendrecv(&recvPtr, sizeof(void*), MPI_BYTE, recvRankX_, 0,
                     &remoteXPtr_, sizeof(void*), MPI_BYTE, sendRankX_, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //Force buffer:
        recvPtr  = static_cast<void*> (d_recvBuf_);
        MPI_Sendrecv(&recvPtr, sizeof(void*), MPI_BYTE, recvRankF_, 0,
                     &remoteFPtr_, sizeof(void*), MPI_BYTE, sendRankF_, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return;
}

void DomdecCuda::Impl::applyXHaloExchange(rvec* d_x_ptr, void *streamNonLocal_ptr)
{
    cudaStream_t streamNonLocal = (cudaStream_t) streamNonLocal_ptr;

    // launch kernel to pack send buffer

    KernelLaunchConfig config;
    config.blockSize[0]     = c_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = (xSendSize_+c_threadsPerBlock-1)/c_threadsPerBlock;
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;
    config.stream           = streamNonLocal;

    rvec             *sendBuf  = d_sendBuf_;
    const rvec       *d_x      =  d_x_ptr;
    const int        *indexMap = d_indexMap_;
    const int         size     = xSendSize_;

    auto              kernelFn = usePBC_ ? packSendBufKernel<true> : packSendBufKernel<false>;

    const auto        kernelArgs   = prepareGpuKernelArguments(kernelFn, config, &sendBuf, &d_x, &indexMap,
                                                               &size, &xShift_);

    launchGpuKernel(kernelFn, config, nullptr, "Domdec GPU Apply X Halo Exchange", kernelArgs);

    haloDataTransfer(d_x_ptr, streamNonLocal, HaloXorF::HaloX);

    return;
}

void DomdecCuda::Impl::applyFHaloExchange(rvec* d_f_ptr, rvec *fshift, void *stream_ptr)
{

    haloDataTransfer(d_f_ptr, stream_ptr, HaloXorF::HaloF);

    cudaStream_t       stream = (cudaStream_t) stream_ptr;

    KernelLaunchConfig config;
    config.blockSize[0]     = c_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = (fRecvSize_+c_threadsPerBlock-1)/c_threadsPerBlock;
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;
    config.stream           = stream;

    rvec            *d_f        = d_f_ptr;
    const rvec      *recvBuf    = d_recvBuf_;
    const int       *indexMap   = d_indexMap_;
    const int        size       = fRecvSize_;

    auto             kernelFn = bShiftForcesNeedPbc_ ? unpackRecvBufKernel<true> : unpackRecvBufKernel<false>;


    ivec    vis;
    clear_ivec(vis);
    vis[0] = 1;
    int    is          = IVEC2IS(vis);
    rvec** d_fShiftPtr = static_cast<rvec**> (static_cast<void*> (&d_fShift_));

    if (bShiftForcesNeedPbc_)     // update shift force buffer if required
    {
        copyToDeviceBuffer(d_fShiftPtr, &fshift[is], 0, 1, stream, GpuApiCallBehavior::Async, nullptr);
    }

    const auto       kernelArgs   = prepareGpuKernelArguments(kernelFn, config, &d_f, &d_fShift_,
                                                              &recvBuf, &indexMap,
                                                              &size);

    launchGpuKernel(kernelFn, config, nullptr, "Domdec GPU Apply F Halo Exchange", kernelArgs);

    if (bShiftForcesNeedPbc_)
    {
        copyFromDeviceBuffer(&fshift[is], d_fShiftPtr, 0, 1, stream, GpuApiCallBehavior::Async, nullptr);
    }
}


void DomdecCuda::Impl::haloDataTransfer(rvec* d_ptr, void *streamNonLocal_ptr, HaloXorF haloXorF)
{

    void                        * sendPtr;
    void                        * recvPtr;
    int                           sendSize;
    int                gmx_unused recvSize;
    void * remotePtr;
    int    sendRank;
    int    recvRank;
    if (haloXorF == HaloXorF::HaloX)
    {
        sendPtr   = static_cast<void*> (d_sendBuf_);
        recvPtr   = static_cast<void*> (&d_ptr[localOffset_][0]);
        sendSize  = xSendSize_;
        recvSize  = xRecvSize_;
        remotePtr = remoteXPtr_;
        sendRank  = sendRankX_;
        recvRank  = recvRankX_;
    }
    else
    {
        sendPtr   = static_cast<void*> (&(d_ptr[localOffset_][0]));
        recvPtr   = static_cast<void*> (d_recvBuf_);
        sendSize  = fSendSize_;
        recvSize  = fRecvSize_;
        remotePtr = remoteFPtr_;
        sendRank  = sendRankF_;
        recvRank  = recvRankF_;
    }

    if (GMX_THREAD_MPI)
    {
        haloDataTransferCudaDirect(sendPtr, sendSize, sendRank, remotePtr, recvRank, streamNonLocal_ptr);
    }
    else
    {
        haloDataTransferCudaMPI(sendPtr, sendSize, sendRank, recvPtr, recvSize, recvRank, streamNonLocal_ptr);
    }
}

void DomdecCuda::Impl::haloDataTransferCudaMPI(void *sendPtr, int sendSize, int sendRank,
                                               void *recvPtr, int recvSize, int recvRank,
                                               void *streamNonLocal_ptr)
{

    cudaStream_t streamNonLocal = (cudaStream_t) streamNonLocal_ptr;

    // wait for non local stream to complete all outstanding
    // activities, to ensure that buffer is up-to-date in GPU memory
    // before transferring to remote rank
    GpuEventSynchronizer syncStreamForHalo;
    syncStreamForHalo.markEvent(streamNonLocal);
    syncStreamForHalo.waitForEvent();

    // perform halo exchange directly in device buffers

    MPI_Request request[2];
    MPI_Status  status[2];

    // recv remote data into halo region
    MPI_Irecv(recvPtr, recvSize*3, MPI_FLOAT, recvRank, 0,
              MPI_COMM_WORLD, &request[1]);

    // send data to remote halo region
    MPI_Isend(sendPtr, sendSize*3, MPI_FLOAT, sendRank, 0,
              MPI_COMM_WORLD, &request[0]);

    MPI_Waitall(2, request, status);

}


void DomdecCuda::Impl::haloDataTransferCudaDirect(void *sendPtr, int sendSize, int sendRank,
                                                  void* remotePtr, int recvRank, void *streamNonLocal_ptr)
{

    cudaError_t  stat;
    cudaStream_t streamNonLocal = (cudaStream_t) streamNonLocal_ptr;

    // We asynchronously push data to remote rank. The remote
    // destination pointer has already been set in the init fn.  We
    // don't need to worry about overwriting data the remote ranks
    // still needs since the halo exchange is just done once per
    // timestep, for each of X and F.

    //push data to neighbor
    stat = cudaMemcpyAsync(remotePtr, sendPtr, sendSize*3*sizeof(float), cudaMemcpyDeviceToDevice, streamNonLocal);
    CU_RET_ERR(stat, "cudaMemcpyAsync on GPU Domdec CUDA direct data transfer failed");

    //ensure pushed data has arrived before remote rank progresses
    // This rank records an event and sends it to the remote rank which has just been pushed data.
    // This rank recieves event from remote rank which has pushed data here, and enqueues that event to
    // its stream.
    GpuEventSynchronizer *haloDataTransferRemote;

    haloDataTransferLaunched_->markEvent(streamNonLocal);

    MPI_Sendrecv(&haloDataTransferLaunched_, sizeof(GpuEventSynchronizer*), MPI_BYTE, sendRank, 0,
                 &haloDataTransferRemote, sizeof(GpuEventSynchronizer*), MPI_BYTE, recvRank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    haloDataTransferRemote->enqueueWaitEvent(streamNonLocal);

}

void DomdecCuda::Impl::setupCudaDevicePeerAccess()
{

    cudaError_t stat;

    // take a note of currently-set GPU
    int gpuRank0;
    stat = cudaGetDevice(&gpuRank0);
    CU_RET_ERR(stat, "cudaGetDevice in GPU Domdec init failed");

    int deviceCount;
    stat = cudaGetDeviceCount(&deviceCount);
    CU_RET_ERR(stat, "cudaGetDeviceCount in GPU Domdec init failed");

    for (int gpuA = 0; gpuA < deviceCount; gpuA++)
    {
        stat = cudaSetDevice(gpuA);
        CU_RET_ERR(stat, "cudaSetDevice in GPU Domdec init failed");
        for (int gpuB = gpuA+1; gpuB < deviceCount; gpuB++)
        {
            int canAccessPeer = 0;
            stat = cudaDeviceCanAccessPeer (&canAccessPeer, gpuA, gpuB);
            CU_RET_ERR(stat, "cudaDeviceCanAccessPeer in GPU Domdec init failed");
            if (canAccessPeer)
            {
                stat = cudaDeviceEnablePeerAccess(gpuB, 0);
                CU_RET_ERR(stat, "cudaDeviceEnablePeerAccess in GPU Domdec init failed");
            }
        }
    }

    //re-set GPU to that originally set
    stat = cudaSetDevice(gpuRank0);
    CU_RET_ERR(stat, "cudaSetDevice in GPU Domdec init failed");
}


/*! \brief Create Domdec GPU object */
DomdecCuda::Impl::Impl(gmx_domdec_t *dd)
{

    int                           nzone          = 1;
    gmx_domdec_comm_t            *comm           = dd->comm;;
    gmx_domdec_comm_dim_t        *cd             = &comm->cd[0];
    bool                          usePBC         = (dd->ci[dd->dim[0]] == 0);
    const gmx_domdec_ind_t       &ind            = cd->ind[0];
    int                           size           = ind.nsend[nzone+1];

    if (dd->ndim > 1)
    {
        gmx_fatal(FARGS, "Error: dd->ndim > 1 is not yet supported in domdec_gpu");
    }

    if (!cd->receiveInPlace)
    {
        gmx_fatal(FARGS, "Error: out-of-place recieve is not yet supported in domdec_gpu");
    }

    if (usePBC && dd->bScrewPBC && dd->dim[0] == XX)
    {
        gmx_fatal(FARGS, "Error: screw is not yet supported in domdec_gpu\n");
    }

    changePinningPolicy(&h_indexMap_, gmx::PinningPolicy::PinnedIfSupported);
    changePinningPolicy(&h_sendBuf_, gmx::PinningPolicy::PinnedIfSupported);
    changePinningPolicy(&h_recvBuf_, gmx::PinningPolicy::PinnedIfSupported);
    h_indexMap_.resize(size);
    h_sendBuf_.resize(size);
    h_recvBuf_.resize(size);

    allocateDeviceBuffer(&d_indexMap_, size, nullptr);
    allocateDeviceBuffer(&d_sendBuf_, size, nullptr);
    allocateDeviceBuffer(&d_recvBuf_, size, nullptr);

    allocateDeviceBuffer(&d_fShift_, 1, nullptr);

    sendRankX_       = dd->neighbor[0][1]; //rank to send data to for X
    recvRankX_       = dd->neighbor[0][0]; //rank to recv data from for X
    sendRankF_       = dd->neighbor[0][0]; //rank to send data to for F
    recvRankF_       = dd->neighbor[0][1]; //rank to recv data from for F
    usePBC_          = usePBC;             // Period Boundary Conditions for this rank

    haloDataTransferLaunched_ = new GpuEventSynchronizer();

    if (GMX_THREAD_MPI && (dd->rank == 0))
    {
        // rank 0 enables peer access between all pairs of GPUs where available
        setupCudaDevicePeerAccess();
    }
}

DomdecCuda::Impl::~Impl()
{
    freeDeviceBuffer(&d_indexMap_);
    freeDeviceBuffer(&d_sendBuf_);
    freeDeviceBuffer(&d_recvBuf_);
    freeDeviceBuffer(&d_fShift_);
    delete haloDataTransferLaunched_;
}

DomdecCuda::DomdecCuda(gmx_domdec_t *dd)
    : impl_(new Impl(dd))
{
}

DomdecCuda::~DomdecCuda() = default;

void DomdecCuda::initHaloExchange(gmx_domdec_t *dd, matrix box, rvec *d_x_ptr, void *stream_ptr)
{
    impl_->initHaloExchange(dd, box, d_x_ptr, stream_ptr);
}

void DomdecCuda::applyXHaloExchange(rvec* d_x_ptr, void *streamNonLocal_ptr)
{
    impl_->applyXHaloExchange(d_x_ptr, streamNonLocal_ptr);
}

void DomdecCuda::applyFHaloExchange(rvec* d_f_ptr, rvec *fshift, void *streamNonLocal_ptr)
{
    impl_->applyFHaloExchange(d_f_ptr, fshift, streamNonLocal_ptr);
}

} //namespace gmx
