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
 * \brief Implements backend-specific part of GPU halo exchange, namely pack and unpack
 * kernels and the host code for scheduling them, using CUDA.
 *
 * NVSHMEM Kernels:
 * Implements NVSHMEM enabled kernel initiated transfers of coordinates and forces.
 *
 * For coordinates, the packSendBufAndPutNvshmemKernel() function packs the coordinates
 * and communicates the packed data to the sendRank. It also signals readiness to receive
 * data from the recvRank using signal counters, ensuring that the recvRank can safely
 * transfer the packed data into this process's buffer.
 *
 * For forces, the putPackedDataUnpackRecvBufNvshmemKernel() function transfers the packed
 * forces to the sendRank and waits for a signal from the recvRank if there is any packed
 * data received from the recvRank needs to be unpacked.
 *
 * \author Alan Gray <alang@nvidia.com>
 * \author Mahesh Doijade <mdoijade@nvidia.com>
 *
 * \ingroup module_domdec
 */
#include "gmxpre.h"

#include "gpuhaloexchange_impl_gpu.h"

#include "config.h"

#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/gpu_utils/vectype_ops_cuda.h"

#include "domdec_struct.h"

#if GMX_NVSHMEM
#    include <nvshmem.h>
#    include <nvshmemx.h>

#    include <cuda/barrier>
#endif

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

#if GMX_NVSHMEM
/*! \brief This kernel is responsible for packing coordinate data into a buffer
 * on the GPU using pre-populated "map" containing index info and then sending this
 * packed data to a remote processing element (PE) or rank using NVSHMEM. It also manages
 * synchronization signals to ensure safe and efficient data transfer.
 *
 * Execution Flow is as below:
 * 1.) Set up signal pointers based on signal offsets which depends on pulse/dim.
 * 2.) If this rank expects to receive data, signal readiness using nvshmemx_signal_op to recvRank.
 * 3.) Loop through assigned indices to pack data into dataPacked.
 * 4.) Use device scoped arrive-wait barrier to ensure all threads have completed their tasks and
       and wait on signal (signalSenderRankX) before proceeding with data transfer.
 * 5.) Perform non-blocking put operations to send packed data to sendRank.
 * 6.) Wait for signal (signalSenderRankX) indicating that all required data has been received
 *     from recvRank.
 * It should be noted that the kernel is launched if either sendSize > 0 or recvSize > 0, as in
 * each of these case we need to synchronize with sendRank or recvRank accordingly.
 *
 * \param[in] data        full array of coordinates
 * \param[out] dataPacked  packed array of coordinates to be transferred
 * \param[in] map         array of indices defining mapping from full to packed array
 * \param[in] sendSize     number of elements in map array to pack and put
 * \param[in] coordinateShift     coordinate shift to be used when usePBC is true
 * \param[in] sendRank    rank at which the packed data needs to be put
 * \param[in] recvRank    rank from which this rank will receive the packed data puts.
 * \param[in] atomOffsetXInSendRank offset in the remote sendRank from where the packed data should be put
 * \param[in] signalReceiverRankX  signal object used for notifying the sendRank completion of puts
 * \param[in] signalSenderRankX signal object used for notifying the recvRank about readiness to receive the puts
 * \param[in] signalXCounter   signal counter value used for notifying put completions for this timestep
 * \param[in] signalObjOffset  Offset to get current pulse/dim signal object
 * \param[in] bar         Device scoped arrive wait barrier used to track all the CTAs data pack completion
 * \param[in] recvSize   number of elements this rank shall receive via puts from recvRank
 */
template<bool usePBC>
__global__ void packSendBufAndPutNvshmemKernel(float3*        data,
                                               float3*        dataPacked,
                                               const int*     map,
                                               const int      sendSize,
                                               const float3   coordinateShift,
                                               const int      sendRank,
                                               const int      recvRank,
                                               const int      atomOffsetXInSendRank,
                                               uint64_t*      signalReceiverRankX,
                                               uint64_t*      signalSenderRankX,
                                               const uint64_t signalXCounter,
                                               const int      signalObjOffset,
                                               cuda::barrier<cuda::thread_scope_device>* bar,
                                               const int                                 recvSize)
{
    using barrier   = cuda::barrier<cuda::thread_scope_device>;
    int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    signalReceiverRankX += signalObjOffset;
    signalSenderRankX += signalObjOffset;

    /* Notify the neighbor about arrival, readiness to receive the packed data via nvshmem_put */
    if (threadIndex == 0 && recvSize > 0)
    {
        nvshmemx_signal_op(signalSenderRankX, signalXCounter, NVSHMEM_SIGNAL_SET, recvRank);
    }

    if (sendSize > 0)
    {
        const int gridSize = blockDim.x * gridDim.x;

        for (int idx = threadIndex; idx < sendSize; idx += gridSize)
        {
            int           atomIndex   = map[idx];
            const float3* gm_dataSrc  = &data[atomIndex];
            float3*       gm_dataDest = &dataPacked[idx];
            if (usePBC)
            {
                *gm_dataDest = *gm_dataSrc + coordinateShift;
            }
            else
            {
                *gm_dataDest = *gm_dataSrc;
            }
        }
        __syncthreads();

        // More study is needed to understand why more than 1 block nvshmem_put performs bad
        constexpr int          numFinalBlocks = 1;
        barrier::arrival_token token;
        if (threadIdx.x == 0)
        {
            token = bar->arrive();
        }

        if (blockIdx.x < numFinalBlocks)
        {
            const int chunkSize   = sendSize / numFinalBlocks;
            const int blockOffset = blockIdx.x * chunkSize;
            const int dataSize    = blockIdx.x != (numFinalBlocks - 1)
                                            ? chunkSize
                                            : max(sendSize - blockOffset, chunkSize);

            float* recvPtr = reinterpret_cast<float*>(&data[atomOffsetXInSendRank + blockOffset]);

            if (threadIdx.x == 0)
            {
                bar->wait(std::move(token));
                nvshmem_signal_wait_until(signalSenderRankX, NVSHMEM_CMP_EQ, signalXCounter);
            }
            __syncthreads();

            nvshmemx_float_put_signal_nbi_block(recvPtr,
                                                reinterpret_cast<const float*>(&dataPacked[blockOffset]),
                                                dataSize * 3,
                                                signalReceiverRankX,
                                                signalXCounter,
                                                NVSHMEM_SIGNAL_SET,
                                                sendRank);
        }
    }

    if (threadIndex == 0 && recvSize > 0)
    {
        // We need to wait for the signal completion before exiting the kernel
        // in order to make sure the packed data received from recvRank's is complete.
        nvshmem_signal_wait_until(signalReceiverRankX, NVSHMEM_CMP_EQ, signalXCounter);
    }
}
#endif


/*! \brief This kernel is responsible for transferring packed force data to a remote
 * processing element(PE) or rank and unpacking received non-local forces data, using NVSHMEM.
 * It also handles synchronization through signaling ensuring safe data transfers.
 *
 * Execution Flow is as below:
 * 1.) Set up signal pointers based on signal offsets which depends on pulse/dim.
 * 2.) If sendSize is greater than zero, perform non-blocking put operations to send packed
 *     data to sendRank.
 * 3.) If this rank expects to receive data, wait for the signal indicating that all required
 *     data has been received.
 * 4.) And then loop through assigned indices to unpack the received dataPacked into data.
 * It should be noted that the kernel is launched if either sendSize > 0 or recvSize > 0, as in
 * each of these case we need to synchronize with sendRank or recvRank accordingly.
 *
 * \param[out] data        full array of force values
 * \param[in]  dataPacked  packed array of force values to be transferred
 * \param[in]  map         array of indices defining mapping from full to packed array
 * \param[in]  recvSize    number of elements received to unpack
 * \param[in]  sendSize    number of elements of the packed data to send
 * \param[in]  atomOffset  Offset from where the packed data starts
 * \param[in]  sendRank    PP rank where the packed data is supposed to be put
 * \param[in]  signalReceiverRankF  signal object used for notifying put PP rank
 * \param[in]  signalReceiverRankFCounter signal counter value used for notifying put completion to
 * PP rank \param[in]  signalObjOffset  Offset to get current pulse/dim signal object
 */
template<bool accumulate>
__global__ void putPackedDataUnpackRecvBufNvshmemKernel(float3* __restrict__ data,
                                                        float3* __restrict__ dataPacked,
                                                        const int* __restrict__ map,
                                                        const int recvSize,
                                                        const int sendSize,
                                                        const int atomOffset,
                                                        const int sendRank,
                                                        uint64_t* signalReceiverRankF,
                                                        uint64_t  signalReceiverRankFCounter,
                                                        int       signalObjOffset)
{
#if GMX_NVSHMEM
    // put approach
    signalReceiverRankF = signalReceiverRankF + signalObjOffset;
    if (blockIdx.x == 0 && sendSize > 0)
    {
        nvshmemx_float_put_signal_nbi_block(reinterpret_cast<float*>(dataPacked),
                                            reinterpret_cast<float*>(&data[atomOffset]),
                                            sendSize,
                                            signalReceiverRankF,
                                            signalReceiverRankFCounter,
                                            NVSHMEM_SIGNAL_SET,
                                            sendRank);
    }

    if (recvSize == 0)
    {
        // when recvSize is 0 there is no packed data received by this rank to unpack
        // hence we exit.
        return;
    }

    const int     threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    const float3* gm_dataSrc;
    float3*       gm_dataDest;
    if (threadIndex < recvSize)
    {
        gm_dataSrc  = &dataPacked[threadIndex];
        gm_dataDest = &data[map[threadIndex]];
    }

    if (threadIdx.x == 0)
    {
        nvshmem_signal_wait_until(signalReceiverRankF, NVSHMEM_CMP_EQ, signalReceiverRankFCounter);
    }
    __syncthreads();

    if (threadIndex < recvSize)
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
#endif
}

void GpuHaloExchange::Impl::launchPackXKernel(const matrix box)
{
    // launch kernel to pack send buffer
    KernelLaunchConfig config;
    config.blockSize[0]     = c_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = gmx::divideRoundUp(xSendSize_, c_threadsPerBlock);
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    const float3* sendBuf  = asFloat3(d_sendBuf_);
    const float3* d_x      = asFloat3(d_x_);
    const int*    indexMap = d_indexMap_;
    // The coordinateShift changes between steps when we have
    // performed a DD partition, or have updated the box e.g. when
    // performing pressure coupling. So, for simplicity, the box
    // is used every step to pass the shift vector as an argument of
    // the packing kernel.
    const int    boxDimensionIndex = dd_->dim[dimIndex_];
    const float3 coordinateShift{ box[boxDimensionIndex][XX],
                                  box[boxDimensionIndex][YY],
                                  box[boxDimensionIndex][ZZ] };
    if (useNvshmem_)
    {
#if GMX_NVSHMEM
        // Launch the coordinates pack & put kernel only if either we have some buffer to pack (xSendSize_ > 0)
        // and put it to sendRankX_ or we have some data received via puts (xRecvSize_ > 0) which needs to be waited on.
        if ((xSendSize_ > 0) || (xRecvSize_ > 0))
        {
            config.blockSize[0] = nvshmemHaloExchange_.c_nvshmemThreadsPerBlock;
            config.gridSize[0]  = nvshmemHaloExchange_.gridDimX_;

            auto       packNvShmemkernelFn = usePBC_ ? packSendBufAndPutNvshmemKernel<true>
                                                     : packSendBufAndPutNvshmemKernel<false>;
            const auto kernelArgs =
                    prepareGpuKernelArguments(packNvShmemkernelFn,
                                              config,
                                              &d_x,
                                              &sendBuf,
                                              &indexMap,
                                              &xSendSize_,
                                              &coordinateShift,
                                              &sendRankX_,
                                              &recvRankX_,
                                              &nvshmemHaloExchange_.putAtomOffsetInReceiverRankXBuf_,
                                              &nvshmemHaloExchange_.d_signalReceiverRankX_,
                                              &nvshmemHaloExchange_.d_signalSenderRankX_,
                                              &nvshmemHaloExchange_.signalReceiverRankXCounter_,
                                              &nvshmemHaloExchange_.signalObjOffset_,
                                              &nvshmemHaloExchange_.d_arriveWaitBarrier_,
                                              &xRecvSize_);
            launchGpuKernel(packNvShmemkernelFn,
                            config,
                            *haloStream_,
                            nullptr,
                            "Domdec GPU Apply X Halo Exchange",
                            kernelArgs);
        }
        nvshmemHaloExchange_.signalReceiverRankXCounter_++;
#endif
    }
    else
    {
        // Avoid launching kernel when there is no work to do
        if (xSendSize_ > 0)
        {
            auto       kernelFn   = usePBC_ ? packSendBufKernel<true> : packSendBufKernel<false>;
            const auto kernelArgs = prepareGpuKernelArguments(
                    kernelFn, config, &sendBuf, &d_x, &indexMap, &xSendSize_, &coordinateShift);
            launchGpuKernel(
                    kernelFn, config, *haloStream_, nullptr, "Domdec GPU Apply X Halo Exchange", kernelArgs);
        }
    }
}

// The following method should be called after non-local buffer operations,
// and before the local buffer operations.
void GpuHaloExchange::Impl::launchUnpackFKernel(bool accumulateForces)
{
    KernelLaunchConfig config;
    config.blockSize[0]     = c_threadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = divideRoundUp(fRecvSize_, c_threadsPerBlock);
    config.gridSize[1]      = 1;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    float3*       d_f      = asFloat3(d_f_);
    const float3* recvBuf  = asFloat3(d_recvBuf_);
    const int*    indexMap = d_indexMap_;

    if (useNvshmem_ && receiveInPlace_)
    {
        // Launch the force put & unpack kernel only if either we have some packed buffer to put
        // (fSendSize_ > 0) or we have some data received (fRecvSize_ > 0) to unpack.
        if ((fSendSize_ > 0) || (fRecvSize_ > 0))
        {
            config.blockSize[0] = nvshmemHaloExchange_.c_nvshmemThreadsPerBlock;
            // Add 1 to fRecvSize_ to take care of case when fRecvSize_ == 0 in such case
            // we still launch the kernel to signal the remote rank recvRankX_
            // from which this rank receives packed data as (fSendSize_ > 0)
            config.gridSize[0] =
                    gmx::divideRoundUp(fRecvSize_ + 1, nvshmemHaloExchange_.c_nvshmemThreadsPerBlock);
            const int  fSendSize = fSendSize_ * DIM;
            auto       kernelFn  = accumulateForces ? putPackedDataUnpackRecvBufNvshmemKernel<true>
                                                    : putPackedDataUnpackRecvBufNvshmemKernel<false>;
            const auto kernelArgs =
                    prepareGpuKernelArguments(kernelFn,
                                              config,
                                              &d_f,
                                              &recvBuf,
                                              &indexMap,
                                              &fRecvSize_,
                                              &fSendSize,
                                              &atomOffset_,
                                              &sendRankF_,
                                              &nvshmemHaloExchange_.d_signalReceiverRankF_,
                                              &nvshmemHaloExchange_.signalReceiverRankFCounter_,
                                              &nvshmemHaloExchange_.signalObjOffset_);
            launchGpuKernel(
                    kernelFn, config, *haloStream_, nullptr, "Domdec GPU Apply F Halo Exchange", kernelArgs);
        }
        nvshmemHaloExchange_.signalReceiverRankFCounter_++;
    }
    else
    {
        if (fRecvSize_ > 0)
        {
            auto kernelFn = accumulateForces ? unpackRecvBufKernel<true> : unpackRecvBufKernel<false>;

            const auto kernelArgs =
                    prepareGpuKernelArguments(kernelFn, config, &d_f, &recvBuf, &indexMap, &fRecvSize_);

            launchGpuKernel(
                    kernelFn, config, *haloStream_, nullptr, "Domdec GPU Apply F Halo Exchange", kernelArgs);
        }
    }
}

} // namespace gmx
