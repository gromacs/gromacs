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
 * For coordinates, the fusedPulsesPackXAndSendKernel() function packs the coordinates
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
#include "fused_gpuhaloexchange.h"
#if GMX_NVSHMEM
#    include <nvshmem.h>
#    include <nvshmemx.h>
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

/*!
 * \brief GPU-scoped acquire load of a 32-bit value from global memory.
 *
 * PTX: ld.acquire.gpu.global.u32
 *
 * Semantics:
 * - Acquire ordering at GPU scope: subsequent global reads/writes by this thread
 *   cannot be reordered before this load. Ensures visibility of prior writes by
 *   other thread blocks on the same GPU (e.g., grid-synchronization counters)
 *   before consuming dependent data.
 *
 * \param[in] ptr  Address of the 32-bit value in global memory
 * \return The loaded 32-bit value
 */
inline __device__ uint32_t loadAcquireGpuAsm(const uint32_t* ptr)
{
    uint32_t retval = 0;
#if __CUDA_ARCH__ >= 700
    asm("ld.acquire.gpu.global.u32 %0, [%1];" : "=r"(retval) : "l"(ptr) : "memory");
#endif
    return retval;
}

/*!
 * \brief System-scoped acquire load of a 64-bit value from global memory.
 *  we need system scoped ld.acquire as the signal reads though local to
 *  the gpu but are written by remote gpu.
 * PTX: ld.acquire.sys.global.u64
 *
 * Semantics:
 * - Acquire ordering at system scope: subsequent global reads/writes by this thread
 *   cannot be reordered before this load. Ensures visibility of remote GPU writes
 *   (e.g., NVSHMEM signals) before consuming dependent data.
 */
inline __device__ uint64_t loadAcquireSysAsm(const uint64_t* ptr)
{
    uint64_t retval = 0;
#if __CUDA_ARCH__ >= 700
    asm("ld.acquire.sys.global.u64 %0, [%1];" : "=l"(retval) : "l"(ptr) : "memory");
#endif
    return retval;
}

/*!
 * \brief System-scoped relaxed load of a 64-bit value from global memory.
 *
 * PTX: ld.relaxed.sys.global.u64
 *  Relaxed loads at system scope
 */
inline __device__ uint64_t loadRelaxedSysAsm(const uint64_t* ptr)
{
    uint64_t retval = 0;
#if __CUDA_ARCH__ >= 700
    asm("ld.relaxed.sys.global.u64 %0, [%1];" : "=l"(retval) : "l"(ptr) : "memory");
#endif
    return retval;
}

/*!
 * \brief System-scoped release store of a 64-bit value to global memory.
 *
 * PTX: st.release.sys.global.u64
 *
 * Semantics:
 * - Release ordering at system scope: all prior global writes by this thread
 *   become visible to system peers before the store is observable. Used to
 *   publish completion signals after making packed data globally visible.
 */
inline __device__ void storeReleaseSysAsm(uint64_t* ptr, const uint64_t val)
{
#if __CUDA_ARCH__ >= 700
    asm("st.release.sys.global.u64 [%0], %1;" : : "l"(ptr), "l"(val) : "memory");
#endif
}

/*!
 * \brief System-scoped relaxed store of a 64-bit value to global memory.
 *
 * PTX: st.relaxed.sys.global.u64
 *
 * Semantics:
 * - Relaxed ordering at system scope: does not impose ordering on prior global writes
 *   from this thread; only this store is guaranteed visible at system scope.
 * - Appropriate when no earlier data from this thread needs to be ordered with the signal.
 */
inline __device__ void stRelaxedSysAsm(uint64_t* ptr, const uint64_t val)
{
#if __CUDA_ARCH__ >= 700
    asm("st.relaxed.sys.global.u64 [%0], %1;" : : "l"(ptr), "l"(val) : "memory");
#endif
}


/*!
 * \brief Atomically increment counter with release semantics.
 *
 * This wraps the PTX instruction:
 *   atom.inc.release.gpu.global.u32 old, [addr], modMinusOne
 *
 * - The counter at addr is incremented modulo (modMinusOne + 1) and the OLD value is returned.
 * - The .release qualifier ensures that all prior global writes by this thread are made visible
 *   at GPU scope before the atomic is observed by others. This acts as a cache-flushing
 *   fence for data the threadBlock produced before signalling its arrival.
 * - We use modulo (numBlocks - 1) + 1 so that the last arriving block can be detected by
 *   comparing the returned old value against (numBlocks - 1).
 */
inline __device__ uint32_t atomicIncReleaseGpu(uint32_t* addr, int32_t modMinusOne)
{
    uint32_t old = 0;
#if __CUDA_ARCH__ >= 700
    asm("atom.inc.release.gpu.global.u32 %0,[%1],%2;"
        : "=r"(old)
        : "l"(addr), "r"(modMinusOne)
        : "memory");
#endif
    return old;
}

/*!
 * \brief Pack a subset of indices into the destination buffer.
 *
 * Packs elements selected by \p map from \p data into \p gm_dataDest using a
 * grid-stride loop. Depending on the template parameters, either all elements
 * are packed unconditionally or only a conditional subset is packed based on
 * a threshold comparison.
 *
 * Template parameters:
 * - packAll: When true, packs all indices without any threshold checks.
 * - packLessThan: When packAll is false, selects elements using
 *   (atomIndex < threshold) if true, or (atomIndex >= threshold) if false.
 *
 * \param[in,out] gm_dataDest        Destination buffer for packed values
 * \param[in,out] data               Full source array of float3
 * \param[in]     map                Index map selecting elements to pack
 * \param[in]     sendSize           Number of elements to consider
 * \param[in]     gridStride         Stride for grid-stride loop
 * \param[in]     threadIndex        Linearized thread index
 * \param[in]     usePBC             Whether to apply periodic shift
 * \param[in]     coordinateShift    Shift to apply when usePBC is true
 * \param[in]     threshold          Threshold atom index used for selection
 * \param[in,out] hasDependencyAtoms Incremented by the number of elements that
 *                                   did not match the predicate in the first pass
 *                                   used only when packLessThan is true.
 */
template<bool packAll, bool packLessThan>
inline __device__ void packSubset(float3*       gm_dataDest,
                                  float3*       data,
                                  const int*    map,
                                  int           sendSize,
                                  int           gridStride,
                                  int           threadIndex,
                                  bool          usePBC,
                                  const float3& coordinateShift,
                                  int           threshold,
                                  int&          hasDependencyAtoms)
{
    for (int idx = threadIndex; idx < sendSize; idx += gridStride)
    {
        int    atomIndex = map[idx];
        float3 srcVal    = data[atomIndex];
        if constexpr (packAll)
        {
            gm_dataDest[idx] = usePBC ? srcVal + coordinateShift : srcVal;
        }
        else
        {
            bool doPackNow = packLessThan ? (atomIndex < threshold) : (atomIndex >= threshold);
            if (doPackNow)
            {
                gm_dataDest[idx] = usePBC ? srcVal + coordinateShift : srcVal;
            }
            else if (packLessThan)
            {
                // Track presence of deferred atoms in the first pass
                hasDependencyAtoms++;
            }
        }
    }
}

/*!
 * \brief Pack X coordinates into a contiguous buffer (non-TMA).
 *
 * Packs selected atoms from the full coordinate array \p data into the destination
 * buffer \p gm_dataDest using the index \p map. When \p usePBC is true, applies the
 * provided \p coordinateShift to each packed value. For pulses beyond the first,
 * the function enforces ordering between pulses by:
 * - Packing all atoms whose indices are below \p dependencyAtomOffset immediately.
 * - Waiting on completion signals from earlier pulses using \p signalReceiverRankXCurr
 *   and \p signalXCounter.
 * - Packing remaining atoms (indices >= \p dependencyAtomOffset) only after the
 *   dependency condition is met across the grid.
 *
 * Notes:
 * - This routine performs only packing and intra-grid dependency coordination. The
 *   signaling to peers (inter- or intra-node) happens separately after packing.
 * - The destination pointer can be a local send buffer or a peer-accessible pointer
 *   (NVLINK) computed by the caller.
 *
 * \param[in,out] gm_dataDest                 Destination buffer for packed coordinates
 * \param[in]     sendSize                    Number of elements to pack
 * \param[in,out] data                        Full coordinate array (source)
 * \param[in]     map                         Index map selecting elements to pack
 * \param[in]     usePBC                      Whether to apply periodic shift
 * \param[in]     coordinateShift             Shift vector to add when usePBC is true
 * \param[in]     currPulse                   Current pulse/dimension index
 * \param[in]     signalXCounter              Signal value identifying this timestep/pulse
 * \param[in,out] xGridSync                   Per-pulse grid sync counter (not modified here)
 * \param[in]     signalReceiverRankXCurr     Pointer to receiver-rank X signal for this pulse
 * \param[in]     sendRank                    Destination rank id for X
 * \param[in]     recvSize                    Number of elements expected to be received
 * \param[in]     totalNumPulses              Total number of pulses across all dimensions
 * \param[in]     threadIndex                 Linearized thread index for grid-stride loops
 * \param[in]     numOfThreadBlocksRequired   Number of CTAs launched for this pulse
 * \param[in]     dependencyAtomOffset        Threshold index separating dependency subset
 */
inline __device__ void packPeerPutXNonTma(float3*         gm_dataDest,
                                          const int&      sendSize,
                                          float3*         data,
                                          const int*      map,
                                          const bool&     usePBC,
                                          const float3&   coordinateShift,
                                          const int&      currPulse,
                                          const uint64_t& signalXCounter,
                                          uint32_t*       xGridSync,
                                          uint64_t*       signalReceiverRankXCurr,
                                          const int&      sendRank,
                                          const int&      recvSize,
                                          const int&      totalNumPulses,
                                          const int       threadIndex,
                                          const int&      numOfThreadBlocksRequired,
                                          const int&      dependencyAtomOffset)
{
    const int  gridStride = blockDim.x * gridDim.x;
    const bool isThread0  = (threadIdx.x == 0);

    if (currPulse > 0)
    {
        // First pass: pack elements below dependency threshold
        int hasDependencyAtoms = 0;
        packSubset</*packAll*/ false, /*packLessThan*/ true>(
                gm_dataDest, data, map, sendSize, gridStride, threadIndex, usePBC, coordinateShift, dependencyAtomOffset, hasDependencyAtoms);

        if (isThread0)
        {
            // wait for the previous pulse/dimension data to be received.
            for (int i = currPulse; i > 0; i--)
            {
                while (loadRelaxedSysAsm(signalReceiverRankXCurr - i) != signalXCounter)
                    ;
            }
        }
        const int need_to_wait = __any_sync(c_fullWarpMask, hasDependencyAtoms > 0);
        __syncthreads();
        if (need_to_wait)
        {
            // Second pass: pack the deferred elements (>= threshold)
            int ignored = 0;
            packSubset</*packAll*/ false, /*packLessThan*/ false>(
                    gm_dataDest, data, map, sendSize, gridStride, threadIndex, usePBC, coordinateShift, dependencyAtomOffset, ignored);
        }
    }
    else
    {
        // No dependency between pulses; pack all without threshold checks
        int ignored = 0;
        packSubset</*packAll*/ true, /*packLessThan*/ true>(
                gm_dataDest, data, map, sendSize, gridStride, threadIndex, usePBC, coordinateShift, 0 /*unused*/, ignored);
    }
    __syncthreads();
}

// Helper function to handle signal synchronization
// Note: this function is supposed to be called by a single thread in the threadblock.
template<bool isRemoteInterNode>
inline __device__ void manageXHaloExchangeSync(uint32_t*       xGridSync,
                                               uint64_t*       signalReceiverRankXCurr,
                                               const uint64_t& signalXCounter,
                                               const int&      sendRank,
                                               const int&      recvSize,
                                               const int&      sendSize,
                                               const int&      currPulse,
                                               const int&      totalNumPulses,
                                               const int&      numOfThreadBlocksRequired,
                                               float*          recvPtr,
                                               const float*    dataPacked)
{
#if GMX_NVSHMEM
    // Increment per-pulse arrive counter with release semantics which tracks how many
    // threadBlocks are done; returns old value.
    uint32_t old_arrive = atomicIncReleaseGpu(xGridSync, (numOfThreadBlocksRequired - 1));
    // Last arriving block will signal completion for this pulse.
    if (old_arrive == (numOfThreadBlocksRequired - 1))
    {
        if constexpr (isRemoteInterNode)
        {

            nvshmem_float_put_signal_nbi(recvPtr,
                                         dataPacked,
                                         sendSize * 3,
                                         signalReceiverRankXCurr,
                                         signalXCounter,
                                         NVSHMEM_SIGNAL_SET,
                                         sendRank);
        }
        else
        {
            uint64_t* sendRankSignalPtr =
                    reinterpret_cast<uint64_t*>(nvshmem_ptr(signalReceiverRankXCurr, sendRank));
            // Use system-scoped release store to publish the completion signal after ensuring
            // all prior global writes from all threadBlocks are visible at system scope to
            // the remote peer.
            storeReleaseSysAsm(sendRankSignalPtr, signalXCounter);
        }

        // We need to wait for the signal completion before exiting the kernel
        // in order to make sure the packed data received from recvRank's is complete.
        if (recvSize > 0 && currPulse == totalNumPulses - 1)
        {
            while (loadAcquireSysAsm(signalReceiverRankXCurr) != signalXCounter)
                ;
        }
    }
#endif
}

/*! \brief Fused pack-and-send kernel for X halo exchange (all pulses in one launch)
 *
 * Processes all pulse/dimension entries in a single kernel launch. Packs X coordinates
 * into the per-pulse send buffer or writes directly to the peer buffer when a peer
 * pointer is available (via NVLink). Performs inter/intra-node
 * synchronization via NVSHMEM signals. For inter-node paths, a put+signal is
 * issued once all of the grid's thread blocks have finished writing the packed
 * data on per pulse basis. For intra-node paths (NVLink), the remote
 * signal location is written directly after ensuring system-wide visibility of
 * prior writes. The kernel iterates over pulse/dimension entries assigned to
 * the current grid.y index.
 *
 * Mapping of blocks/threads to work:
 * - Pulses/dimensions: gridDim.y == totalNumPulses. Each block row (blockIdx.y)
 *   selects a single pulse; all CTAs with the same blockIdx.y cooperate on that pulse.
 * - CTAs within a pulse: gridDim.x CTAs partition the work for that pulse via
 *   grid-stride loops; numOfThreadBlocksRequired == gridDim.x.
 * - Threads to atoms: each thread uses a linear index
 *   threadIndex = blockIdx.x * blockDim.x + threadIdx.x and iterates
 *   idx += blockDim.x * gridDim.x. Each iteration packs atom map[idx], where
 *   map is the per-pulse device index map (HaloExchangeData::d_indexMap).
 *
 * \param[in]  data                        Full array of coordinates on this rank
 * \param[in]  dataPacked                  Array of halo exchange entries (one per pulse)
 * \param[in]  signalReceiverRankX         Base pointer to coordinates NVSHMEM signal objects
 * \param[in]  signalXCounter              Signal counter value to write per entry
 * \param[in,out] xGridSync_               Per-pulse device counters used to coordinate CTAs
 * \param[in]  totalNumPulses          Total number of pulse/dimension entries
 * \param[in]  coordinateShiftFirstBoxVec  Shift vector for box dimension 0
 * \param[in]  coordinateShiftSecondBoxVec Shift vector for box dimension 1
 * \param[in]  coordinateShiftThirdBoxVec  Shift vector for box dimension 2
 * \param[in]  dependencyAtomOffset        Atom index boundary for dependency ordering between pulses
 */
__global__ void fusedPulsesPackXAndSendKernel(float3* __restrict__ data,
                                              FusedGpuHaloExchange::HaloExchangeData* __restrict__ dataPacked,
                                              uint64_t* __restrict__ signalReceiverRankX,
                                              const uint64_t signalXCounter,
                                              uint32_t* __restrict__ xGridSync_,
                                              const int    totalNumPulses,
                                              const float3 coordinateShiftFirstBoxVec,
                                              const float3 coordinateShiftSecondBoxVec,
                                              const float3 coordinateShiftThirdBoxVec,
                                              const int    dependencyAtomOffset)
{
#if GMX_NVSHMEM
    namespace ptx          = cuda::ptx;
    int        threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
    const bool isThread0   = (threadIdx.x == 0);

    for (int currPulse = blockIdx.y; currPulse < totalNumPulses; currPulse += gridDim.y)
    {
        FusedGpuHaloExchange::HaloExchangeData haloExchangeData = dataPacked[currPulse];
        const int                              sendRank         = haloExchangeData.sendRankX;
        const int                              recvRank         = haloExchangeData.recvRankX;
        const int                              sendSize         = haloExchangeData.xSendSize;
        const int                              recvSize         = haloExchangeData.xRecvSize;
        const int atomOffsetXInSendRank = haloExchangeData.nvshmemData.putAtomOffsetInReceiverRankXBuf_;
        const int    boxDimensionIndex       = haloExchangeData.boxDimensionIndex;
        const float3 coordinateShift         = (boxDimensionIndex == 0) ? coordinateShiftFirstBoxVec
                                               : (boxDimensionIndex == 1) ? coordinateShiftSecondBoxVec
                                                                          : coordinateShiftThirdBoxVec;
        const int*   map                     = haloExchangeData.d_indexMap;
        uint64_t*    signalReceiverRankXCurr = signalReceiverRankX + currPulse;
        uint32_t*    xGridSync               = xGridSync_ + currPulse;
        float3* gm_remotePtr = reinterpret_cast<float3*>(haloExchangeData.nvshmemData.remoteXPeerPutPtr);
        // When the remote is directly accessible we query and store the remote pointer at setup time.
        const bool remoteIsDirectlyAccessible = (gm_remotePtr != nullptr);
        const int  numOfThreadBlocksRequired  = gridDim.x;

        if (sendSize > 0)
        {
            const bool usePBC = haloExchangeData.usePBC;
            float3*    dest   = remoteIsDirectlyAccessible ? (gm_remotePtr + atomOffsetXInSendRank)
                                                           : (float3*)haloExchangeData.d_sendBuf;
            packPeerPutXNonTma(dest,
                               sendSize,
                               data,
                               map,
                               usePBC,
                               coordinateShift,
                               currPulse,
                               signalXCounter,
                               xGridSync,
                               signalReceiverRankXCurr,
                               sendRank,
                               recvSize,
                               totalNumPulses,
                               threadIndex,
                               numOfThreadBlocksRequired,
                               dependencyAtomOffset);
            if (isThread0)
            {
                if (!remoteIsDirectlyAccessible)
                {
                    float*       recvPtr = reinterpret_cast<float*>(&data[atomOffsetXInSendRank]);
                    const float* dataPackedPtr = reinterpret_cast<const float*>(dest);
                    manageXHaloExchangeSync<true>(xGridSync,
                                                  signalReceiverRankXCurr,
                                                  signalXCounter,
                                                  sendRank,
                                                  recvSize,
                                                  sendSize,
                                                  currPulse,
                                                  totalNumPulses,
                                                  numOfThreadBlocksRequired,
                                                  recvPtr,
                                                  dataPackedPtr);
                }
                else
                {
                    manageXHaloExchangeSync<false>(xGridSync,
                                                   signalReceiverRankXCurr,
                                                   signalXCounter,
                                                   sendRank,
                                                   recvSize,
                                                   sendSize,
                                                   currPulse,
                                                   totalNumPulses,
                                                   numOfThreadBlocksRequired,
                                                   nullptr,
                                                   nullptr);
                }
            }
        }
        else if (threadIdx.x == 0 && blockIdx.x == (gridDim.x - 1) && recvSize > 0)
        {
            // We need to wait for the signal completion before exiting the kernel
            // in order to make sure the packed data received from recvRank's is complete.
            // nvshmem_signal_wait_until(signalReceiverRankX, NVSHMEM_CMP_EQ, signalXCounter);
            while (loadAcquireSysAsm(signalReceiverRankXCurr) != signalXCounter)
                ;
        }
    }
#endif
}

/*!
 * \brief Waits for NVSHMEM signal and then unpacks received forces.
 *
 * Details:
 * - Spins (thread 0 only) until the signal at \p signalReceiverRankFCurr
 *   reaches \p signalReceiverRankFCounter, establishing visibility of the remote writes.
 * - Uses a grid-stride loop to copy/accumulate values from \p gm_dataSrc into
 *   \p gm_dataDest using the provided \p map.
 * - When \p accumulate is true, atomically accumulates the forces else overwrites.
 * - light-weight L1 prefetch for indexMap until spinning on signal wait.
 */
__device__ void unpackForcesAfterSignal(const int       threadIndex,
                                        const float3*   gm_dataSrc,
                                        float3*         gm_dataDest,
                                        const int*      map,
                                        const int       recvSize,
                                        uint64_t*       signalReceiverRankFCurr,
                                        uint64_t        signalReceiverRankFCounter,
                                        const bool      accumulate,
                                        const int&      currPulse,
                                        const int&      totalNumPulses,
                                        const uint32_t* d_fGridSync_)
{

    int gridStride = blockDim.x * gridDim.x;
    if (threadIdx.x == 0)
    {
        while (loadRelaxedSysAsm(signalReceiverRankFCurr) != signalReceiverRankFCounter)
            ;
    }

    asm("prefetch.global.L1 [%0];" ::"l"(map + threadIndex) : "memory");
    asm("prefetch.global.L1 [%0];" ::"l"(map + threadIndex + gridStride) : "memory");

    __syncthreads();

    for (int i = threadIndex; i < recvSize; i += gridStride)
    {
        int    mapIndex = map[i];
        float3 srcVal   = gm_dataSrc[i];
        if (accumulate)
        {
            float* gm_dataDest_ptr = (float*)(gm_dataDest + mapIndex);
            atomicAdd(gm_dataDest_ptr, srcVal.x);
            atomicAdd(gm_dataDest_ptr + 1, srcVal.y);
            atomicAdd(gm_dataDest_ptr + 2, srcVal.z);
        }
        else
        {
            gm_dataDest[mapIndex] = srcVal;
        }
    }
}


/*!
 * \brief Fused NVSHMEM force unpack kernel for F halo exchange.
 *
 * Algorithm overview:
 * - The kernel runs with grid.y = total number of pulses. Each
 *   block in grid.y processes one pulse identified by \p currPulse, iterating in
 *   descending order to honor dependencies on later pulses.
 * - For each pulse:
 *   1) Determine communication path:
 *      - Inter-node (NIC): when the peer's destination pointer is not directly accessible
 *        (remoteDstIsDirectlyAccessible == false), this rank issues a PUT of the packed
 *        force slice to the peer and signals it (nvshmem_float_put_signal_nbi).
 *      - Intra-node (NVLink / same node): when the peer's destination is directly accessible,
 *        this rank signals the peer to GET (nvshmem_ptr() non-null) and the peer pulls the data.
 *   2) Locally unpack the received forces for this pulse:
 *      - If a peer GET source pointer is available, read directly from that remote source + offset.
 *      - Otherwise read from this rank's receive buffer.
 *      - Call unpackForcesAfterSignal() which first waits on the pulse-specific NVSHMEM signal and
 *        then performs a grid-stride accumulate-or-overwrite into the global force array.
 *   3) If there exists a dependent earlier pulse (nextPulse = currPulse - 1), the last arriving
 *      thread block for the current pulse signals readiness for nextPulse. This is implemented with
 *      a release-scoped atomicInc on d_fGridSync_[currPulse] and a spin with acquire loads on
 *      the nextPulse side.
 *
 * Synchronization and satisfying pulse ordering:
 * - NVSHMEM signaling between ranks:
 *   - First pulse that signals the peer within a dim/pulse sequence can use relaxed system store
 *     since prior data for it is visible at sys scope from previous NBF non-local kernel; subsequent signaling uses release system
 *     stores to ensure all prior writes are globally visible before peers observe the signal.
 * - Intra-kernel inter-pulse ordering:
 *   - The last thread block to finish a pulse performs atom.inc.release.gpu on d_fGridSync_.
 *   - Consumers of the next pulse spin using load.acquire.gpu until the expected count is reached.
 *
 * Path selection booleans:
 *   - remoteDstIsDirectlyAccessible: peer destination reachable via NVLink/same node.
 *   - remoteSrcIsDirectlyAccessible: peer source pointer available for direct GET.
 *
 * \param[out] data                         Full array of forces to accumulate into
 * \param[in]  dataPacked                   Array of halo exchange entries (one per pulse)
 * \param[in]  signalReceiverRankF          Base pointer for forces NVSHMEM signal objects
 * \param[in]  signalReceiverRankFCounter   Monotonic signal counter value to write per pulse
 * \param[in]  totalNumPulses              Total number of pulses
 * \param[in]  accumulateForces             Whether to accumulate (atomic) or overwrite forces
 * \param[in,out] d_fGridSync_              Per-pulse device counters used for inter-pulse sync
 */
__global__ void fusedUnPackFRecvBufNvshmemKernel(float3* __restrict__ data,
                                                 FusedGpuHaloExchange::HaloExchangeData* __restrict__ dataPacked,
                                                 uint64_t* __restrict__ signalReceiverRankF,
                                                 const uint64_t signalReceiverRankFCounter,
                                                 int            totalNumPulses,
                                                 bool           accumulateForces,
                                                 uint32_t* __restrict__ d_fGridSync_)
{
#if GMX_NVSHMEM
    const int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;

    for (int currPulse = blockIdx.y; currPulse >= 0; currPulse -= gridDim.y)
    {
        FusedGpuHaloExchange::HaloExchangeData haloExchangeData = dataPacked[currPulse];
        // Reverse the send and recv ranks/sizes of X is send and recv for F
        const int sendRank = haloExchangeData.recvRankX;
        const int recvRank = haloExchangeData.sendRankX;
        const int sendSize = haloExchangeData.xRecvSize;
        const int recvSize = haloExchangeData.xSendSize;

        const int atomOffset = haloExchangeData.atomOffset;
        const int remoteGetAtomOffset = haloExchangeData.nvshmemData.putAtomOffsetInReceiverRankXBuf_;
        const float3* remoteForcePeerGetPtr =
                (const float3*)haloExchangeData.nvshmemData.remoteForcePeerGetPtr;
        float3* remoteForcePeerPutPtr = (float3*)haloExchangeData.nvshmemData.remoteForcePeerPutPtr;
        const bool remoteDstIsDirectlyAccessible = (remoteForcePeerPutPtr != nullptr);
        float3*    recvBuf                       = (float3*)haloExchangeData.d_recvBuf;

        // Accumulation behavior:
        // - Accumulate when there are CPU-side force contributions.
        // - If there are no CPU-side contributions, accumulate only for 2D/3D domain decomposition for all pulses,
        //   and for 1D only from the second pulse onward (not on pulse 0 in 1D).
        // const bool accumulate = haloExchangeData.accumulateForces ? true : accumulateForces;
        const bool accumulate              = accumulateForces && (currPulse == 0) ? accumulateForces
                                                                                  : haloExchangeData.accumulateForces;
        const int* map                     = haloExchangeData.d_indexMap;
        uint64_t*  signalReceiverRankFCurr = signalReceiverRankF + currPulse;

        if (threadIndex == 0 && sendSize > 0 && currPulse == (totalNumPulses - 1))
        {
            // Inter-node path (NIC): use nvshmem put+signal to transfer and notify the peer.
            if (!remoteDstIsDirectlyAccessible)
            {
                nvshmem_float_put_signal_nbi(reinterpret_cast<float*>(recvBuf),
                                             reinterpret_cast<float*>(&data[atomOffset]),
                                             sendSize * 3,
                                             signalReceiverRankFCurr,
                                             signalReceiverRankFCounter,
                                             NVSHMEM_SIGNAL_SET,
                                             sendRank);
            }
            else
            {
                // Intra-node path (NVLink): signal the peer to pull (GET) the data directly.
                uint64_t* remoteSignalReceiverRankFCurr =
                        reinterpret_cast<uint64_t*>(nvshmem_ptr(signalReceiverRankFCurr, sendRank));
                GMX_ASSERT(remoteSignalReceiverRankFCurr != nullptr,
                           "nvshmem_ptr returned null for peer signal pointer");
                //  As this is the first dim/pulse we can make use of relaxed memory order to signal the sendRank.
                //  As there is no prior data written by this rank to be made visible at system scope.
                stRelaxedSysAsm(remoteSignalReceiverRankFCurr, signalReceiverRankFCounter);
            }
        }

        if (recvSize > 0)
        {
            const bool    remoteSrcIsDirectlyAccessible = (remoteForcePeerGetPtr != nullptr);
            const float3* gm_dataSrc                    = remoteSrcIsDirectlyAccessible
                                                                  ? (remoteForcePeerGetPtr + remoteGetAtomOffset)
                                                                  : recvBuf;
            unpackForcesAfterSignal(threadIndex,
                                    gm_dataSrc,
                                    data,
                                    map,
                                    recvSize,
                                    signalReceiverRankFCurr,
                                    signalReceiverRankFCounter,
                                    accumulate,
                                    currPulse,
                                    totalNumPulses,
                                    d_fGridSync_);
        }

        if (currPulse > 0)
        {
            // we skip the last grid.y as the last that is 0th pulse/dim doesn't have any dependent threadblocks
            __syncthreads();
            if (threadIdx.x == 0)
            {
                uint32_t old_val = atomicIncReleaseGpu(d_fGridSync_ + currPulse, gridDim.x);

                if (old_val == gridDim.x - 1)
                {
                    int nextPulse = currPulse - 1;
                    if (nextPulse >= 0)
                    {
                        FusedGpuHaloExchange::HaloExchangeData haloExchangeData = dataPacked[nextPulse];
                        // Reverse the send and recv ranks/sizes of X is send and recv for F
                        const int sendRank   = haloExchangeData.recvRankX;
                        const int sendSize   = haloExchangeData.xRecvSize;
                        const int atomOffset = haloExchangeData.atomOffset;

                        Float3* remoteForcePeerPutPtr = haloExchangeData.nvshmemData.remoteForcePeerPutPtr;
                        const bool remoteDstIsDirectlyAccessible = (remoteForcePeerPutPtr != nullptr);
                        Float3*   recvBuf                 = haloExchangeData.d_recvBuf;
                        uint64_t* signalReceiverRankFNext = signalReceiverRankF + nextPulse;

                        if (sendSize > 0)
                        {
                            uint64_t* remoteSignalReceiverRankFCurr = reinterpret_cast<uint64_t*>(
                                    nvshmem_ptr(signalReceiverRankFNext, sendRank));
                            GMX_ASSERT(remoteSignalReceiverRankFCurr != nullptr,
                                       "nvshmem_ptr returned null for peer signal pointer");
                            for (int pulseId = currPulse + 1; pulseId < totalNumPulses; pulseId++)
                            {
                                const uint32_t* syncOnPrevPulse = d_fGridSync_ + pulseId;
                                while (loadAcquireGpuAsm(syncOnPrevPulse) != gridDim.x)
                                    ;
                            }
                            if (!remoteDstIsDirectlyAccessible)
                            {
                                // for inter-node comm over NIC issue nvshmem_float_put_signal_nbi
                                nvshmem_float_put_signal_nbi(reinterpret_cast<float*>(recvBuf),
                                                             reinterpret_cast<float*>(&data[atomOffset]),
                                                             sendSize * 3,
                                                             signalReceiverRankFNext,
                                                             signalReceiverRankFCounter,
                                                             NVSHMEM_SIGNAL_SET,
                                                             sendRank);
                            }
                            else
                            {
                                // For NVLink comm we just signal the sender rank using system-scoped
                                // release store to get(pull) the forces as they are ready.
                                storeReleaseSysAsm(remoteSignalReceiverRankFCurr, signalReceiverRankFCounter);
                            }
                        }
                    }
                    if (currPulse == 1)
                    {
                        // reset the grid sync before kernel exit
                        for (int i = 0; i < totalNumPulses; i++)
                        {
                            d_fGridSync_[i] = 0;
                        }
                    }
                }
            }
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
    // Avoid launching kernel when there is no work to do
    if (xSendSize_ > 0)
    {
        auto       kernelFn   = usePBC_ ? packSendBufKernel<true> : packSendBufKernel<false>;
        const auto kernelArgs = prepareGpuKernelArguments(
                kernelFn, config, &sendBuf, &d_x, &indexMap, &xSendSize_, &coordinateShift);
        launchGpuKernel(kernelFn, config, *haloStream_, nullptr, "Domdec GPU Apply X Halo Exchange", kernelArgs);
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

    if (fRecvSize_ > 0)
    {
        auto kernelFn = accumulateForces ? unpackRecvBufKernel<true> : unpackRecvBufKernel<false>;
        const auto kernelArgs =
                prepareGpuKernelArguments(kernelFn, config, &d_f, &recvBuf, &indexMap, &fRecvSize_);

        launchGpuKernel(kernelFn, config, *haloStream_, nullptr, "Domdec GPU Apply F Halo Exchange", kernelArgs);
    }
}

void FusedGpuHaloExchange::launchPackXKernel(const matrix box)
{
#if GMX_NVSHMEM
    // Configure kernel launch parameters
    KernelLaunchConfig config;
    config.blockSize[0]     = c_fusedKernelsThreadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = maxGridXSize_;
    config.gridSize[1]      = totalNumPulses_;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    float3* d_x                 = asFloat3(sharedBuffers_.d_x);
    auto    packNvShmemkernelFn = fusedPulsesPackXAndSendKernel;

    // The coordinateShift changes between steps when we have
    // performed a DD partition, or have updated the box e.g. when
    // performing pressure coupling. So, for simplicity, the box
    // is used every step to pass the shift vector as an argument of
    // the packing kernel.
    const float3 coordinateShiftFirstBoxVec{ box[0][XX], box[0][YY], box[0][ZZ] };
    const float3 coordinateShiftSecondBoxVec{ box[1][XX], box[1][YY], box[1][ZZ] };
    const float3 coordinateShiftThirdBoxVec{ box[2][XX], box[2][YY], box[2][ZZ] };
    const auto   kernelArgs = prepareGpuKernelArguments(packNvShmemkernelFn,
                                                      config,
                                                      &d_x,
                                                      &d_haloExchangeData_,
                                                      &sharedBuffers_.d_signalReceiverRankX_,
                                                      &signalReceiverRankXCounter_,
                                                      &d_xGridSync_,
                                                      &totalNumPulses_,
                                                      &coordinateShiftFirstBoxVec,
                                                      &coordinateShiftSecondBoxVec,
                                                      &coordinateShiftThirdBoxVec,
                                                      &haloExchangeData_[0].atomOffset);

    launchGpuKernel(
            packNvShmemkernelFn, config, *haloStream_, nullptr, "Domdec GPU Apply X Halo Exchange", kernelArgs);
    signalReceiverRankXCounter_++;
#endif
}

void FusedGpuHaloExchange::launchUnpackFKernel(bool accumulateForces)
{
#if GMX_NVSHMEM
    KernelLaunchConfig config;
    config.blockSize[0]     = c_fusedKernelsThreadsPerBlock;
    config.blockSize[1]     = 1;
    config.blockSize[2]     = 1;
    config.gridSize[0]      = maxGridFSize_;
    config.gridSize[1]      = totalNumPulses_;
    config.gridSize[2]      = 1;
    config.sharedMemorySize = 0;

    float3* d_f = asFloat3(sharedBuffers_.d_f);

    auto kernelFn = fusedUnPackFRecvBufNvshmemKernel;

    const auto kernelArgs = prepareGpuKernelArguments(kernelFn,
                                                      config,
                                                      &d_f,
                                                      &d_haloExchangeData_,
                                                      &sharedBuffers_.d_signalReceiverRankF_,
                                                      &signalReceiverRankFCounter_,
                                                      &totalNumPulses_,
                                                      &accumulateForces,
                                                      &d_fGridSync_);

    launchGpuKernel(kernelFn, config, *haloStream_, nullptr, "Domdec GPU Apply F Halo Exchange", kernelArgs);

    signalReceiverRankFCounter_++;
#endif
}

} // namespace gmx
