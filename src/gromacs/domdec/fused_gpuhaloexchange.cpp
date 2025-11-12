/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief Implements NVSHMEM based Fused PP GPU halo exchange.
 *
 *
 * \author Mahesh Doijade<mdoijade@nvidia.com>
 *
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "fused_gpuhaloexchange.h"

#include "config.h"

#include <cstdint>

#include <numeric>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/commrec.h" // t_commrec definition
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "domdec_internal.h"
#include "domdec_struct.h" // gmx_domdec_t, gmx_domdec_comm_t
#include "halo_plan_utils.h"

#if GMX_NVSHMEM
#    include <nvshmem.h>
#    include <nvshmemx.h>
#endif

namespace gmx
{

namespace detail
{

template<std::size_t alignment>
inline bool is_sufficiently_aligned(const void* ptr)
{
    // alignment must be power of two
    return (reinterpret_cast<std::uintptr_t>(ptr) & (alignment - 1U)) == 0U;
}
} // namespace detail

FusedGpuHaloExchange::FusedGpuHaloExchange(const DeviceContext& deviceContext,
                                           gmx_wallcycle*       wcycle,
                                           MPI_Comm             mpi_comm_mysim,
                                           MPI_Comm             mpi_comm_mysim_world) :
    haloStream_(new DeviceStream(deviceContext, DeviceStreamPriority::High, false)),
    deviceContext_(deviceContext),
    wcycle_(wcycle),
    signalReceiverRankXCounter_(0),
    signalReceiverRankFCounter_(0),
    enableFusedForceKernelSync_(false),
    mpi_comm_mysim_(mpi_comm_mysim),
    mpi_comm_mysim_world_(mpi_comm_mysim_world)
{
    enableFusedForceKernelSync_ = (getenv("GMX_NVSHMEM_FUSED_FORCE_KERNEL_SYNC") != nullptr);
}

FusedGpuHaloExchange::~FusedGpuHaloExchange()
{
    // Free grid sync buffers
    freeDeviceBuffer(&d_fGridSync_);
    d_fGridSyncSize_      = -1;
    d_fGridSyncSizeAlloc_ = -1;

    freeDeviceBuffer(&d_xGridSync_);
    d_xGridSyncSize_      = -1;
    d_xGridSyncSizeAlloc_ = -1;
    freeDeviceBuffer(&d_haloExchangeData_);
}

GpuEventSynchronizer* FusedGpuHaloExchange::launchAllCoordinateExchanges(const matrix box,
                                                                         GpuEventSynchronizer* dependencyEvent)
{
    wallcycle_start(wcycle_, WallCycleCounter::LaunchGpuPp);

    dependencyEvent->enqueueWaitEvent(*haloStream_);

    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuMoveX);
    // relies on device-side kernels from existing GPU halo exchange
    // to be wrapped by private helpers in this class; not shown here
    // as only the provided interface is required
    launchPackXKernel(box);
    coordinateHaloLaunched_.markEvent(*haloStream_);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuMoveX);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);

    return &coordinateHaloLaunched_;
}

GpuEventSynchronizer* FusedGpuHaloExchange::launchAllForceExchanges(
        bool                                           accumulateForces,
        FixedCapacityVector<GpuEventSynchronizer*, 2>* dependencyEvents)
{
    while (!dependencyEvents->empty())
    {
        auto* dependency = dependencyEvents->back();
        dependency->enqueueWaitEvent(*haloStream_);
        dependencyEvents->pop_back();
    }

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuMoveF);

    launchUnpackFKernel(accumulateForces);
    forceHaloLaunched_.markEvent(*haloStream_);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuMoveF);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);

    return &forceHaloLaunched_;
}

void FusedGpuHaloExchange::destroyAllHaloExchangeBuffers()
{
    // TODO: Re-enable collective free of symmetric buffers once freeing hang with PME is resolved
    // freeDeviceBuffer(&d_unifiedRecvBuf_);
    // freeDeviceBuffer(&d_unifiedSendBuf_);
}

void FusedGpuHaloExchange::allocateAndCopyHaloExchangeData()
{
    reallocateDeviceBuffer(&d_haloExchangeData_,
                           haloExchangeData_.size(),
                           &d_haloExchangeDataSize_,
                           &d_haloExchangeDataSizeAlloc_,
                           deviceContext_);

    copyToDeviceBuffer(&d_haloExchangeData_,
                       haloExchangeData_.data(),
                       0,
                       haloExchangeData_.size(),
                       *haloStream_,
                       GpuApiCallBehavior::Async,
                       nullptr);
}

void FusedGpuHaloExchange::allocateUnifiedHaloBuffers()
{
    constexpr int alignmentBytes = FusedGpuHaloExchange::c_haloEntryAlignBytes;
    static_assert((alignmentBytes & (alignmentBytes - 1)) == 0,
                  "alignmentBytes must be a power of two");
    constexpr int bytesPerElem = sizeof(Float3);
    auto          alignUpBytes = [](size_t value, size_t alignment) -> size_t
    { return (value + alignment - 1) / alignment * alignment; };

    // Find a multiple of c_haloEntryAlignBytes that is also a multiple of element size
    constexpr std::uint64_t elemAlignBytes = std::lcm(alignmentBytes, bytesPerElem);

    std::uint64_t maxPerPulseSendBytes = 0;
    std::uint64_t maxPerPulseRecvBytes = 0;
    for (const auto& e : haloExchangeData_)
    {
        maxPerPulseSendBytes =
                std::max(maxPerPulseSendBytes, static_cast<size_t>(e.xSendSize) * bytesPerElem);
        maxPerPulseRecvBytes =
                std::max(maxPerPulseRecvBytes, static_cast<size_t>(e.xRecvSize) * bytesPerElem);
    }
    maxPerPulseSendBytes = alignUpBytes(maxPerPulseSendBytes, elemAlignBytes);
    maxPerPulseRecvBytes = alignUpBytes(maxPerPulseRecvBytes, elemAlignBytes);
    //! max-reduced per-pulse extent across ranks (elements), unified for send/recv
    std::uint64_t localSendRecvMax =
            gmx::divideRoundUp(std::max(maxPerPulseSendBytes, maxPerPulseRecvBytes),
                               static_cast<std::uint64_t>(bytesPerElem));
    // Use max across all ranks to ensure symmetric allocation sizes
    std::uint64_t globalSendRecvMax = 0;
#if GMX_MPI
    MPI_Allreduce(&localSendRecvMax, &globalSendRecvMax, 1, MPI_UINT64_T, MPI_MAX, mpi_comm_mysim_world_);
#endif
    std::uint64_t totalSendRecvSize = globalSendRecvMax;

    constexpr bool symmetricAlloc = true;
    // We allocate unified send/recv buffers for totalNumPulses_ entries to ensure
    // symmetric per-pulse/per-dimension layout across all ranks. The per-pulse pointers
    // are also assigned here. Without symmetric sizing here,
    // NVSHMEM pointer arithmetic on the remote peer can address beyond the intended
    // region and overwrite memory on the receiver. This approach lets us call
    // MPI_Allreduce and perform symmetric allocation only once, at the cost of some
    // extra total allocation.
    reallocateDeviceBuffer(&d_unifiedSendBuf_,
                           totalSendRecvSize * totalNumPulses_,
                           &unifiedSendSize_,
                           &unifiedSendCapacity_,
                           deviceContext_,
                           symmetricAlloc);
    reallocateDeviceBuffer(&d_unifiedRecvBuf_,
                           totalSendRecvSize * totalNumPulses_,
                           &unifiedRecvSize_,
                           &unifiedRecvCapacity_,
                           deviceContext_,
                           symmetricAlloc);

    // Point each pulse's send/recv buffer into unified buffers using aligned offsets
    int               idx = 0;
    const std::string msg =
            gmx::formatString("Unified send entry pointer is not %d-byte aligned", alignmentBytes);
    for (auto& e : haloExchangeData_)
    {
        e.d_sendBuf = d_unifiedSendBuf_ + idx * totalSendRecvSize;
        GMX_RELEASE_ASSERT(gmx::detail::is_sufficiently_aligned<alignmentBytes>(e.d_sendBuf), msg.c_str());
        e.d_recvBuf = d_unifiedRecvBuf_ + idx * totalSendRecvSize;
        GMX_RELEASE_ASSERT(gmx::detail::is_sufficiently_aligned<alignmentBytes>(e.d_recvBuf), msg.c_str());
        idx++;
    }
}

GpuEventSynchronizer* FusedGpuHaloExchange::getForcesReadyOnDeviceEvent()
{
    return &forceHaloLaunched_;
}

void FusedGpuHaloExchange::reinitAllHaloExchanges(const t_commrec&       cr,
                                                  DeviceBuffer<RVec>     d_coordinatesBuffer,
                                                  DeviceBuffer<RVec>     d_forcesBuffer,
                                                  DeviceBuffer<uint64_t> d_syncBase,
                                                  int                    totalNumPulses)
{
    wallcycle_start_nocount(wcycle_, WallCycleCounter::Domdec);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::DDGpu);
    // Set shared buffers
    sharedBuffers_.d_x = d_coordinatesBuffer;
    sharedBuffers_.d_f = d_forcesBuffer;

    // Set NVSHMEM signals
    d_ppHaloExSyncBase_     = d_syncBase;
    ppHaloExPerSyncBufSize_ = totalNumPulses;
    totalNumPulses_         = totalNumPulses;

    maxGridXSize_ = 0;
    maxGridFSize_ = 0;
    // Build per-dimension/pulse entries mirroring Impl::reinitHalo
    haloExchangeData_.resize(totalNumPulses_);

    const gmx_domdec_comm_t& comm = *cr.dd->comm;
    // Host-side storage for device barriers without requiring copy/move
    DeviceBarrier* barriers               = nullptr;
    int            numBarriersConstructed = 0;
    GMX_RELEASE_ASSERT(totalNumPulses_ > 0, "Total number of pulses must be > 0");
    barriers = static_cast<DeviceBarrier*>(
            ::operator new(static_cast<size_t>(totalNumPulses_) * sizeof(DeviceBarrier)));

    int idxEntry = 0;
    for (int d = 0; d < cr.dd->ndim; d++)
    {
        const int  dimIndex  = d;
        const int  sendRankX = cr.dd->neighbor[dimIndex][1];
        const int  recvRankX = cr.dd->neighbor[dimIndex][0];
        const bool usePBC    = (cr.dd->ci[cr.dd->dim[dimIndex]] == 0);

        if (usePBC && cr.dd->unitCellInfo.haveScrewPBC)
        {
            gmx_fatal(FARGS, "Error: screw is not yet supported in GPU halo exchange\n");
        }

        const gmx_domdec_comm_dim_t& cd = comm.cd[dimIndex];

        for (int pulse = 0; pulse < cd.numPulses(); pulse++)
        {
            const auto plan =
                    computeHaloPlan(comm, dimIndex, pulse, mpi_comm_mysim_, sendRankX, recvRankX);
            const int               atomOffset = plan.atomOffset;
            const int               xSendSize  = plan.xSendSize;
            const int               xRecvSize  = plan.xRecvSize;
            const gmx_domdec_ind_t* ind        = plan.ind;
            GMX_RELEASE_ASSERT(ind->index.get_allocator().pinningPolicy() == PinningPolicy::PinnedIfSupported,
                               "Array of communication indices must have been pinned");

            auto& data             = haloExchangeData_[idxEntry];
            data.xSendSize         = xSendSize;
            data.xRecvSize         = xRecvSize;
            data.atomOffset        = atomOffset;
            data.sendRankX         = sendRankX;
            data.recvRankX         = recvRankX;
            data.boxDimensionIndex = cr.dd->dim[dimIndex];
            data.usePBC            = usePBC;
            data.accumulateForces  = (pulse > 0 || cr.dd->ndim > 1);

            // Copy index map to device; set pinning policy on the original allocation
            const int mapSize = xSendSize;
            reallocateDeviceBuffer(
                    &data.d_indexMap, mapSize, &data.indexMapSize, &data.indexMapCapacity, deviceContext_);
            if (mapSize > 0)
            {
                copyToDeviceBuffer(&data.d_indexMap,
                                   ind->index.data(),
                                   0,
                                   mapSize,
                                   *haloStream_,
                                   GpuApiCallBehavior::Async,
                                   nullptr);
            }

            // Remote pointers and offsets for NVSHMEM
            data.nvshmemData.remoteXPeerPutPtr = (Float3*)nvshmem_ptr(sharedBuffers_.d_x, sendRankX);
            const int sendRankF = recvRankX;
            const int recvRankF = sendRankX;
            data.nvshmemData.remoteForcePeerGetPtr =
                    (const Float3*)nvshmem_ptr(sharedBuffers_.d_f, recvRankF);
            data.nvshmemData.remoteForcePeerPutPtr = (Float3*)nvshmem_ptr(sharedBuffers_.d_f, sendRankF);

// Exchange atom offset with neighbors for coordinate PUT destination
#if GMX_MPI
            MPI_Sendrecv(&atomOffset,
                         sizeof(int),
                         MPI_BYTE,
                         recvRankX,
                         0,
                         &data.nvshmemData.putAtomOffsetInReceiverRankXBuf_,
                         sizeof(int),
                         MPI_BYTE,
                         sendRankX,
                         0,
                         mpi_comm_mysim_,
                         MPI_STATUS_IGNORE);

#endif
            // Add 1 to size to take care of case when xSendSize_ == 0 in such case
            // we still launch the kernel to signal the remote rank recvRankX_
            // from which this rank receives packed data as (xRecvSize_ > 0)
            const int gridDimX = gmx::divideRoundUp(data.xSendSize + 1,
                                                    FusedGpuHaloExchange::c_nvshmemThreadsPerBlock);

            // Construct in-place to avoid copy/move of cuda::barrier
            if (barriers != nullptr)
            {
                new (&barriers[numBarriersConstructed]) DeviceBarrier(gridDimX);
                numBarriersConstructed++;
            }
            idxEntry++;
        }
    }

    // Allocate unified send/recv buffers and point entries into them
    allocateUnifiedHaloBuffers();

    // Copy entry metadata to device
    allocateAndCopyHaloExchangeData();

    // Initialize shared signal buffers
    GMX_RELEASE_ASSERT(d_ppHaloExSyncBase_ != nullptr,
                       "NVSHMEM PP Halo exchange requires valid signal buffer");
    // NVSHMEM signal buffer layout (contiguous, symmetric across ranks):
    //   [0 .. N-1]     -> sender-rank X signals
    //   [N .. 2N-1]    -> receiver-rank X signals
    //   [2N .. 3N-1]   -> receiver-rank F signals
    // where N == ppHaloExPerSyncBufSize_ (totalNumPulses_).
#if GMX_GPU_CUDA
    sharedBuffers_.d_signalSenderRankX_   = d_ppHaloExSyncBase_;
    sharedBuffers_.d_signalReceiverRankX_ = d_ppHaloExSyncBase_ + totalNumPulses_;
    sharedBuffers_.d_signalReceiverRankF_ = d_ppHaloExSyncBase_ + 2 * totalNumPulses_;
#endif

    reallocateDeviceBuffer(
            &d_fGridSync_, totalNumPulses_, &d_fGridSyncSize_, &d_fGridSyncSizeAlloc_, deviceContext_);
    reallocateDeviceBuffer(
            &d_xGridSync_, totalNumPulses_, &d_xGridSyncSize_, &d_xGridSyncSizeAlloc_, deviceContext_);

    // Initialize arrive-wait barriers once using the constructed array
    if (numBarriersConstructed > 0)
    {
        copyToDeviceBuffer(
                &d_xGridSync_, barriers, 0, numBarriersConstructed, *haloStream_, GpuApiCallBehavior::Sync, nullptr);
        clearDeviceBufferAsync(&d_fGridSync_, 0, totalNumPulses_, *haloStream_);
        // Destroy and release host-side barriers storage
        for (int i = 0; i < numBarriersConstructed; ++i)
        {
            barriers[i].~DeviceBarrier();
        }
        ::operator delete(barriers);
    }

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::DDGpu);
    wallcycle_stop(wcycle_, WallCycleCounter::Domdec);
}

} // namespace gmx
