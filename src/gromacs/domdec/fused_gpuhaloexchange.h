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
 * \brief Implements NVSHMEM based Fused GPU halo exchange.
 *
 *
 * \author Mahesh Doijade<mdoijade@nvidia.com>
 *
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_FUSED_GPUHALOEXCHANGE_H
#define GMX_DOMDEC_FUSED_GPUHALOEXCHANGE_H

#include <vector>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/utility/fixedcapacityvector.h"
#include "gromacs/utility/gmxmpi.h"

#if GMX_NVSHMEM
#    include <cuda/barrier>
#endif

struct gmx_wallcycle;
// Forward decls for types used in signatures without including heavy headers
struct t_commrec;
class DeviceStream;
struct gmx_domdec_t;
struct gmx_wallcycle;

namespace gmx
{

/*! \brief NVSHMEM-based fused GPU halo exchange
 *
 * Aggregates all pulses into a single structure for a fused GPU-initiated exchange
 * path using NVSHMEM. Combining coordinate (X) and force (F) metadata is a
 * by-product of this aggregation. The class:
 * - Builds per-dimension/pulse metadata and aggregates entries into one struct
 * - Allocates unified device buffers for send/recv regions and computes per-pulse
 *   offsets aligned to c_haloEntryAlignBytes for coalesced access and high-throughput
 *   transfers (e.g., NVSHMEM device puts/gets, CUDA vector/TMA async copies)
 * - Exposes device pointers per entry for pack/unpack kernels
 * - Manages NVSHMEM signal buffers via a single symmetric allocation and
 *   pointer-offset views (see implementation for exact layout)
 * - Reduces MPI collectives: unified send/recv sizes enable a single
 *   MPI_Allreduce for max sizing instead of per-pulse max-reductions
 * - Unified buffer sizes are max-reduced across ranks to guarantee identical sizes
 *
 */
class FusedGpuHaloExchange
{
public:
    /*! \brief Creates NVSHMEM-based fused GPU halo exchange object.
     *
     * \param [in] deviceContext        GPU device context
     * \param [in] wcycle               Wallclock cycle accounting
     * \param [in] mpi_comm_mysim       MPI communicator used for simulation
     * \param [in] mpi_comm_mysim_world MPI communicator involving PP + PME
     */
    FusedGpuHaloExchange(const DeviceContext& deviceContext,
                         gmx_wallcycle*       wcycle,
                         MPI_Comm             mpi_comm_mysim,
                         MPI_Comm             mpi_comm_mysim_world);
    /*! \brief Destructor. */
    ~FusedGpuHaloExchange();

    /*! \brief Launch fused coordinate (X) exchanges across all pulses.
     * \param [in] box              Coordinate box (for shifts)
     * \param [in] dependencyEvent  Dependency event for this operation
     * \returns                     Event recorded when this operation has been launched
     */
    GpuEventSynchronizer* launchAllCoordinateExchanges(const matrix          box,
                                                       GpuEventSynchronizer* dependencyEvent);
    /*! \brief Launch fused force (F) exchanges across all pulses.
     * \param [in] accumulateForces True if forces should accumulate, otherwise they are set
     * \param [in] dependencyEvents Dependency events for this operation
     * \returns                     Event recorded when this operation has been launched
     */
    GpuEventSynchronizer* launchAllForceExchanges(bool accumulateForces,
                                                  FixedCapacityVector<GpuEventSynchronizer*, 2>* dependencyEvents);

    /*! \brief (Re-)initialize fused halo exchanges for all dimensions and pulses.
     * Builds per-pulse entries, sets shared buffers, NVSHMEM signals, and prepares metadata.
     * \param [in] cr                  Communication record
     * \param [in] d_coordinatesBuffer Pointer to coordinates buffer in GPU memory
     * \param [in] d_forcesBuffer      Pointer to forces buffer in GPU memory
     * \param [in] d_syncBase          Base device pointer for NVSHMEM signal buffer
     * \param [in] totalNumPulses      Total number of pulses across all dimensions
     */
    void reinitAllHaloExchanges(const t_commrec&       cr,
                                DeviceBuffer<RVec>     d_coordinatesBuffer,
                                DeviceBuffer<RVec>     d_forcesBuffer,
                                DeviceBuffer<uint64_t> d_syncBase,
                                int                    totalNumPulses);

    /*! \brief Destroy unified/symmetric buffers used by fused halo exchange (if any). */
    void destroyAllHaloExchangeBuffers();

    // Use cuda::barrier when available, otherwise fall back to an empty type
#if GMX_NVSHMEM
    using DeviceBarrier = cuda::barrier<cuda::thread_scope_device>;
#else
    struct DeviceBarrier
    {
    };
#endif
    struct HaloExchangeData
    {
        //! The dimension index corresponding to this pulse
        int boxDimensionIndex = 0;
        //! send copy size from this pulse for X
        int xSendSize = 0;
        //! recv copy size to this pulse for X
        int xRecvSize = 0;
        //! The atom offset for receive (x) or send (f) for this pulse
        int atomOffset = 0;
        //! rank to send data to for X
        int sendRankX = 0;
        //! rank to recv data from for X
        int recvRankX = 0;
        //! whether PBC applies for this pulse
        bool usePBC = false;
        //! true when forces should accumulate for this pulse
        bool accumulateForces = false;
        //! device copy of index map for this pulse
        DeviceBuffer<int> d_indexMap = nullptr;
        //! size of per-pulse device index map
        int indexMapSize = -1;
        //! capacity of per-pulse device index map
        int indexMapCapacity = -1;
        //! device receive buffer pointer for this pulse, points into d_unifiedRecvBuf_
        DeviceBuffer<Float3> d_recvBuf = nullptr;
        //! device send buffer pointer for this pulse, points into d_unifiedSendBuf_
        DeviceBuffer<Float3> d_sendBuf = nullptr;
        struct NvshmemData
        {
            //! remote coordinates buffer pointer on peer used for NVSHMEM put
            Float3* remoteXPeerPutPtr = nullptr;
            //! remote forces buffer pointer on peer used for NVSHMEM get
            const Float3* remoteForcePeerGetPtr = nullptr;
            //! remote forces buffer pointer on peer used for NVSHMEM put
            Float3* remoteForcePeerPutPtr = nullptr;
            //! atom offset into receiver-rank X buffer used as destination for PUT
            int putAtomOffsetInReceiverRankXBuf_ = 0;
        } nvshmemData;
    };

    struct SharedBuffers
    {
        //! full forces buffer in GPU memory, symmetrically allocated using nvshmem_malloc
        DeviceBuffer<Float3> d_f = nullptr;
        //! full coordinates buffer in GPU memory, symmetrically allocated using nvshmem_malloc
        DeviceBuffer<Float3> d_x = nullptr;
        //! base of receiver-rank X signals
        DeviceBuffer<uint64_t> d_signalReceiverRankX_ = nullptr;
        //! base of sender-rank X signals
        DeviceBuffer<uint64_t> d_signalSenderRankX_ = nullptr;
        //! base of receiver-rank F signals
        DeviceBuffer<uint64_t> d_signalReceiverRankF_ = nullptr;
    };
    SharedBuffers sharedBuffers_;

    // NVSHMEM sync base pointer and per-sync size
    DeviceBuffer<uint64_t> d_ppHaloExSyncBase_     = nullptr;
    int                    ppHaloExPerSyncBufSize_ = 0;

    /*! \brief Copy per-pulse metadata to device. */
    void allocateAndCopyHaloExchangeData();
    /*! \brief Get the event synchronizer for the fused forces ready on device.
     *  \returns  The event to synchronize the stream that consumes forces on device.
     */
    GpuEventSynchronizer* getForcesReadyOnDeviceEvent();
    /*! \brief Allocate unified send/recv buffers and point entries into them.
     * Uses max-reduced per-pulse extents to size symmetric buffers and assigns
     * aligned per-entry pointers into those buffers.
     */
    void allocateUnifiedHaloBuffers();


private:
    // Threads per block used by NVSHMEM kernels
    static constexpr int c_nvshmemThreadsPerBlock = 1024; // fix this for fused kernel
    /*! \brief Minimum Alignment (bytes) for per-entry regions in unified buffers.
     *
     * Chosen as 256 B to be safe and performant:
     * - Matches CUDA global allocation guarantees (cuda device mem allocators returns 256 B-aligned addresses).
     * - Equals PTX max L2 prefetch size 256 B-aligned bases tend to avoid split transactions.
     * - Fits well with NVSHMEM puts/gets and vectorized LD/ST; future-friendly for TMA usage (min 16B but 128B to be performant).
     */
    static constexpr int c_haloEntryAlignBytes = 256;
    //! Device stream for this halo exchange
    std::unique_ptr<DeviceStream> haloStream_;
    //! GPU context object
    const DeviceContext& deviceContext_;
    //! Wallclock cycle accounting
    gmx_wallcycle* wcycle_ = nullptr;
    //! Coordinates communication signal counter, used to synchronize with the send-to(put)
    // and recv-from ranks for the current timestep. Typically the value for the
    // signal counter is the timestep value and is same for all the pulse/dim in same
    // timestep.
    uint64_t signalReceiverRankXCounter_ = 0;
    //! Forces communication signal counter,  used to synchronize with the send-to(put)
    // and recv-from ranks for the current timestep. Typically the value for the
    // signal counter is the timestep value and is same for all the pulse/dim in same
    // timestep.
    uint64_t signalReceiverRankFCounter_ = 0;
    //! Whether to insert a synchronization at the fused force kernel (debug/knob).
    //! Toggle via env var: GMX_NVSHMEM_FUSED_FORCE_KERNEL_SYNC
    bool enableFusedForceKernelSync_ = false;
    //! Total number of pulses across all dimensions
    int totalNumPulses_ = 0;
    //! maximum grid size across all X pulses
    int32_t maxGridXSize_ = 0;
    //! maximum grid size across all F pulses
    int32_t maxGridFSize_ = 0;

    //! per-pulse device-side grid sync counters for F kernel
    DeviceBuffer<uint32_t> d_fGridSync_ = nullptr;
    //! number of F grid sync pulses
    int d_fGridSyncSize_ = -1;
    //! capacity of F grid sync pulses
    int d_fGridSyncSizeAlloc_ = -1;

    //! per-pulse device-side arrive-wait barriers for X kernel
    DeviceBuffer<DeviceBarrier> d_xGridSync_ = nullptr;
    //! number of X grid sync pulses
    int d_xGridSyncSize_ = -1;
    //! capacity of X grid sync pulses
    int d_xGridSyncSizeAlloc_ = -1;

    //! device array of per-pulse metadata
    DeviceBuffer<HaloExchangeData> d_haloExchangeData_ = nullptr;
    //! number of per-pulse metadata records
    int d_haloExchangeDataSize_ = -1;
    //! capacity of per-pulse metadata records
    int d_haloExchangeDataSizeAlloc_ = -1;
    //! host-side vector of per-pulse metadata
    gmx::HostVector<HaloExchangeData> haloExchangeData_;
    //! Event triggered when coordinate halo has been launched
    GpuEventSynchronizer coordinateHaloLaunched_;
    //! Event triggered when force halo has been launched
    GpuEventSynchronizer forceHaloLaunched_;
    // MPI communicator used for symmetric allocations sizing (NVSHMEM path)
    //! MPI communicator used for simulation
    MPI_Comm mpi_comm_mysim_ = MPI_COMM_NULL;
    //! MPI communicator involving PP + PME.
    MPI_Comm mpi_comm_mysim_world_ = MPI_COMM_NULL;

    // Unified send/recv buffers across all dims/pulses (optional fused path)
    //! unified send buffer across all pulses
    DeviceBuffer<Float3> d_unifiedSendBuf_ = nullptr;
    //! all pulses extent in unified send buffer (elements)
    int unifiedSendSize_ = -1;
    //! capacity for all pulses extent in unified send buffer
    int unifiedSendCapacity_ = -1;
    //! unified recv buffer across all pulses
    DeviceBuffer<Float3> d_unifiedRecvBuf_ = nullptr;
    //! all pulses extent in unified recv buffer (elements)
    int unifiedRecvSize_ = -1;
    //! capacity for all pulses extent in unified recv buffer
    int unifiedRecvCapacity_ = -1;

private:
    /*! \brief Backend-specific function for launching coordinate packing and send kernel. */
    void launchPackXKernel(const matrix box);
    /*! \brief Backend-specific function for launching force unpacking and recv kernel. */
    void launchUnpackFKernel(bool accumulateForces);
};

} // namespace gmx

#endif // GMX_DOMDEC_FUSED_GPUHALOEXCHANGE_H
