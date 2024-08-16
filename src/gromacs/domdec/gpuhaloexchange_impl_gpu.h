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
 * \brief Declares the implementation of GPU Halo Exchange.
 *
 * This header file is needed to include from both the device-side
 * kernels file, and the host-side management code.
 *
 * \author Alan Gray <alang@nvidia.com>
 * \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDEC_GPUHALOEXCHANGE_IMPL_GPU_H
#define GMX_DOMDEC_GPUHALOEXCHANGE_IMPL_GPU_H

#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/utility/gmxmpi.h"

#include "domdec_internal.h"

#if GMX_NVSHMEM
#    include <cuda/barrier>
#endif

#if GMX_GPU_HIP
#    include "gromacs/gpu_utils/gputraits_hip.h"
#endif

#if GMX_GPU_SYCL
#    include "gromacs/gpu_utils/gputraits_sycl.h"
#endif

struct gmx_wallcycle;

namespace gmx
{

/*! \brief Whether the halo exchange is of coordinates or forces. */
enum class HaloType
{
    Coordinates,
    Forces,
};

/*! \internal \brief Struct encapsulating data for NVSHMEM Enabled GPU Halo Exchange */
struct nvshmemHaloExchangeOps
{
    //! offset of remote PP rank coord buffer used for nvshmem_put
    int putAtomOffsetInReceiverRankXBuf_ = 0;
    //! signal object used for notifying rank which will produce/put the halo coordinates
    // about readiness of the coordinates buffer.
    DeviceBuffer<uint64_t> d_signalSenderRankX_ = nullptr;
    //! signal object used for notifying rank which will consume the halo forces
    // about completion of the puts.
    DeviceBuffer<uint64_t> d_signalReceiverRankF_ = nullptr;
    //! Forces communication signal counter,  used to synchronize with the send-to(put)
    // and recv-from ranks for the current timestep. Typically the value for the
    // signal counter is the timestep value and is same for all the pulse/dim in same
    // timestep.
    uint64_t signalReceiverRankFCounter_ = 0;
    //! signal object used for notifying rank which will consume the halo coordinates
    // about completion of the puts.
    DeviceBuffer<uint64_t> d_signalReceiverRankX_ = nullptr;
    //! Coordinates communication signal counter, used to synchronize with the send-to(put)
    // and recv-from ranks for the current timestep. Typically the value for the
    // signal counter is the timestep value and is same for all the pulse/dim in same
    // timestep.
    uint64_t signalReceiverRankXCounter_ = 0;
    //! value to initialize the arrive-wait barrier, this is grid size of X kernel.
    uint32_t arriveWaitBarrierVal_ = 0;
    //! Grid size of X kernel.
    uint32_t gridDimX_ = 0;
#if GMX_NVSHMEM
    //! device scoped arrive-wait barrier to sync threadblocks for signalling writes
    cuda::barrier<cuda::thread_scope_device>* d_arriveWaitBarrier_ = nullptr;
#endif
    //! offset from which signal object for the given pulse/dim should be used.
    int signalObjOffset_ = 0;
    //! Thread block size for NVSHMEM based X and F kernels.
    const static int c_nvshmemThreadsPerBlock = 256 * 4;
};

/*! \internal \brief Class with interfaces and data for GPU Halo Exchange */
class GpuHaloExchange::Impl
{

public:
    /*! \brief Creates GPU Halo Exchange object.
     *
     * \param [inout] dd                       domdec structure
     * \param [in]    dimIndex                 the dimension index for this instance
     * \param [in]    mpi_comm_mysim           communicator used for simulation
     * \param [in]    mpi_comm_mysim_world     communicator used for simulation with PP + PME.
     * \param [in]    deviceContext            GPU device context
     * \param [in]    pulse                    the communication pulse for this instance
     * \param [in]    useNvshmem               use NVSHMEM for communication
     * \param [in]    wcycle                   The wallclock counter
     */
    Impl(gmx_domdec_t*        dd,
         int                  dimIndex,
         MPI_Comm             mpi_comm_mysim,
         MPI_Comm             mpi_comm_mysim_world,
         const DeviceContext& deviceContext,
         int                  pulse,
         bool                 useNvshmem,
         gmx_wallcycle*       wcycle);
    ~Impl();

    /*! \brief
     * (Re-) Initialization for GPU halo exchange
     * \param [in] d_coordinatesBuffer  pointer to coordinates buffer in GPU memory
     * \param [in] d_forcesBuffer       pointer to forces buffer in GPU memory
     */
    void reinitHalo(DeviceBuffer<Float3> d_coordinatesBuffer, DeviceBuffer<Float3> d_forcesBuffer);

    /*! \brief
     * (Re-) Initialization for NVSHMEM Signal objects
     * \param [in] cr  Communication structure ref.
     * \param [in] signalObjOffset  offset of the signal object corresponding to given pulse/dim.
     */
    void reinitNvshmemSignal(const t_commrec& cr, int signalObjOffset);

    /*! \brief
     * GPU halo exchange of coordinates buffer
     * \param [in] box  Coordinate box (from which shifts will be constructed)
     * \param [in] dependencyEvent   Dependency event for this operation
     * \returns                      Event recorded when this operation has been launched
     */
    GpuEventSynchronizer* communicateHaloCoordinates(const matrix box, GpuEventSynchronizer* dependencyEvent);

    /*! \brief  GPU halo exchange of force buffer
     * \param [in] accumulateForces  True if forces should accumulate, otherwise they are set
     * \param [in] dependencyEvents  Dependency events for this operation
     */
    void communicateHaloForces(bool                                           accumulateForces,
                               FixedCapacityVector<GpuEventSynchronizer*, 2>* dependencyEvents);

    /*! \brief Get the event synchronizer for the forces ready on device.
     *  \returns  The event to synchronize the stream that consumes forces on device.
     */
    GpuEventSynchronizer* getForcesReadyOnDeviceEvent();

    /*! \brief Destructor for symmetric d_recvBuf used by NVSHMEM.
     */
    void destroyGpuHaloExchangeNvshmemBuf();

private:
    /*! \brief Data transfer wrapper for GPU halo exchange
     * \param [in] sendPtr      send buffer address
     * \param [in] sendSize     number of elements to send
     * \param [in] sendRank     rank of destination
     * \param [in] recvPtr      receive buffer address
     * \param [in] recvSize     number of elements to receive
     * \param [in] recvRank     rank of source
     * \param [in] haloType     whether halo exchange is of coordinates or forces
     */
    void communicateHaloData(Float3*  sendPtr,
                             int      sendSize,
                             int      sendRank,
                             Float3*  recvPtr,
                             int      recvSize,
                             int      recvRank,
                             HaloType haloType);

    /*! \brief Data transfer for GPU halo exchange using peer-to-peer copies
     * \param [inout] sendPtr    address to send data from
     * \param [in] sendSize      number of atoms to be sent
     * \param [in] sendRank      rank to send data to
     * \param [in] remotePtr     remote address to recv data
     * \param [in] recvRank      rank to recv data from
     * \param [in] haloType      whether halo exchange is of coordinates or forces
     */
    void communicateHaloDataPeerToPeer(Float3*  sendPtr,
                                       int      sendSize,
                                       int      sendRank,
                                       Float3*  remotePtr,
                                       int      recvRank,
                                       HaloType haloType);

    /*! \brief Data transfer for GPU halo exchange using GPU-aware MPI
     * \param [in] sendPtr      send buffer address
     * \param [in] sendSize     number of elements to send
     * \param [in] sendRank     rank of destination
     * \param [in] recvPtr      receive buffer address
     * \param [in] recvSize     number of elements to receive
     * \param [in] recvRank     rank of source
     */
    void communicateHaloDataGpuAwareMpi(Float3* sendPtr,
                                        int     sendSize,
                                        int     sendRank,
                                        Float3* recvPtr,
                                        int     recvSize,
                                        int     recvRank);

    //! Backend-specific function for launching coordinate packing kernel
    void launchPackXKernel(const matrix box);

    //! Backend-specific function for launching force unpacking kernel
    void launchUnpackFKernel(bool accumulateForces);

    /*! \brief Transfer GPU coordinate halo data out of place
     * \param [in] sendPtr      send buffer address
     * \param [in] sendSize     number of elements to send
     * \param [in] sendRank     rank of destination
     * \param [in] recvSize     number of elements to receive
     * \param [in] recvRank     rank of source
     */
    void communicateHaloCoordinatesOutOfPlace(DeviceBuffer<Float3> sendPtr,
                                              int                  sendSize,
                                              int                  sendRank,
                                              int                  recvSize,
                                              int                  recvRank);

    /*! \brief Transfer GPU force halo data out of place
     * \param [in] sendPtr      send buffer address
     * \param [in] sendSize     number of elements to send
     * \param [in] sendRank     rank of destination
     * \param [in] recvSize     number of elements to receive
     * \param [in] recvRank     rank of source
     */
    void communicateHaloForcesOutOfPlace(DeviceBuffer<Float3> sendPtr,
                                         int                  sendSize,
                                         int                  sendRank,
                                         int                  recvSize,
                                         int                  recvRank);

    /*! \brief Exchange coordinate-ready event with neighbor ranks and enqueue wait in halo stream
     * \param [in] coordinatesReadyOnDeviceEvent event recorded when coordinates/forces are ready to device.
     */
    void enqueueWaitRemoteCoordinatesReadyEvent(GpuEventSynchronizer* coordinatesReadyOnDeviceEvent);

    //! NVSHMEM-specific function for setting the grid size and
    // initialize the arrive wait barrier for X kernel.
    void reinitXGridSizeAndDevBarrier();

    //! Domain decomposition object
    gmx_domdec_t* dd_ = nullptr;
    //! map of indices to be sent from this rank
    gmx::HostVector<int> h_indexMap_;
    //! device copy of index map
    DeviceBuffer<int> d_indexMap_ = nullptr;
    //! number of elements in index map array
    int indexMapSize_ = -1;
    //! number of elements allocated in index map array
    int indexMapSizeAlloc_ = -1;
    //! device buffer for sending packed data
    DeviceBuffer<Float3> d_sendBuf_ = nullptr;
    //! number of atoms in \c sendbuf array
    int sendBufSize_ = -1;
    //! number of atoms allocated in \c sendbuf array
    int sendBufSizeAlloc_ = -1;
    //! device buffer for receiving packed data
    DeviceBuffer<Float3> d_recvBuf_ = nullptr;
    //! maximum size of packed buffer
    int maxPackedBufferSize_ = 0;
    //! number of atoms in \c recvbuf array
    int recvBufSize_ = -1;
    //! number of atoms allocated in \c recvbuf array
    int recvBufSizeAlloc_ = -1;
    //! rank to send data to for X
    int sendRankX_ = 0;
    //! rank to recv data from for X
    int recvRankX_ = 0;
    //! rank to send data to for F
    int sendRankF_ = 0;
    //! rank to recv data from for F
    int recvRankF_ = 0;
    //! send copy size from this rank for X
    int xSendSize_ = 0;
    //! recv copy size to this rank for X
    int xRecvSize_ = 0;
    //! send copy size from this rank for F
    int fSendSize_ = 0;
    //! recv copy size to this rank for F
    int fRecvSize_ = 0;
    //! number of home atoms - offset of local halo region
    int numHomeAtoms_ = 0;
    //! remote GPU coordinates buffer pointer for pushing data
    DeviceBuffer<Float3> remoteXPtr_ = nullptr;
    //! remote GPU force buffer pointer for pushing data
    DeviceBuffer<Float3> remoteFPtr_ = nullptr;
    //! Periodic Boundary Conditions for this rank
    bool usePBC_ = false;
    //! force shift buffer on device
    DeviceBuffer<Float3> d_fShift_ = nullptr;
    //! Event triggered when X halo transfer has been launched with peer-to-peer memory copy
    std::unique_ptr<GpuEventSynchronizer> haloXDataTransferLaunched_;
    //! Event triggered when F halo transfer has been launched with peer-to-peer memory copy
    std::unique_ptr<GpuEventSynchronizer> haloFDataTransferLaunched_;
    //! MPI communicator used for simulation
    MPI_Comm mpi_comm_mysim_;
    //! MPI communicator involving PP + PME.
    MPI_Comm mpi_comm_mysim_world_;
    //! GPU context object
    const DeviceContext& deviceContext_;
    //! Device stream for this halo exchange
    std::unique_ptr<DeviceStream> haloStream_;
    //! full coordinates buffer in GPU memory
    DeviceBuffer<Float3> d_x_ = nullptr;
    //! full forces buffer in GPU memory
    DeviceBuffer<Float3> d_f_ = nullptr;
    //! An event recorded once the exchanged forces are ready on the GPU
    GpuEventSynchronizer fReadyOnDevice_;
    //! The dimension index corresponding to this halo exchange instance
    int dimIndex_ = 0;
    //! Number of halo zones associated with this halo exchange instance
    int numZone_ = 0;
    //! The pulse corresponding to this halo exchange instance
    int pulse_ = 0;
    //! The wallclock counter
    gmx_wallcycle* wcycle_ = nullptr;
    //! The atom offset for receive (x) or send (f) for dimension index and pulse corresponding to this halo exchange instance
    int atomOffset_ = 0;
    //! whether nvshmem should be used.
    bool useNvshmem_ = false;
    // Contains all the objects and metadata required for NVSHMEM enabled
    // PP Halo exchange.
    nvshmemHaloExchangeOps nvshmemHaloExchange_;
    //! Event triggered when coordinate halo has been launched
    GpuEventSynchronizer coordinateHaloLaunched_;
    //! flag on whether the recieve for this halo exchange is performed in-place
    bool receiveInPlace_ = true;
    //! The indices to communicate for this halo exchange
    const gmx_domdec_ind_t* ind_ = nullptr;
    //! Buffer to stage out-of-place send data
    gmx::HostVector<Float3> h_outOfPlaceSendBuffer_;
    //! Buffer to stage out-of-place receive data
    gmx::HostVector<Float3> h_outOfPlaceRecvBuffer_;
};

} // namespace gmx

#endif // GMX_DOMDEC_GPUHALOEXCHANGE_IMPL_GPU_H
