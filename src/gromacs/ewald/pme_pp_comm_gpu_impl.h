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
 * \brief Declares GPU implementation class for PME-PP communications
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_ewald
 */
#ifndef GMX_PME_PP_COMM_GPU_IMPL_H
#define GMX_PME_PP_COMM_GPU_IMPL_H

#include <atomic>

#include "gromacs/ewald/pme_pp_comm_gpu.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxmpi.h"


namespace gmx
{

/*! \internal \brief Class with interfaces and data for GPU version of PME-PP Communication */
class PmePpCommGpu::Impl
{

public:
    /*! \brief Creates PME-PP GPU communication object.
     *
     * \param[in] comm              Communicator used for simulation
     * \param[in] pmeRank           Rank of PME task
     * \param[in] pmeCpuForceBuffer Buffer for PME force in CPU memory
     * \param[in] deviceContext     GPU context.
     * \param[in] deviceStream      Device stream.
     * \param[in] useNvshmem        NVSHMEM enable/disable for GPU comm.
     */
    Impl(MPI_Comm                    comm,
         int                         pmeRank,
         gmx::HostVector<gmx::RVec>* pmeCpuForceBuffer,
         const DeviceContext&        deviceContext,
         const DeviceStream&         deviceStream,
         bool                        useNvshmem);
    ~Impl();

    /*! \brief Perform steps required when buffer size changes
     * \param[in]  size   Number of elements in buffer
     */
    void reinit(int size);

    /*! \brief Pull force buffer directly from GPU memory on PME
     * rank to either GPU or CPU memory on PP task using CUDA
     * Memory copy or GPU-aware MPI.
     *
     * recvPtr should be in GPU or CPU memory if recvPmeForceToGpu
     * is true or false, respectively. If receiving to GPU, this
     * method should be called before the local GPU buffer
     * operations. If receiving to CPU it should be called
     * before forces are reduced with the other force
     * contributions on the CPU. It will automatically wait for
     * remote PME force data to be ready.
     *
     * \param[out] recvPtr CPU or GPU buffer to receive PME force data
     * \param[in] recvSize Number of elements to receive
     * \param[in] receivePmeForceToGpu Whether receive is to GPU, otherwise CPU
     */
    void receiveForceFromPme(Float3* recvPtr, int recvSize, bool receivePmeForceToGpu);

    /*! \brief Push coordinates buffer directly to GPU memory on PME
     * task, from either GPU or CPU memory on PP task using CUDA
     * Memory copy or GPU-aware MPI. If sending from GPU, this method should
     * be called after the local GPU coordinate buffer operations.
     * The remote PME task will automatically wait for data to be copied
     * before commencing PME force calculations.
     * \param[in] sendPtr Buffer with coordinate data
     * \param[in] sendSize Number of elements to send
     * \param[in] coordinatesReadyOnDeviceEvent Event recorded when coordinates are available on device
     */
    void sendCoordinatesToPme(Float3* sendPtr, int sendSize, GpuEventSynchronizer* coordinatesReadyOnDeviceEvent);

    /*! \brief
     * Return pointer to buffer used for staging PME force on GPU
     */
    DeviceBuffer<Float3> getGpuForceStagingPtr();

    /*! \brief
     * Return pointer to event recorded when forces are ready
     */
    GpuEventSynchronizer* getForcesReadySynchronizer();

    /*! \brief
     * Return pointer to NVSHMEM sync object used for staging PME force on GPU
     */
    DeviceBuffer<uint64_t> getGpuForcesSyncObj();

private:
    /*! \brief Receive buffer from GPU memory on PME rank to either
     * GPU or CPU memory on PP rank. Data is pushed from PME force
     * sender object using CUDA memory copy functionality, and this
     * method performs the necessary synchronization on that
     * communication. This method is used with thread-MPI.
     * \param[in] receivePmeForceToGpu Whether receive is to GPU, otherwise CPU
     */
    void receiveForceFromPmePeerToPeer(bool receivePmeForceToGpu);

    /*! \brief Receive buffer from GPU memory on PME rank to either
     * GPU or CPU memory on PP rank using GPU-aware MPI. This method
     * is used with process-MPI.
     * \param[out] recvPtr CPU or GPU buffer to receive PME force data into
     * \param[in] recvSize Number of elements to receive
     */
    void receiveForceFromPmeGpuAwareMpi(Float3* recvPtr, int recvSize);

    /*! \brief Push coordinates buffer directly to GPU memory on PME
     * task, from either GPU or CPU memory on PP task using CUDA Memory copy.
     * This method is used with Thread-MPI.
     * \param[in] sendPtr Buffer with coordinate data
     * \param[in] sendSize Number of elements to send
     * \param[in] coordinatesReadyOnDeviceEvent Event recorded when coordinates are available on device
     */
    void sendCoordinatesToPmePeerToPeer(Float3*               sendPtr,
                                        int                   sendSize,
                                        GpuEventSynchronizer* coordinatesReadyOnDeviceEvent);

    /*! \brief Push coordinates buffer directly to GPU memory on PME
     * task, from either GPU or CPU memory on PP task using GPU-aware MPI.
     * This method is used with process-MPI.
     * \param[in] sendPtr Buffer with coordinate data
     * \param[in] sendSize Number of elements to send
     * \param[in] coordinatesReadyOnDeviceEvent Event recorded when coordinates are available on device
     */
    void sendCoordinatesToPmeGpuAwareMpi(Float3*               sendPtr,
                                         int                   sendSize,
                                         GpuEventSynchronizer* coordinatesReadyOnDeviceEvent);

    //! Device context handle
    const DeviceContext& deviceContext_;
    //! Handle for the device stream used for the communication operations in this class
    const DeviceStream& pmePpCommStream_;
    //! Remote location of PME coordinate data buffer
    Float3* remotePmeXBuffer_ = nullptr;
    //! communicator for simulation
    MPI_Comm comm_;
    //! Rank of PME task
    int pmeRank_ = -1;
    //! Buffer for PME force on CPU
    gmx::HostVector<gmx::RVec>* pmeCpuForceBuffer_;
    //! Buffer for staging PME force on GPU
    DeviceBuffer<gmx::RVec> d_pmeForces_;
    //! number of atoms in PME force staging array
    int d_pmeForcesSize_ = -1;
    //! number of atoms allocated in recvbuf array
    int d_pmeForcesSizeAlloc_ = -1;
    //! PME force synchronization NVSHMEM object
    DeviceBuffer<uint64_t> forcesReadyNvshmemFlags;
    //! PME force synchronization NVSHMEM object size tracker
    int forcesReadyNvshmemFlagsSize_ = -1;
    //! PME force synchronization NVSHMEM object size tracker
    int forcesReadyNvshmemFlagsSizeAlloc_ = -1;
    //! Event recorded when PME forces are ready on PME task
    GpuEventSynchronizer forcesReadySynchronizer_;
    //! Event recorded when coordinates have been transferred to PME task
    GpuEventSynchronizer pmeCoordinatesSynchronizer_;
    //! Event recorded by remote PME task when forces have been transferred
    GpuEventSynchronizer* remotePmeForceSendEvent_;
    //! Flag to track when remote PP event has been recorded, ready for enqueueing
    volatile std::atomic<bool>* remotePmeForceSendEventRecorded_;
    //! Whether GPU to CPU communication should staged through GPU
    //! memory rather than performed directly, for lib-MPI. Staging is
    //! expected to have significant benefits for systems with servers
    //! with direct links between GPUs, because it allows the device
    //! to host transfer to be split across multiple PCIe buses, thus
    //! accessing more bandwidth. Direct communication may have
    //! benefits on servers with only PCIe connectivity, and/or for
    //! small atom counts where latency is more important than
    //! bandwidth.
    bool stageLibMpiGpuCpuComm_ = true;
    // MPI Request associated with non-blocking coordinate send
    MPI_Request coordinateSendRequest_;
    // Flag on whether a non-blocking coordinate send is active
    bool coordinateSendRequestIsActive_ = false;
    // Flag on whether to use NVSHMEM for GPU communication
    bool useNvshmem_ = false;
};

} // namespace gmx

#endif
