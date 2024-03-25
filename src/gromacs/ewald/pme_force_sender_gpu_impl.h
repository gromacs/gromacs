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
 * \brief Declaration of class which sends PME Force from GPU memory to PP task
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_ewald
 */
#ifndef GMX_PMEFORCESENDERGPU_IMPL_H
#define GMX_PMEFORCESENDERGPU_IMPL_H

#include <atomic>
#include <new>

#include "gromacs/ewald/pme_force_sender_gpu.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/utility/arrayref.h"

// Portable definition of cache line size
#ifdef __cpp_lib_hardware_interference_size
using std::hardware_destructive_interference_size;
#else
constexpr std::size_t hardware_destructive_interference_size = 64;
#endif

class GpuEventSynchronizer;

namespace gmx
{

/*! \internal \brief Class with interfaces and data for CUDA version of PME Force sending functionality*/

typedef struct CacheLineAlignedFlag
{
    alignas(hardware_destructive_interference_size) bool flag;
} CacheLineAlignedFlag;

/*! \internal
 *  \brief Object to manage communications with a specific PP rank
 */
struct PpForceCommManager
{
    //! Stream used communication with remote PP rank
    std::unique_ptr<DeviceStream> stream;
    //! Event used for manging sync with remote PP rank
    std::unique_ptr<GpuEventSynchronizer> event;
    //! Flag to track when PP transfer event has been recorded
    std::unique_ptr<std::atomic<CacheLineAlignedFlag>> eventRecorded;
    //! Address of local force buffer to send to remote PP rank
    DeviceBuffer<RVec> localForcePtr;
    //! CPU force buffer pointer for remote PP rank
    Float3* pmeRemoteCpuForcePtr;
    //! GPU force buffer pointers for remote PP rank
    Float3* pmeRemoteGpuForcePtr;
};

class PmeForceSenderGpu::Impl
{

public:
    /*! \brief Creates PME GPU Force sender object
     * \param[in] pmeForcesReady  Event synchronizer marked when PME forces are ready on the GPU
     * \param[in] comm            Communicator used for simulation
     * \param[in] deviceContext   GPU context
     * \param[in] ppRanks         List of PP ranks
     */
    Impl(GpuEventSynchronizer*  pmeForcesReady,
         MPI_Comm               comm,
         const DeviceContext&   deviceContext,
         gmx::ArrayRef<PpRanks> ppRanks);
    // NOLINTNEXTLINE(performance-trivially-destructible)
    ~Impl();

    /*! \brief
     * Sets location of force to be sent to each PP rank
     * \param[in] d_f   force buffer in GPU memory
     */
    void setForceSendBuffer(DeviceBuffer<Float3> d_f);

    /*! \brief
     * Send force to PP rank (used with Thread-MPI)
     * \param[in] ppRank                   PP rank to receive data
     * \param[in] numAtoms                 number of atoms to send
     * \param[in] sendForcesDirectToPpGpu  whether forces are transferred direct to remote GPU memory
     */
    void sendFToPpPeerToPeer(int ppRank, int numAtoms, bool sendForcesDirectToPpGpu);

    /*! \brief
     * Send force to PP rank (used with Lib-MPI)
     * \param[in] sendbuf  force buffer in GPU memory
     * \param[in] offset   starting element in buffer
     * \param[in] numBytes number of bytes to transfer
     * \param[in] ppRank   PP rank to receive data
     * \param[in] request  MPI request to track asynchronous MPI call status
     */
    void sendFToPpGpuAwareMpi(DeviceBuffer<RVec> sendbuf, int offset, int numBytes, int ppRank, MPI_Request* request);

    void waitForEvents();

private:
    //! Event indicating when PME forces are ready on the GPU in order for PP stream to sync with the PME stream
    GpuEventSynchronizer* pmeForcesReady_;
    //! communicator for simulation
    MPI_Comm comm_;
    //! list of PP ranks
    gmx::ArrayRef<PpRanks> ppRanks_;
    //! Whether GPU to CPU communication should be staged as GPU to
    //! GPU via P2P cudaMemcpy, then local D2H, for thread-MPI This
    //! may be beneficial when using servers with direct links between
    //! GPUs, but direct communication is expected to be advantageous
    //! for PCIe-only servers or for low atom counts (for which
    //! latency is important).
    bool stageThreadMpiGpuCpuComm_ = false;
    //! Communication manager objects corresponding to multiple receiving PP ranks
    std::vector<PpForceCommManager> ppCommManagers_;
};

} // namespace gmx

#endif
