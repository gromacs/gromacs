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
 * \brief Implememnts backend-specific code for PME-PP communication using CUDA.
 *
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/utility/gmxmpi.h"

#include "pme_force_sender_gpu_impl.h"

namespace gmx
{

/*! \brief Send PME synchronizer directly using CUDA memory copy */
void PmeForceSenderGpu::Impl::sendFToPpPeerToPeer(int ppRank, int numAtoms, bool sendForcesDirectToPpGpu)
{

    GMX_ASSERT(GMX_THREAD_MPI, "sendFToPpCudaDirect is expected to be called only for Thread-MPI");

#if GMX_MPI
    Float3* pmeRemoteForcePtr = (sendForcesDirectToPpGpu || stageThreadMpiGpuCpuComm_)
                                        ? ppCommManagers_[ppRank].pmeRemoteGpuForcePtr
                                        : ppCommManagers_[ppRank].pmeRemoteCpuForcePtr;

    pmeForcesReady_->enqueueWaitEvent(*ppCommManagers_[ppRank].stream);

    // Push data to remote GPU's memory
    cudaError_t stat = cudaMemcpyAsync(asFloat3(pmeRemoteForcePtr),
                                       ppCommManagers_[ppRank].localForcePtr,
                                       numAtoms * sizeof(rvec),
                                       cudaMemcpyDefault,
                                       ppCommManagers_[ppRank].stream->stream());
    CU_RET_ERR(stat, "cudaMemcpyAsync on Recv from PME CUDA direct data transfer failed");

    if (stageThreadMpiGpuCpuComm_ && !sendForcesDirectToPpGpu)
    {
        // Perform local D2H (from remote GPU memory to remote PP rank's CPU memory)
        // to finalize staged data transfer
        stat = cudaMemcpyAsync(ppCommManagers_[ppRank].pmeRemoteCpuForcePtr,
                               ppCommManagers_[ppRank].pmeRemoteGpuForcePtr,
                               numAtoms * sizeof(rvec),
                               cudaMemcpyDefault,
                               ppCommManagers_[ppRank].stream->stream());
        CU_RET_ERR(stat, "cudaMemcpyAsync on local device to host transfer of PME forces failed");
    }

    ppCommManagers_[ppRank].event->markEvent(*ppCommManagers_[ppRank].stream);
    std::atomic<bool>* tmpPpCommEventRecordedPtr =
            reinterpret_cast<std::atomic<bool>*>(ppCommManagers_[ppRank].eventRecorded.get());
    tmpPpCommEventRecordedPtr->store(true, std::memory_order_release);
#else
    GMX_UNUSED_VALUE(ppRank);
    GMX_UNUSED_VALUE(numAtoms);
#endif
}

} // namespace gmx
