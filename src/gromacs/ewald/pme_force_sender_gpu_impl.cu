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
 * \brief Implements PME-PP communication using CUDA
 *
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "pme_force_sender_gpu_impl.h"

#include "config.h"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/gpueventsynchronizer.cuh"
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{

/*! \brief Create PME-PP GPU communication object */
PmeForceSenderGpu::Impl::Impl(GpuEventSynchronizer*  pmeForcesReady,
                              MPI_Comm               comm,
                              gmx::ArrayRef<PpRanks> ppRanks) :
    pmeForcesReady_(pmeForcesReady), comm_(comm), ppRanks_(ppRanks)
{
}

PmeForceSenderGpu::Impl::~Impl() = default;

/*! \brief  sends force buffer address to PP ranks */
void PmeForceSenderGpu::Impl::sendForceBufferAddressToPpRanks(DeviceBuffer<Float3> d_f)
{
    // Need to send address to PP rank only for thread-MPI as PP rank pulls
    // data using cudamemcpy
    if (!GMX_THREAD_MPI)
    {
        return;
    }
#if GMX_MPI
    int ind_start = 0;
    int ind_end   = 0;
    for (const auto& receiver : ppRanks_)
    {
        ind_start = ind_end;
        ind_end   = ind_start + receiver.numAtoms;

        // Data will be transferred directly from GPU.
        Float3* sendBuf = &d_f[ind_start];

        MPI_Send(&sendBuf, sizeof(Float3*), MPI_BYTE, receiver.rankId, 0, comm_);
    }
#else
    GMX_UNUSED_VALUE(d_f);
#endif
}

/*! \brief Send PME synchronizer directly using CUDA memory copy */
void PmeForceSenderGpu::Impl::sendFSynchronizerToPpCudaDirect(int ppRank)
{
    GMX_ASSERT(GMX_THREAD_MPI,
               "sendFSynchronizerToPpCudaDirect is expected to be called only for Thread-MPI");

    // Data will be pulled directly from PP task
#if GMX_MPI
    // TODO Using MPI_Isend would be more efficient, particularly when
    // sending to multiple PP ranks
    MPI_Send(&pmeForcesReady_, sizeof(GpuEventSynchronizer*), MPI_BYTE, ppRank, 0, comm_);
#else
    GMX_UNUSED_VALUE(ppRank);
#endif
}

/*! \brief Send PME data directly using CUDA-aware MPI */
void PmeForceSenderGpu::Impl::sendFToPpCudaMpi(DeviceBuffer<RVec> sendbuf,
                                               int                offset,
                                               int                numBytes,
                                               int                ppRank,
                                               MPI_Request*       request)
{
    GMX_ASSERT(GMX_LIB_MPI, "sendFToPpCudaMpi is expected to be called only for Lib-MPI");

#if GMX_MPI
    // if using GPU direct comm with CUDA-aware MPI, make sure forces are ready on device
    // before sending it to PP ranks
    pmeForcesReady_->waitForEvent();

    MPI_Isend(sendbuf[offset], numBytes, MPI_BYTE, ppRank, 0, comm_, request);

#else
    GMX_UNUSED_VALUE(sendbuf);
    GMX_UNUSED_VALUE(offset);
    GMX_UNUSED_VALUE(numBytes);
    GMX_UNUSED_VALUE(ppRank);
    GMX_UNUSED_VALUE(request);
#endif
}

PmeForceSenderGpu::PmeForceSenderGpu(GpuEventSynchronizer*  pmeForcesReady,
                                     MPI_Comm               comm,
                                     gmx::ArrayRef<PpRanks> ppRanks) :
    impl_(new Impl(pmeForcesReady, comm, ppRanks))
{
}

PmeForceSenderGpu::~PmeForceSenderGpu() = default;

void PmeForceSenderGpu::sendForceBufferAddressToPpRanks(DeviceBuffer<RVec> d_f)
{
    impl_->sendForceBufferAddressToPpRanks(d_f);
}

void PmeForceSenderGpu::sendFSynchronizerToPpCudaDirect(int ppRank)
{
    impl_->sendFSynchronizerToPpCudaDirect(ppRank);
}

void PmeForceSenderGpu::sendFToPpCudaMpi(DeviceBuffer<RVec> sendbuf,
                                         int                offset,
                                         int                numBytes,
                                         int                ppRank,
                                         MPI_Request*       request)
{
    impl_->sendFToPpCudaMpi(sendbuf, offset, numBytes, ppRank, request);
}


} // namespace gmx
