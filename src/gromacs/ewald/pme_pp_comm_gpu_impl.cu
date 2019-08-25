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
 * \brief Implements PME-PP communication using CUDA
 *
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "pme_pp_comm_gpu_impl.h"

#include "config.h"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/gpueventsynchronizer.cuh"
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{

PmePpCommGpu::Impl::Impl(MPI_Comm comm, int pmeRank)
    : comm_(comm),
      pmeRank_(pmeRank)
{
    GMX_RELEASE_ASSERT(GMX_THREAD_MPI, "PME-PP GPU Communication is currently only supported with thread-MPI enabled");
    cudaStreamCreate(&pmePpCommStream_);
}

PmePpCommGpu::Impl::~Impl() = default;

void PmePpCommGpu::Impl::receiveForceBufferAddress()
{
    // This rank will pull data from the PME rank, so needs to recieve the remote PME buffer address.
    MPI_Recv(&remotePmeFBuffer_, sizeof(void**), MPI_BYTE, pmeRank_,
             0, comm_, MPI_STATUS_IGNORE);
    return;
}

void PmePpCommGpu::Impl::receiveForceFromPmeCudaDirect(void *recvPtr, int recvSize)
{

    // Receive event from PME task and add to stream, to ensure pull of data doesn't
    // occur before PME force calc is completed
    GpuEventSynchronizer *pmeSync;
    MPI_Recv(&pmeSync, sizeof(GpuEventSynchronizer*),
             MPI_BYTE, pmeRank_, 0,
             comm_, MPI_STATUS_IGNORE);
    pmeSync->enqueueWaitEvent(pmePpCommStream_);

    // Pull force data from remote GPU
    cudaError_t stat = cudaMemcpyAsync(recvPtr, remotePmeFBuffer_,
                                       recvSize*3*sizeof(float), cudaMemcpyDefault,
                                       pmePpCommStream_);
    CU_RET_ERR(stat, "cudaMemcpyAsync on Recv from PME CUDA direct data transfer failed");

    // Ensure CPU waits for PME forces to be copied before reducing
    // them with other forces on the CPU
    cudaStreamSynchronize(pmePpCommStream_);
}

PmePpCommGpu::PmePpCommGpu(MPI_Comm comm, int pmeRank)
    : impl_(new Impl(comm,  pmeRank))
{
}

PmePpCommGpu::~PmePpCommGpu() = default;

void PmePpCommGpu::receiveForceBufferAddress()
{
    impl_->receiveForceBufferAddress();
}

void PmePpCommGpu::receiveForceFromPmeCudaDirect(void *recvPtr, int recvSize)
{
    impl_->receiveForceFromPmeCudaDirect(recvPtr, recvSize);
}

} //namespace gmx
