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
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.cuh"
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{

PmePpCommGpu::Impl::Impl(MPI_Comm comm, int pmeRank) : comm_(comm), pmeRank_(pmeRank)
{
    GMX_RELEASE_ASSERT(
            GMX_THREAD_MPI,
            "PME-PP GPU Communication is currently only supported with thread-MPI enabled");
    cudaStreamCreate(&pmePpCommStream_);
}

PmePpCommGpu::Impl::~Impl() = default;

void PmePpCommGpu::Impl::reinit(int size)
{
    // This rank will access PME rank memory directly, so needs to receive the remote PME buffer addresses.
#if GMX_MPI
    MPI_Recv(&remotePmeXBuffer_, sizeof(void**), MPI_BYTE, pmeRank_, 0, comm_, MPI_STATUS_IGNORE);
    MPI_Recv(&remotePmeFBuffer_, sizeof(void**), MPI_BYTE, pmeRank_, 0, comm_, MPI_STATUS_IGNORE);

    // Reallocate buffer used for staging PME force on GPU
    reallocateDeviceBuffer(&d_pmeForces_, size, &d_pmeForcesSize_, &d_pmeForcesSizeAlloc_, nullptr);
#else
    GMX_UNUSED_VALUE(size);
#endif
    return;
}

// TODO make this asynchronous by splitting into this into
// launchRecvForceFromPmeCudaDirect() and sycnRecvForceFromPmeCudaDirect()
void PmePpCommGpu::Impl::receiveForceFromPmeCudaDirect(void* recvPtr, int recvSize, bool receivePmeForceToGpu)
{
#if GMX_MPI
    // Receive event from PME task and add to stream, to ensure pull of data doesn't
    // occur before PME force calc is completed
    GpuEventSynchronizer* pmeSync;
    MPI_Recv(&pmeSync, sizeof(GpuEventSynchronizer*), MPI_BYTE, pmeRank_, 0, comm_, MPI_STATUS_IGNORE);
    pmeSync->enqueueWaitEvent(pmePpCommStream_);

    // Pull force data from remote GPU
    void*       pmeForcePtr = receivePmeForceToGpu ? static_cast<void*>(d_pmeForces_) : recvPtr;
    cudaError_t stat = cudaMemcpyAsync(pmeForcePtr, remotePmeFBuffer_, recvSize * DIM * sizeof(float),
                                       cudaMemcpyDefault, pmePpCommStream_);
    CU_RET_ERR(stat, "cudaMemcpyAsync on Recv from PME CUDA direct data transfer failed");

    if (receivePmeForceToGpu)
    {
        // Record event to be enqueued in the GPU local buffer operations, to
        // satisfy dependency on receiving the PME force data before
        // reducing it with the other force contributions.
        forcesReadySynchronizer_.markEvent(pmePpCommStream_);
    }
    else
    {
        // Ensure CPU waits for PME forces to be copied before reducing
        // them with other forces on the CPU
        cudaStreamSynchronize(pmePpCommStream_);
    }
#else
    GMX_UNUSED_VALUE(recvPtr);
    GMX_UNUSED_VALUE(recvSize);
    GMX_UNUSED_VALUE(receivePmeForceToGpu);
#endif
}

void PmePpCommGpu::Impl::sendCoordinatesToPmeCudaDirect(void* sendPtr,
                                                        int   sendSize,
                                                        bool gmx_unused sendPmeCoordinatesFromGpu,
                                                        GpuEventSynchronizer* coordinatesReadyOnDeviceEvent)
{
#if GMX_MPI
    // ensure stream waits until coordinate data is available on device
    coordinatesReadyOnDeviceEvent->enqueueWaitEvent(pmePpCommStream_);

    cudaError_t stat = cudaMemcpyAsync(remotePmeXBuffer_, sendPtr, sendSize * DIM * sizeof(float),
                                       cudaMemcpyDefault, pmePpCommStream_);
    CU_RET_ERR(stat, "cudaMemcpyAsync on Send to PME CUDA direct data transfer failed");

    // Record and send event to allow PME task to sync to above transfer before commencing force calculations
    pmeCoordinatesSynchronizer_.markEvent(pmePpCommStream_);
    GpuEventSynchronizer* pmeSync = &pmeCoordinatesSynchronizer_;
    MPI_Send(&pmeSync, sizeof(GpuEventSynchronizer*), MPI_BYTE, pmeRank_, 0, comm_);
#else
    GMX_UNUSED_VALUE(sendPtr);
    GMX_UNUSED_VALUE(sendSize);
    GMX_UNUSED_VALUE(sendPmeCoordinatesFromGpu);
    GMX_UNUSED_VALUE(coordinatesReadyOnDeviceEvent);
#endif
}
void* PmePpCommGpu::Impl::getGpuForceStagingPtr()
{
    return static_cast<void*>(d_pmeForces_);
}

void* PmePpCommGpu::Impl::getForcesReadySynchronizer()
{
    return static_cast<void*>(&forcesReadySynchronizer_);
}

PmePpCommGpu::PmePpCommGpu(MPI_Comm comm, int pmeRank) : impl_(new Impl(comm, pmeRank)) {}

PmePpCommGpu::~PmePpCommGpu() = default;

void PmePpCommGpu::reinit(int size)
{
    impl_->reinit(size);
}

void PmePpCommGpu::receiveForceFromPmeCudaDirect(void* recvPtr, int recvSize, bool receivePmeForceToGpu)
{
    impl_->receiveForceFromPmeCudaDirect(recvPtr, recvSize, receivePmeForceToGpu);
}

void PmePpCommGpu::sendCoordinatesToPmeCudaDirect(void*                 sendPtr,
                                                  int                   sendSize,
                                                  bool                  sendPmeCoordinatesFromGpu,
                                                  GpuEventSynchronizer* coordinatesReadyOnDeviceEvent)
{
    impl_->sendCoordinatesToPmeCudaDirect(sendPtr, sendSize, sendPmeCoordinatesFromGpu,
                                          coordinatesReadyOnDeviceEvent);
}

void* PmePpCommGpu::getGpuForceStagingPtr()
{
    return impl_->getGpuForceStagingPtr();
}

void* PmePpCommGpu::getForcesReadySynchronizer()
{
    return impl_->getForcesReadySynchronizer();
}

} // namespace gmx
