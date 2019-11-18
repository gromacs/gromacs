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
 * \brief Implements class which recieves coordinates to GPU memory on PME task using CUDA
 *
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "pme_coordinate_receiver_gpu_impl.h"

#include "config.h"

#include <assert.h>
#include <stdio.h>

#include "gromacs/ewald/pme.h"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/gpueventsynchronizer.cuh"
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{

PmeCoordinateReceiverGpu::Impl::Impl(void* pmeStream, MPI_Comm comm, gmx::ArrayRef<PpRanks> ppRanks) :
    pmeStream_(*static_cast<cudaStream_t*>(pmeStream)),
    comm_(comm),
    ppRanks_(ppRanks)
{
    GMX_RELEASE_ASSERT(
            GMX_THREAD_MPI,
            "PME-PP GPU Communication is currently only supported with thread-MPI enabled");
}

PmeCoordinateReceiverGpu::Impl::~Impl() = default;

void PmeCoordinateReceiverGpu::Impl::sendCoordinateBufferAddressToPpRanks(rvec* d_x)
{

    int ind_start = 0;
    int ind_end   = 0;
    for (const auto& receiver : ppRanks_)
    {
        ind_start = ind_end;
        ind_end   = ind_start + receiver.numAtoms;

        // Data will be transferred directly from GPU.
        void* sendBuf = reinterpret_cast<void*>(&d_x[ind_start]);

#if GMX_MPI
        MPI_Send(&sendBuf, sizeof(void**), MPI_BYTE, receiver.rankId, 0, comm_);
#endif
    }
}

/*! \brief Receive coordinate data directly using CUDA memory copy */
void PmeCoordinateReceiverGpu::Impl::receiveCoordinatesFromPpCudaDirect(int ppRank)
{
    // Data will be pushed directly from PP task

#if GMX_MPI
    // Receive event from PP task and add to PME stream, to ensure PME calculation doesn't
    // commence until coordinate data has been transferred
    GpuEventSynchronizer* ppSync;
    MPI_Recv(&ppSync, sizeof(GpuEventSynchronizer*), MPI_BYTE, ppRank, 0, comm_, MPI_STATUS_IGNORE);
    ppSync->enqueueWaitEvent(pmeStream_);
#endif
}

PmeCoordinateReceiverGpu::PmeCoordinateReceiverGpu(void*                  pmeStream,
                                                   MPI_Comm               comm,
                                                   gmx::ArrayRef<PpRanks> ppRanks) :
    impl_(new Impl(pmeStream, comm, ppRanks))
{
}

PmeCoordinateReceiverGpu::~PmeCoordinateReceiverGpu() = default;

void PmeCoordinateReceiverGpu::sendCoordinateBufferAddressToPpRanks(rvec* d_x)
{
    impl_->sendCoordinateBufferAddressToPpRanks(d_x);
}

void PmeCoordinateReceiverGpu::receiveCoordinatesFromPpCudaDirect(int ppRank)
{
    impl_->receiveCoordinatesFromPpCudaDirect(ppRank);
}

} // namespace gmx
