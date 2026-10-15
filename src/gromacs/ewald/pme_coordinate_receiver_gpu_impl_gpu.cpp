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
 * \brief Implements class which receives coordinates to GPU memory on PME task using CUDA/SYCL.
 *
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/ewald/pme_force_sender_gpu.h"
#include "gromacs/ewald/pme_pp_communication.h"
#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/utility/gmxmpi.h"

#include "pme_coordinate_receiver_gpu_impl.h"

namespace gmx
{

PmeCoordinateReceiverGpu::Impl::Impl(MPI_Comm                     comm,
                                     const DeviceContext&         deviceContext,
                                     gmx::ArrayRef<const PpRanks> ppRanks) :
    comm_(comm),
#if GMX_MPI
    requests_(ppRanks.size(), MPI_REQUEST_NULL)
#else
    requests_()
#endif
{
    // Create streams to manage pipelining
    ppCommManagers_.reserve(ppRanks.size());
    for (const auto& ppRank : ppRanks)
    {
        ppCommManagers_.emplace_back(PpCommManager{
                ppRank,
                std::make_unique<DeviceStream>(deviceContext, DeviceStreamPriority::High, false),
                nullptr,
                std::make_unique<GpuEventSynchronizer>(),
                { 0, 0 } });
    }
#if !GMX_MPI
    GMX_UNUSED_VALUE(comm_);
#endif
}

PmeCoordinateReceiverGpu::Impl::~Impl() = default;

void PmeCoordinateReceiverGpu::Impl::reinitCoordinateReceiver(DeviceBuffer<RVec> d_x)
{
#if GMX_MPI
    int indEnd = 0;
    for (auto& ppCommManager : ppCommManagers_)
    {
        int indStart = indEnd;
        indEnd       = indStart + ppCommManager.ppRank.numAtoms;

        ppCommManager.atomRange = std::make_tuple(indStart, indEnd);

        // Need to send address to PP rank only for thread-MPI as PP rank pushes data using cudamemcpy
        if (GMX_THREAD_MPI)
        {
            GMX_RELEASE_ASSERT(GpuConfigurationCapabilities::PpPmeDirectComm,
                               "Direct PME-PP communication with threadMPI needs to be supported "
                               "by the backend.");
            // Data will be transferred directly from GPU.
            void* sendBuf = reinterpret_cast<void*>(asRawDevicePointer(d_x) + indStart);
            MPI_Send(&sendBuf,
                     sizeof(void**),
                     MPI_BYTE,
                     ppCommManager.ppRank.rankId,
                     eCommType_COORD_GPU_REMOTE_GPU_PTR,
                     comm_);
        }
    }
#else
    GMX_UNUSED_VALUE(d_x);
#endif
}

void PmeCoordinateReceiverGpu::Impl::launchReceiveCoordinatesFromPp(DeviceBuffer<RVec> recvbuf,
                                                                    int                atomOffset,
                                                                    int                numBytes,
                                                                    int                senderIndex)
{
#if GMX_MPI
    GMX_ASSERT(requests_[senderIndex] == MPI_REQUEST_NULL,
               "The previous coordinate receive must be consumed before posting a new one");
    const int ppRank = ppCommManagers_[senderIndex].ppRank.rankId;
    if (GMX_THREAD_MPI)
    {
        GMX_ASSERT(GpuConfigurationCapabilities::PpPmeDirectComm,
                   "Direct PME-PP communication not supported with this backend and thread-MPI; "
                   "use lib-MPI instead.");
        // The data will be pushed directly from the PP task, so we
        // only receive the event that signals when it has arrived.
        MPI_Irecv(&ppCommManagers_[senderIndex].sync,
                  sizeof(GpuEventSynchronizer*), // NOLINT(bugprone-sizeof-expression)
                  MPI_BYTE,
                  ppRank,
                  eCommType_COORD_GPU_SYNCHRONIZER,
                  comm_,
                  &(requests_[senderIndex]));
    }
    else
    {
        MPI_Irecv(asRawDevicePointer(recvbuf) + atomOffset,
                  numBytes,
                  MPI_BYTE,
                  ppRank,
                  eCommType_COORD_GPU,
                  comm_,
                  &(requests_[senderIndex]));
    }
#else
    GMX_UNUSED_VALUE(recvbuf);
    GMX_UNUSED_VALUE(atomOffset);
    GMX_UNUSED_VALUE(numBytes);
    GMX_UNUSED_VALUE(senderIndex);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
std::tuple<int, GpuEventSynchronizer*> PmeCoordinateReceiverGpu::Impl::waitForCoordinatesFromAnyPpRank()
{
#if GMX_MPI
    // MPI_Waitany returns in senderIndex the index of one of the
    // requests, i.e. the index of the sender within the set of PP
    // ranks that collaborate with this PME rank.
    int senderIndex = -1;
    MPI_Waitany(requests_.size(), requests_.data(), &senderIndex, MPI_STATUS_IGNORE);
    // Note that coordinates are always transferred, even from
    // empty domains. Thus senderIndex must be non-negative after
    // MPI_Waitany returns.
    GMX_RELEASE_ASSERT(senderIndex >= 0, "Sender index must be valid");

    if (GMX_THREAD_MPI)
    {
        // With thread-MPI the returned event must be enqueued by the
        // caller into the stream that will consume the coordinates.
        GpuEventSynchronizer* event = ppCommManagers_[senderIndex].sync;
        GMX_RELEASE_ASSERT(event != nullptr,
                           "Thread-MPI coordinate receives must provide a GPU transfer event");
        return std::make_tuple(senderIndex, event);
    }
    else
    {
        // With lib-MPI the transfer is already complete and there is no event
        return std::make_tuple(senderIndex, nullptr);
    }
#else
    return std::make_tuple(-1, nullptr);
#endif
}

DeviceStream* PmeCoordinateReceiverGpu::Impl::ppCommStream(int senderIndex)
{
    return ppCommManagers_[senderIndex].stream.get();
}

std::tuple<int, int> PmeCoordinateReceiverGpu::Impl::ppCommAtomRange(int senderIndex)
{
    return ppCommManagers_[senderIndex].atomRange;
}

int PmeCoordinateReceiverGpu::Impl::ppCommNumRanks()
{
    return ppCommManagers_.size();
}

void PmeCoordinateReceiverGpu::Impl::insertAsDependencyIntoStream(int senderIndex, const DeviceStream& stream)
{
    GMX_ASSERT(senderIndex >= 0, "Must have valid sender index");
    ppCommManagers_[senderIndex].ready->markEvent(*ppCommManagers_[senderIndex].stream);
    ppCommManagers_[senderIndex].ready->enqueueWaitEvent(stream);
}

PmeCoordinateReceiverGpu::PmeCoordinateReceiverGpu(MPI_Comm               comm,
                                                   const DeviceContext&   deviceContext,
                                                   gmx::ArrayRef<PpRanks> ppRanks) :
    impl_(new Impl(comm, deviceContext, ppRanks))
{
}

PmeCoordinateReceiverGpu::~PmeCoordinateReceiverGpu() = default;

void PmeCoordinateReceiverGpu::reinitCoordinateReceiver(DeviceBuffer<RVec> d_x)
{
    impl_->reinitCoordinateReceiver(d_x);
}

void PmeCoordinateReceiverGpu::launchReceiveCoordinatesFromPp(DeviceBuffer<RVec> recvbuf,
                                                              int                atomOffset,
                                                              int                numBytes,
                                                              int                senderIndex)
{
    impl_->launchReceiveCoordinatesFromPp(recvbuf, atomOffset, numBytes, senderIndex);
}

std::tuple<int, GpuEventSynchronizer*> PmeCoordinateReceiverGpu::waitForCoordinatesFromAnyPpRank()
{
    return impl_->waitForCoordinatesFromAnyPpRank();
}

DeviceStream* PmeCoordinateReceiverGpu::ppCommStream(int senderIndex)
{
    return impl_->ppCommStream(senderIndex);
}

std::tuple<int, int> PmeCoordinateReceiverGpu::ppCommAtomRange(int senderIndex)
{
    return impl_->ppCommAtomRange(senderIndex);
}

int PmeCoordinateReceiverGpu::ppCommNumRanks()
{
    return impl_->ppCommNumRanks();
}

void PmeCoordinateReceiverGpu::insertAsDependencyIntoStream(int senderIndex, const DeviceStream& stream)
{
    impl_->insertAsDependencyIntoStream(senderIndex, stream);
}

} // namespace gmx
