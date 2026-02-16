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

#include <algorithm>

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
            void* sendBuf = reinterpret_cast<void*>(asMpiPointer(d_x) + indStart);
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

/*! \brief Receive coordinate synchronizer pointer from the PP ranks. */
void PmeCoordinateReceiverGpu::Impl::receiveCoordinatesSynchronizerFromPpPeerToPeer(int ppRank)
{
    GMX_ASSERT(GMX_THREAD_MPI,
               "receiveCoordinatesSynchronizerFromPpPeerToPeer is expected to be called only for "
               "Thread-MPI");
    GMX_ASSERT(GpuConfigurationCapabilities::PpPmeDirectComm,
               "Direct PME-PP communication not supported with with the backend and threadMPI; use "
               "libMPI "
               "instead.");

    // Data will be pushed directly from PP task

#if GMX_MPI
    // Receive event from PP task
    MPI_Irecv(&ppCommManagers_[ppRank].sync,
              sizeof(GpuEventSynchronizer*), // NOLINT(bugprone-sizeof-expression)
              MPI_BYTE,
              ppRank,
              eCommType_COORD_GPU_SYNCHRONIZER,
              comm_,
              &(requests_[ppRank]));
#else
    GMX_UNUSED_VALUE(ppRank);
#endif
}

/*! \brief Receive coordinate data using GPU-aware MPI */
void PmeCoordinateReceiverGpu::Impl::launchReceiveCoordinatesFromPpGpuAwareMpi(DeviceBuffer<RVec> recvbuf,
                                                                               int numAtoms,
                                                                               int numBytes,
                                                                               int ppRank,
                                                                               int senderIndex)
{
    GMX_ASSERT(
            GMX_LIB_MPI,
            "launchReceiveCoordinatesFromPpGpuAwareMpi is expected to be called only for Lib-MPI");

#if GMX_MPI
    MPI_Irecv(asMpiPointer(recvbuf) + numAtoms,
              numBytes,
              MPI_BYTE,
              ppRank,
              eCommType_COORD_GPU,
              comm_,
              &(requests_[senderIndex]));
#else
    GMX_UNUSED_VALUE(recvbuf);
    GMX_UNUSED_VALUE(numAtoms);
    GMX_UNUSED_VALUE(numBytes);
    GMX_UNUSED_VALUE(ppRank);
    GMX_UNUSED_VALUE(senderIndex);
#endif
}

std::tuple<int, GpuEventSynchronizer*> PmeCoordinateReceiverGpu::Impl::receivePpCoordinateSendEvent(int senderIndex)
{
#if GMX_MPI
    // Loop until a message is received from a PP rank that transferred
    // a non-zero number of atoms.
    do
    {
        // MPI_Waitany is not available in thread-MPI. However, the
        // MPI_Wait here is not associated with data but is host-side
        // scheduling code to receive a CUDA event, and will be executed
        // in advance of the actual data transfer. Therefore we can
        // receive in order of pipeline stage, still allowing the
        // scheduled GPU-direct comms to initiate out-of-order in their
        // respective streams.

        // Loop until we find a request that has not yet been
        // waited upon.
        while (requests_[senderIndex] == MPI_REQUEST_NULL)
        {
            ++senderIndex;
        }
        MPI_Wait(&(requests_[senderIndex]), MPI_STATUS_IGNORE);
        // Ensure that future calls to this method for later pipeline
        // stages of the same step will not wait upon the same sender.
        requests_[senderIndex] = MPI_REQUEST_NULL;
    } while (ppCommManagers_[senderIndex].ppRank.numAtoms == 0);

    // Return a send event from a PP rank that transferred a non-zero
    // number of atoms.
    return std::make_tuple(senderIndex, ppCommManagers_[senderIndex].sync);
#else
    GMX_UNUSED_VALUE(senderIndex);
    return std::make_tuple(-1, nullptr);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
int PmeCoordinateReceiverGpu::Impl::waitForCoordinatesFromAnyPpRank()
{
#if GMX_LIB_MPI
    // Loop until a message is received from a PP rank that is sending
    // a non-zero number of atoms.
    int senderIndex = -1;
    do
    {
        // Wait on data from any one of the PP sender GPUs.
        //
        // MPI_Waitany returns in senderIndex the index of one of the
        // requests, i.e. the index of the sender within the set of PP
        // ranks that collaborate with this PME rank.
        MPI_Waitany(requests_.size(), requests_.data(), &senderIndex, MPI_STATUS_IGNORE);
        // Note that coordinates are always transferred, even from
        // empty domains. Thus senderIndex must be non-negative after
        // MPI_Waitany returns.
        GMX_ASSERT(senderIndex >= 0, "Sender index must be valid");
    } while (ppCommManagers_[senderIndex].ppRank.numAtoms == 0);
    return senderIndex;
#else
    return -1;
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

int PmeCoordinateReceiverGpu::Impl::ppCommNumRanksSendingParticles()
{
    return std::count_if(ppCommManagers_.begin(),
                         ppCommManagers_.end(),
                         [](const PpCommManager& m) { return m.ppRank.numAtoms > 0; });
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

void PmeCoordinateReceiverGpu::receiveCoordinatesSynchronizerFromPpPeerToPeer(int ppRank)
{
    impl_->receiveCoordinatesSynchronizerFromPpPeerToPeer(ppRank);
}

void PmeCoordinateReceiverGpu::launchReceiveCoordinatesFromPpGpuAwareMpi(DeviceBuffer<RVec> recvbuf,
                                                                         int numAtoms,
                                                                         int numBytes,
                                                                         int ppRank,
                                                                         int senderIndex)
{
    impl_->launchReceiveCoordinatesFromPpGpuAwareMpi(recvbuf, numAtoms, numBytes, ppRank, senderIndex);
}

std::tuple<int, GpuEventSynchronizer*> PmeCoordinateReceiverGpu::receivePpCoordinateSendEvent(int senderIndex)
{
    return impl_->receivePpCoordinateSendEvent(senderIndex);
}

int PmeCoordinateReceiverGpu::waitForCoordinatesFromAnyPpRank()
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

int PmeCoordinateReceiverGpu::ppCommNumRanksSendingParticles()
{
    return impl_->ppCommNumRanksSendingParticles();
}

void PmeCoordinateReceiverGpu::insertAsDependencyIntoStream(int senderIndex, const DeviceStream& stream)
{
    impl_->insertAsDependencyIntoStream(senderIndex, stream);
}

} // namespace gmx
