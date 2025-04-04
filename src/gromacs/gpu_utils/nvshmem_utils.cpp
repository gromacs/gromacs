/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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

/*! \libinternal \file
 *
 * \brief Definitions for NVSHMEM halo exchange helper class.
 *
 * \author Mahesh Doijade <mdoijade@nvidia.com>
 *
 * \ingroup module_gpu_utils
 * \inlibraryapi
 */

#include "nvshmem_utils.h"

#include "config.h"

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/utility/basedefinitions.h"
#if GMX_NVSHMEM
#    include <nvshmem.h>
#endif
#if GMX_GPU
#    include "gromacs/gpu_utils/devicebuffer.h"
#endif

gmxNvshmemHandle::gmxNvshmemHandle() : d_ppHaloExSyncBase_(nullptr) {}

gmxNvshmemHandle::~gmxNvshmemHandle()
{
#if GMX_NVSHMEM
    freeDeviceBuffer(&d_ppHaloExSyncBase_);
#endif
}

// NOLINTNEXTLINE readability-convert-member-functions-to-static
void gmxNvshmemHandle::allocateAndInitSignalBufs(int                  totalDimsAndPulses,
                                                 const DeviceContext& deviceContext,
                                                 const DeviceStream*  localStream)
{
#if GMX_GPU
    int totalSyncBufSize = totalDimsAndPulses * numOfPpHaloExSyncBufs;
    reallocateDeviceBuffer(&d_ppHaloExSyncBase_,
                           totalSyncBufSize,
                           &ppHaloExSyncBufSize_,
                           &ppHaloExSyncBufCapacity_,
                           deviceContext,
                           true);
    // If the num of dims/pulses have changed we initialize the signalling
    // buffer to max val.
    if (ppHaloExPerSyncBufSize_ < totalDimsAndPulses)
    {
        // Initialize the signalling buffer with max value.
        ppHaloExPerSyncBufSize_ = totalDimsAndPulses;

        // TODO A class member variable for hostBuffer will tend
        // to reduce the number of host allocations made
        gmx::HostVector<uint64_t> hostBuffer = {
            {}, gmx::HostAllocationPolicy(gmx::PinningPolicy::PinnedIfSupported)
        };
        hostBuffer.resize(totalSyncBufSize, ~0);
        copyToDeviceBuffer<uint64_t>(&d_ppHaloExSyncBase_,
                                     hostBuffer.data(),
                                     0,
                                     static_cast<size_t>(totalSyncBufSize),
                                     *localStream,
                                     GpuApiCallBehavior::Async,
                                     nullptr);
#    if GMX_NVSHMEM
        nvshmemx_sync_all_on_stream(localStream->stream());
#    endif
    }
#else
    GMX_UNUSED_VALUE(totalDimsAndPulses);
    GMX_UNUSED_VALUE(deviceContext);
    GMX_UNUSED_VALUE(localStream);
#endif
}
