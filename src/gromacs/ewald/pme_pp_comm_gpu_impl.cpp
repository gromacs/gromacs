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
 * \brief May be used to implement PME-PP GPU comm interfaces for non-GPU builds.
 *
 * Currently, reports and exits if any of the interfaces are called.
 * Needed to satisfy compiler on systems, where CUDA is not available.
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "config.h"

#include <cstdint>

#include "gromacs/ewald/pme_pp_comm_gpu.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"

class DeviceContext;
class DeviceStream;
class GpuEventSynchronizer;

#if !GMX_GPU_CUDA && !GMX_GPU_SYCL

namespace gmx
{

/*!\brief \internal Impl class stub. */
class PmePpCommGpu::Impl
{
};

/*!\brief Constructor stub. */
PmePpCommGpu::PmePpCommGpu(MPI_Comm /* comm */,
                           int /* pmeRank */,
                           gmx::HostVector<gmx::RVec>* /* pmeCpuForceBuffer */,
                           const DeviceContext& /* deviceContext */,
                           const DeviceStream& /* deviceStream */,
                           const bool /*useNvshmem*/) :
    impl_(nullptr)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
}

PmePpCommGpu::~PmePpCommGpu() = default;

/*!\brief init PME-PP GPU communication stub */
//NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void PmePpCommGpu::reinit(int /* size */)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication initialization was called instead of the "
               "correct implementation.");
}

//NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void PmePpCommGpu::receiveForceFromPme(RVec* /* recvPtr */, int /* recvSize */, bool /* receivePmeForceToGpu */)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
}

//NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void PmePpCommGpu::sendCoordinatesToPmeFromGpu(DeviceBuffer<RVec> /* sendPtr */,
                                               int /* sendSize */,
                                               GpuEventSynchronizer* /* coordinatesOnDeviceEvent */)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
}

//NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void PmePpCommGpu::sendCoordinatesToPmeFromCpu(RVec* /* sendPtr */, int /* sendSize */)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
}

//NOLINTNEXTLINE(readability-convert-member-functions-to-static)
DeviceBuffer<gmx::RVec> PmePpCommGpu::getGpuForceStagingPtr()
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
    return DeviceBuffer<gmx::RVec>{};
}

//NOLINTNEXTLINE(readability-convert-member-functions-to-static)
GpuEventSynchronizer* PmePpCommGpu::getForcesReadySynchronizer()
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
    return nullptr;
}

DeviceBuffer<uint64_t> PmePpCommGpu::getGpuForcesSyncObj()
{
    GMX_ASSERT(!impl_,
               "A CPU stub for PME-PP GPU communication was called instead of the correct "
               "implementation.");
    return nullptr;
}


} // namespace gmx

#endif // !GMX_GPU_CUDA && !GMX_GPU_SYCL
