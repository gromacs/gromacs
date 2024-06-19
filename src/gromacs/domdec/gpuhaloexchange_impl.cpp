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
 * \brief May be used to implement Domdec CUDA interfaces for non-GPU builds.
 *
 * Currently, reports and exits if any of the interfaces are called.
 * Needed to satisfy compiler on systems, where CUDA is not available.
 *
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_domdec
 */
#include "gmxpre.h"

#include "config.h"

#include <memory>
#include <utility>

#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxassert.h"

class DeviceContext;
class GpuEventSynchronizer;
namespace gmx
{
template<typename T, size_t capacity_>
class FixedCapacityVector;
} // namespace gmx
struct gmx_domdec_t;
struct gmx_wallcycle;

#if !GMX_GPU_CUDA && !GMX_GPU_SYCL

namespace gmx
{

/*!\brief Impl class stub. */
class GpuHaloExchange::Impl
{
};

/*!\brief Constructor stub. */
GpuHaloExchange::GpuHaloExchange(gmx_domdec_t* /* dd */,
                                 int /* dimIndex */,
                                 MPI_Comm /* mpi_comm_mysim */,
                                 const DeviceContext& /* deviceContext */,
                                 int /*pulse */,
                                 gmx_wallcycle* /*wcycle*/) :
    impl_(nullptr)
{
    GMX_ASSERT(false,
               "A CPU stub for GPU Halo Exchange was called insted of the correct implementation.");
}

GpuHaloExchange::~GpuHaloExchange() = default;

GpuHaloExchange::GpuHaloExchange(GpuHaloExchange&&) noexcept = default;

GpuHaloExchange& GpuHaloExchange::operator=(GpuHaloExchange&& other) noexcept
{
    std::swap(impl_, other.impl_);
    return *this;
}

/*!\brief init halo exhange stub. */
void GpuHaloExchange::reinitHalo(DeviceBuffer<RVec> /* d_coordinatesBuffer */,
                                 DeviceBuffer<RVec> /* d_forcesBuffer */)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for GPU Halo Exchange was called insted of the correct implementation.");
}

/*!\brief apply X halo exchange stub. */
GpuEventSynchronizer* GpuHaloExchange::communicateHaloCoordinates(const matrix /* box */,
                                                                  GpuEventSynchronizer* /*dependencyEvent*/)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for GPU Halo Exchange exchange was called insted of the correct "
               "implementation.");
    return nullptr;
}

/*!\brief apply F halo exchange stub. */
void GpuHaloExchange::communicateHaloForces(bool /* accumulateForces */,
                                            FixedCapacityVector<GpuEventSynchronizer*, 2>* /*dependencyEvents*/)
{
    GMX_ASSERT(!impl_,
               "A CPU stub for GPU Halo Exchange was called insted of the correct implementation.");
}

/*!\brief get forces ready on device event stub. */
GpuEventSynchronizer* GpuHaloExchange::getForcesReadyOnDeviceEvent()
{
    GMX_ASSERT(!impl_,
               "A CPU stub for GPU Halo Exchange was called insted of the correct implementation.");
    return nullptr;
}

} // namespace gmx

#endif // !GMX_GPU_CUDA && !GMX_GPU_SYCL
