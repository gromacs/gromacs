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
 * \brief The CPU stub for the state propagator data class.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/mdtypes/state_propagator_data_gpu.h"

#if GMX_GPU == GMX_GPU_NONE
namespace gmx
{

class StatePropagatorDataGpu::Impl
{
};

StatePropagatorDataGpu::StatePropagatorDataGpu(const void* /* pmeStream       */,
                                               const void* /* localStream     */,
                                               const void* /* nonLocalStream  */,
                                               const void* /* deviceContext   */,
                                               GpuApiCallBehavior /* transferKind    */,
                                               int /* paddingSize     */,
                                               gmx_wallcycle* /*   wcycle */) :
    impl_(nullptr)
{
}

StatePropagatorDataGpu::StatePropagatorDataGpu(const void* /* pmeStream       */,
                                               const void* /* deviceContext   */,
                                               GpuApiCallBehavior /* transferKind    */,
                                               int /* paddingSize     */,
                                               gmx_wallcycle* /*   wcycle */) :
    impl_(nullptr)
{
}

StatePropagatorDataGpu::StatePropagatorDataGpu(StatePropagatorDataGpu&& /* other */) noexcept = default;

StatePropagatorDataGpu& StatePropagatorDataGpu::operator=(StatePropagatorDataGpu&& /* other */) noexcept = default;

StatePropagatorDataGpu::~StatePropagatorDataGpu() = default;

void StatePropagatorDataGpu::reinit(int /* numAtomsLocal */, int /* numAtomsAll   */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
}

std::tuple<int, int> StatePropagatorDataGpu::getAtomRangesFromAtomLocality(AtomLocality /* atomLocality */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
    return std::make_tuple(0, 0);
}

DeviceBuffer<float> StatePropagatorDataGpu::getCoordinates()
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
    return DeviceBuffer<float>{};
}

GpuEventSynchronizer* StatePropagatorDataGpu::getCoordinatesReadyOnDeviceEvent(
        AtomLocality /* atomLocality */,
        const SimulationWorkload& /* simulationWork */,
        const StepWorkload& /* stepWork       */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
    return nullptr;
}

void StatePropagatorDataGpu::waitCoordinatesCopiedToDevice(AtomLocality /* atomLocality */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
}

GpuEventSynchronizer* StatePropagatorDataGpu::xUpdatedOnDevice()
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
    return nullptr;
}

void StatePropagatorDataGpu::copyCoordinatesToGpu(const gmx::ArrayRef<const gmx::RVec> /* h_x */,
                                                  AtomLocality /* atomLocality */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
}

void StatePropagatorDataGpu::waitCoordinatesReadyOnHost(AtomLocality /* atomLocality */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
}

void StatePropagatorDataGpu::copyCoordinatesFromGpu(gmx::ArrayRef<gmx::RVec> /* h_x          */,
                                                    AtomLocality /* atomLocality */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
}


DeviceBuffer<float> StatePropagatorDataGpu::getVelocities()
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
    return DeviceBuffer<float>{};
}

void StatePropagatorDataGpu::copyVelocitiesToGpu(const gmx::ArrayRef<const gmx::RVec> /* h_v */,
                                                 AtomLocality /* atomLocality */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
}

GpuEventSynchronizer* StatePropagatorDataGpu::getVelocitiesReadyOnDeviceEvent(AtomLocality /* atomLocality */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
    return nullptr;
}

void StatePropagatorDataGpu::copyVelocitiesFromGpu(gmx::ArrayRef<gmx::RVec> /* h_v          */,
                                                   AtomLocality /* atomLocality */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
}

void StatePropagatorDataGpu::waitVelocitiesReadyOnHost(AtomLocality /* atomLocality */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
}


DeviceBuffer<float> StatePropagatorDataGpu::getForces()
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
    return DeviceBuffer<float>{};
}

void StatePropagatorDataGpu::copyForcesToGpu(const gmx::ArrayRef<const gmx::RVec> /* h_f          */,
                                             AtomLocality /* atomLocality */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
}

GpuEventSynchronizer* StatePropagatorDataGpu::getForcesReadyOnDeviceEvent(AtomLocality /* atomLocality */,
                                                                          bool /* useGpuFBufferOps */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
    return nullptr;
}

GpuEventSynchronizer* StatePropagatorDataGpu::fReducedOnDevice()
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
    return nullptr;
}

void StatePropagatorDataGpu::copyForcesFromGpu(gmx::ArrayRef<gmx::RVec> /* h_f          */,
                                               AtomLocality /* atomLocality */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
}

void StatePropagatorDataGpu::waitForcesReadyOnHost(AtomLocality /* atomLocality */)
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
}


void* StatePropagatorDataGpu::getUpdateStream()
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
    return nullptr;
}

int StatePropagatorDataGpu::numAtomsLocal()
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
    return 0;
}

int StatePropagatorDataGpu::numAtomsAll()
{
    GMX_ASSERT(false,
               "A CPU stub method from GPU state propagator data was called instead of one from "
               "GPU implementation.");
    return 0;
}

} // namespace gmx

#endif // GMX_GPU == GMX_GPU_NONE
