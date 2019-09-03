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

StatePropagatorDataGpu::StatePropagatorDataGpu(const void *       /* commandStream */,
                                               const void *       /* gpuContext    */,
                                               GpuApiCallBehavior /* transferKind  */,
                                               int                /* paddingSize   */)
    : impl_(nullptr)
{
}

StatePropagatorDataGpu::~StatePropagatorDataGpu()
{
}

void StatePropagatorDataGpu::reinit(int  /* numAtomsLocal */,
                                    int  /* numAtomsAll   */)
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
}

std::tuple<int, int> StatePropagatorDataGpu::getAtomRangesFromAtomLocality(AtomLocality  /* atomLocality */)
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
    return std::make_tuple(0, 0);
}

DeviceBuffer<float> StatePropagatorDataGpu::getCoordinates()
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
    return DeviceBuffer<float> {};
}

void StatePropagatorDataGpu::copyCoordinatesToGpu(const gmx::ArrayRef<const gmx::RVec>  /* h_x          */,
                                                  AtomLocality                          /* atomLocality */)
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
}

void StatePropagatorDataGpu::copyCoordinatesFromGpu(gmx::ArrayRef<gmx::RVec>  /* h_x          */,
                                                    AtomLocality              /* atomLocality */)
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
}


DeviceBuffer<float> StatePropagatorDataGpu::getVelocities()
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
    return DeviceBuffer<float> {};
}

void StatePropagatorDataGpu::copyVelocitiesToGpu(const gmx::ArrayRef<const gmx::RVec>  /* h_v          */,
                                                 AtomLocality                          /* atomLocality */)
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
}

void StatePropagatorDataGpu::copyVelocitiesFromGpu(gmx::ArrayRef<gmx::RVec>  /* h_v          */,
                                                   AtomLocality              /* atomLocality */)
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
}


DeviceBuffer<float> StatePropagatorDataGpu::getForces()
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
    return DeviceBuffer<float> {};
}

void StatePropagatorDataGpu::copyForcesToGpu(const gmx::ArrayRef<const gmx::RVec>  /* h_f          */,
                                             AtomLocality                          /* atomLocality */)
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
}

void StatePropagatorDataGpu::copyForcesFromGpu(gmx::ArrayRef<gmx::RVec>  /* h_f          */,
                                               AtomLocality              /* atomLocality */)
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
}

void StatePropagatorDataGpu::synchronizeStream()
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
}

int StatePropagatorDataGpu::numAtomsLocal()
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
    return 0;
}

int StatePropagatorDataGpu::numAtomsAll()
{
    GMX_ASSERT(false, "A CPU stub method from GPU state propagator data was called insted of one from GPU implementation.");
    return 0;
}

}      // namespace gmx

#endif // GMX_GPU == GMX_GPU_NONE
