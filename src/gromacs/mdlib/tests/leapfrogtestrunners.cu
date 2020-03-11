/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
 * \brief Runner for CUDA version of the integrator
 *
 * Handles GPU data management and actual numerical integration.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "leapfrogtestrunners.h"

#include <assert.h>

#include <cmath>

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "gromacs/gpu_utils/devicebuffer.cuh"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/leapfrog_gpu.cuh"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdtypes/group.h"

namespace gmx
{
namespace test
{

void integrateLeapFrogGpu(LeapFrogTestData* testData, int numSteps)
{
    DeviceInformation   deviceInfo;
    const DeviceContext deviceContext(deviceInfo);
    const DeviceStream  deviceStream(deviceContext, DeviceStreamPriority::Normal, false);

    int numAtoms = testData->numAtoms_;

    float3* h_x  = reinterpret_cast<float3*>(testData->x_.data());
    float3* h_xp = reinterpret_cast<float3*>(testData->xPrime_.data());
    float3* h_v  = reinterpret_cast<float3*>(testData->v_.data());
    float3* h_f  = reinterpret_cast<float3*>(testData->f_.data());

    float3 *d_x, *d_xp, *d_v, *d_f;

    allocateDeviceBuffer(&d_x, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_xp, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_v, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_f, numAtoms, deviceContext);

    copyToDeviceBuffer(&d_x, h_x, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_xp, h_xp, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_v, h_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_f, h_f, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    auto integrator = std::make_unique<LeapFrogGpu>(deviceContext, deviceStream);

    integrator->set(testData->mdAtoms_, testData->numTCoupleGroups_, testData->mdAtoms_.cTC);

    bool doTempCouple = testData->numTCoupleGroups_ > 0;
    for (int step = 0; step < numSteps; step++)
    {
        // This follows the logic of the CPU-based implementation
        bool doPressureCouple = testData->doPressureCouple_
                                && do_per_step(step + testData->inputRecord_.nstpcouple - 1,
                                               testData->inputRecord_.nstpcouple);
        integrator->integrate(d_x, d_xp, d_v, d_f, testData->timestep_, doTempCouple,
                              testData->kineticEnergyData_.tcstat, doPressureCouple,
                              testData->dtPressureCouple_, testData->velocityScalingMatrix_);
    }

    copyFromDeviceBuffer(h_xp, &d_x, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyFromDeviceBuffer(h_v, &d_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    freeDeviceBuffer(&d_x);
    freeDeviceBuffer(&d_xp);
    freeDeviceBuffer(&d_v);
    freeDeviceBuffer(&d_f);
}

} // namespace test
} // namespace gmx
