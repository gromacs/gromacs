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
 * \brief Runner for GPU version of the integrator
 *
 * Handles GPU data management and actual numerical integration.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "config.h"

#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"

#include "leapfrogtestrunners.h"

#if GPU_LEAPFROG_SUPPORTED
#    include "gromacs/gpu_utils/devicebuffer.h"
#endif
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/mdlib/leapfrog_gpu.h"
#include "gromacs/mdlib/stat.h"

namespace gmx
{
namespace test
{
class LeapFrogTestData;

#if GPU_LEAPFROG_SUPPORTED

void LeapFrogDeviceTestRunner::integrate(LeapFrogTestData* testData, int numSteps)
{
    const DeviceContext& deviceContext = testDevice_.deviceContext();
    const DeviceStream&  deviceStream  = testDevice_.deviceStream();
    deviceContext.activate();

    int numAtoms = testData->numAtoms_;

    Float3* h_x  = gmx::asGenericFloat3Pointer(testData->x_);
    Float3* h_xp = gmx::asGenericFloat3Pointer(testData->xPrime_);
    Float3* h_v  = gmx::asGenericFloat3Pointer(testData->v_);
    Float3* h_f  = gmx::asGenericFloat3Pointer(testData->f_);

    DeviceBuffer<Float3> d_x, d_xp, d_v, d_f;

    allocateDeviceBuffer(&d_x, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_xp, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_v, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_f, numAtoms, deviceContext);

    copyToDeviceBuffer(&d_x, h_x, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_xp, h_xp, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_v, h_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_f, h_f, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    auto integrator =
            std::make_unique<LeapFrogGpu>(deviceContext, deviceStream, testData->numTCoupleGroups_);

    integrator->set(numAtoms, testData->inverseMasses_, testData->mdAtoms_.cTC);

    bool doTempCouple = testData->numTCoupleGroups_ > 0;
    for (int step = 0; step < numSteps; step++)
    {
        // This follows the logic of the CPU-based implementation
        bool doPressureCouple =
                testData->doPressureCouple_
                && do_per_step(step + testData->inputRecord_.pressureCouplingOptions.nstpcouple - 1,
                               testData->inputRecord_.pressureCouplingOptions.nstpcouple);
        integrator->integrate(d_x,
                              d_xp,
                              d_v,
                              d_f,
                              testData->timestep_,
                              doTempCouple,
                              testData->kineticEnergyData_.tcstat,
                              doPressureCouple,
                              testData->dtPressureCouple_,
                              testData->velocityScalingMatrix_);
    }

    copyFromDeviceBuffer(h_xp, &d_x, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyFromDeviceBuffer(h_v, &d_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);

    freeDeviceBuffer(&d_x);
    freeDeviceBuffer(&d_xp);
    freeDeviceBuffer(&d_v);
    freeDeviceBuffer(&d_f);
}

#else // GPU_LEAPFROG_SUPPORTED

void LeapFrogDeviceTestRunner::integrate(LeapFrogTestData* /* testData */, int /* numSteps */)
{
    GMX_UNUSED_VALUE(testDevice_);
    FAIL() << "Dummy Leap-Frog GPU function was called instead of the real one.";
}

#endif // GPU_LEAPFROG_SUPPORTED

} // namespace test
} // namespace gmx
