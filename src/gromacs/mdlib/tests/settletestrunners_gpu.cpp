/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Defines the runner for CUDA version of SETTLE.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "config.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/basedefinitions.h"

#include "settletestrunners.h"

#if GPU_SETTLE_SUPPORTED
#    include "gromacs/gpu_utils/devicebuffer.h"
#endif
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/mdlib/settle_gpu.h"

namespace gmx
{
namespace test
{
class SettleTestData;

#if GPU_SETTLE_SUPPORTED

void SettleDeviceTestRunner::applySettle(SettleTestData* testData,
                                         const t_pbc     pbc,
                                         const bool      updateVelocities,
                                         const bool      calcVirial,
                                         const std::string& /* testDescription */)
{
    const DeviceContext& deviceContext = testDevice_.deviceContext();
    const DeviceStream&  deviceStream  = testDevice_.deviceStream();
    deviceContext.activate();

    auto settleGpu = std::make_unique<SettleGpu>(testData->mtop_, deviceContext, deviceStream);

    settleGpu->set(*testData->idef_);
    PbcAiuc pbcAiuc;
    setPbcAiuc(pbc.ndim_ePBC, pbc.box, &pbcAiuc);

    int numAtoms = testData->numAtoms_;

    DeviceBuffer<Float3> d_x, d_xp, d_v;

    Float3* h_x  = gmx::asGenericFloat3Pointer(testData->x_);
    Float3* h_xp = gmx::asGenericFloat3Pointer(testData->xPrime_);
    Float3* h_v  = gmx::asGenericFloat3Pointer(testData->v_);

    allocateDeviceBuffer(&d_x, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_xp, numAtoms, deviceContext);
    allocateDeviceBuffer(&d_v, numAtoms, deviceContext);

    copyToDeviceBuffer(&d_x, h_x, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyToDeviceBuffer(&d_xp, h_xp, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    if (updateVelocities)
    {
        copyToDeviceBuffer(&d_v, h_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    }
    settleGpu->apply(
            d_x, d_xp, updateVelocities, d_v, testData->reciprocalTimeStep_, calcVirial, testData->virial_, pbcAiuc);

    copyFromDeviceBuffer(h_xp, &d_xp, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    if (updateVelocities)
    {
        copyFromDeviceBuffer(h_v, &d_v, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    }

    freeDeviceBuffer(&d_x);
    freeDeviceBuffer(&d_xp);
    freeDeviceBuffer(&d_v);
}

#else // GPU_SETTLE_SUPPORTED

void SettleDeviceTestRunner::applySettle(SettleTestData* /* testData */,
                                         const t_pbc /* pbc */,
                                         const bool /* updateVelocities */,
                                         const bool /* calcVirial */,
                                         const std::string& /* testDescription */)
{
    GMX_UNUSED_VALUE(testDevice_);
    FAIL() << "Dummy SETTLE GPU function was called instead of the real one in the SETTLE test.";
}

#endif // GPU_SETTLE_SUPPORTED

} // namespace test
} // namespace gmx
