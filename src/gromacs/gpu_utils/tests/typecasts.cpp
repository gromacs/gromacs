/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief
 * Tests for CUDA float3 type layout.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "gmxpre.h"

#include "config.h"

#if GMX_GPU_CUDA || GMX_GPU_HIP


#    include "gromacs/gpu_utils/device_context.h"
#    include "gromacs/hardware/device_management.h"
#    include "gromacs/utility/arrayref.h"

#    include "testutils/test_hardware_environment.h"
#    include "testutils/testasserts.h"
#    include "testutils/testmatchers.h"

#    include "typecasts_runner.h"

namespace gmx
{

namespace test
{

//! Test data in RVec format
static const std::vector<RVec> rVecInput = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 } };

TEST(GpuDataTypesCompatibilityTest, RVecAndFloat3Host)
{
    std::vector<RVec> rVecOutput(rVecInput.size());
    convertRVecToFloat3OnHost(rVecOutput, rVecInput);
    EXPECT_THAT(rVecInput, testing::Pointwise(RVecEq(ulpTolerance(0)), rVecOutput));
}

TEST(GpuDataTypesCompatibilityTest, RVecAndFloat3Device)
{
    for (const auto& testDevice : getTestHardwareEnvironment()->getTestDeviceList())
    {
        testDevice->deviceContext().activate();
        std::vector<RVec> rVecOutput(rVecInput.size());
        convertRVecToFloat3OnDevice(rVecOutput, rVecInput, testDevice.get());
        EXPECT_THAT(rVecInput, testing::Pointwise(RVecEq(ulpTolerance(0)), rVecOutput));
    }
}

} // namespace test
} // namespace gmx

#endif // GMX_GPU_CUDA || GMX_GPU_HIP
