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
 * Check if GPUs are available when they should be.
 *
 * This is to test that CI can detect and use GPUs, when they are available.
 * Driver and CUDA mismatch can lead to the tests falling back to the CPU
 * code path quietly, leaving the GPU path untested. This is designed to
 * fail when GPUs should be available, indicated by setting the environment
 * variable \c GMX_TEST_REQUIRED_NUMBER_OF_DEVICES to 1 or greater.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Andrey Alekseenko <al42and@gmail.com>
 */
#include "gmxpre.h"

#include "config.h"

#include <cerrno>
#include <cstdlib>

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/hardware/device_management.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/gmxassert.h"

#include "testutils/test_device.h"
#include "testutils/test_hardware_environment.h"

namespace gmx
{

namespace test
{

static unsigned long int getRequiredNumberOfDevices()
{
    const char* required_num_of_devices_str = getenv("GMX_TEST_REQUIRED_NUMBER_OF_DEVICES");
    if (required_num_of_devices_str == nullptr)
    {
        return 0;
    }
    else
    {
        errno = 0;
        char*          strtol_end;
        const long int required_num_of_devices =
                std::strtol(required_num_of_devices_str, &strtol_end, 10);
        GMX_RELEASE_ASSERT(errno == 0, "Invalid value of GMX_TEST_REQUIRED_NUMBER_OF_DEVICES.");
        GMX_RELEASE_ASSERT(*strtol_end == '\0',
                           "Trailing characters in GMX_TEST_REQUIRED_NUMBER_OF_DEVICES.");
        GMX_RELEASE_ASSERT(required_num_of_devices >= 0,
                           "GMX_TEST_REQUIRED_NUMBER_OF_DEVICES can not be negative.");
        return required_num_of_devices;
    }
}

TEST(DevicesAvailable, ShouldBeAbleToRunOnDevice)
{
    const unsigned long int required_num_of_devices = getRequiredNumberOfDevices();

    if (required_num_of_devices != 0)
    {
        ASSERT_TRUE(GMX_GPU) << "GROMACS was compiled without GPU support, yet "
                                "GMX_TEST_REQUIRED_NUMBER_OF_DEVICES is set.";

        std::string platformString(getGpuImplementationString());

        std::string errorMessage = "Can't perform device detection in " + platformString + ":\n";
        std::string detectionErrorMessage;

        // Test if the detection is working
        ASSERT_TRUE(isDeviceDetectionEnabled())
                << errorMessage << "GPU device detection is disabled by "
                << "GMX_DISABLE_GPU_DETECTION environment variable.";
        ASSERT_TRUE(isDeviceDetectionFunctional(&detectionErrorMessage))
                << errorMessage
                << "GPU device checks failed with the following message: " << detectionErrorMessage;

        // This is to test that the GPU detection test environment is working
        const std::vector<std::unique_ptr<TestDevice>>& testDeviceList =
                getTestHardwareEnvironment()->getTestDeviceList();
        ASSERT_FALSE(testDeviceList.empty())
                << errorMessage << "The test environment has an empty list of capable devices.";
        ASSERT_GE(testDeviceList.size(), required_num_of_devices)
                << errorMessage << "Requested " << required_num_of_devices << " device(s), "
                << "but only " << testDeviceList.size() << " detected in the test environment.";
    }
}

} // namespace test
} // namespace gmx
