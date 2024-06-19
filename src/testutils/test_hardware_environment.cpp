/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
 * Implements test environment class which performs hardware enumeration for unit tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_testutils
 */

#include "gmxpre.h"

#include "testutils/test_hardware_environment.h"

#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/physicalnodecommunicator.h"

#include "testutils/test_device.h"

struct DeviceInformation;

namespace gmx
{
namespace test
{

//! Mutex for making the test hardware environment
static std::mutex s_testHardwareEnvironmentMutex;
//! The test hardware environment
static std::unique_ptr<TestHardwareEnvironment> s_testHardwareEnvironment;

const TestHardwareEnvironment* getTestHardwareEnvironment()
{
    if (!s_testHardwareEnvironment)
    {
        // Construct and fill the environment
        std::lock_guard<std::mutex> lock(s_testHardwareEnvironmentMutex);
        s_testHardwareEnvironment = std::make_unique<TestHardwareEnvironment>();
    }
    return s_testHardwareEnvironment.get();
}

TestHardwareEnvironment::TestHardwareEnvironment() :
    hardwareInfo_(gmx_detect_hardware(PhysicalNodeCommunicator{ MPI_COMM_WORLD, gmx_physicalnode_id_hash() },
                                      MPI_COMM_WORLD))
{
    // Following the ::testing::Environment protocol
    this->SetUp();

    // Constructing contexts for all compatible GPUs - will be empty on non-GPU builds
    for (const DeviceInformation& compatibleDeviceInfo : getCompatibleDevices(hardwareInfo_->deviceInfoList))
    {
        setActiveDevice(compatibleDeviceInfo);
        std::string description = getDeviceInformationString(compatibleDeviceInfo);
        testDeviceList_.emplace_back(std::make_unique<TestDevice>(description.c_str(), compatibleDeviceInfo));
    }
}

// static
void TestHardwareEnvironment::gmxSetUp()
{
    // Ensure the environment has been constructed
    getTestHardwareEnvironment();
}

// static
void TestHardwareEnvironment::gmxTearDown()
{
    std::lock_guard<std::mutex> lock(s_testHardwareEnvironmentMutex);
    if (!s_testHardwareEnvironment)
    {
        return;
    }
    s_testHardwareEnvironment->testDeviceList_.clear();
    s_testHardwareEnvironment->hardwareInfo_.reset();
}

} // namespace test
} // namespace gmx
