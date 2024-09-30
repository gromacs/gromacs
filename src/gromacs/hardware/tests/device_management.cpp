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
 * Tests for DevicesManager
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_hardware
 */
#include "gmxpre.h"

#include "gromacs/hardware/device_management.h"

#include "config.h"

#include <cstdint>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/hardware/device_information.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(DevicesManagerTest, Serialization)
{
    if (canPerformDeviceDetection(nullptr) && c_canSerializeDeviceInformation)
    {
        std::vector<std::unique_ptr<DeviceInformation>> deviceInfoListIn = findDevices();
        gmx::InMemorySerializer                         writer;
        serializeDeviceInformations(deviceInfoListIn, &writer);
        auto buffer = writer.finishAndGetBuffer();

        gmx::InMemoryDeserializer                       reader(buffer, false);
        std::vector<std::unique_ptr<DeviceInformation>> deviceInfoListOut =
                deserializeDeviceInformations(&reader);

        EXPECT_EQ(deviceInfoListOut.size(), deviceInfoListIn.size())
                << "Number of accessible devices changed after serialization/deserialization.";

        for (int deviceId = 0; deviceId < static_cast<int>(deviceInfoListIn.size()); deviceId++)
        {
            EXPECT_FALSE(deviceInfoListIn[deviceId] == nullptr) << gmx::formatString(
                    "Device #%d information is nullptr before serialization.", deviceId);
            EXPECT_FALSE(deviceInfoListOut[deviceId] == nullptr) << gmx::formatString(
                    "Device #%d information is nullptr after serialization.", deviceId);

            const DeviceInformation& deviceInfoIn  = *deviceInfoListIn[deviceId];
            const DeviceInformation& deviceInfoOut = *deviceInfoListOut[deviceId];
            EXPECT_EQ(deviceInfoIn.status, deviceInfoOut.status) << gmx::formatString(
                    "Device status changed after serialization/deserialization for device #%d.", deviceId);

            EXPECT_EQ(deviceInfoIn.id, deviceInfoOut.id) << gmx::formatString(
                    "Device id changed after serialization/deserialization for device #%d.", deviceId);

#if GMX_GPU_OPENCL
            EXPECT_EQ(deviceInfoIn.oclPlatformId, deviceInfoOut.oclPlatformId) << gmx::formatString(
                    "Device OpenCL platform ID changed after serialization/deserialization for "
                    "device "
                    "#%d.",
                    deviceId);

#endif // GMX_GPU_OPENCL
        }
    }
}

std::string uuidToString(const std::array<std::byte, 16>& uuid)
{
    // Write a string in the frequently used 8-4-4-4-12 format,
    // xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx, where every x represents 4 bits
    constexpr char const* digits = "0123456789ABCDEF";
    std::string           result;
    for (std::size_t i = 0; i < uuid.size(); ++i)
    {
        const std::uint8_t uuidByte = static_cast<std::uint8_t>(uuid[i]);

        result.append(1, digits[(uuidByte >> 4) & 0x0F]);
        result.append(1, digits[uuidByte & 0x0F]);

        if (i == 3 || i == 5 || i == 7 || i == 9)
        {
            result.append(1, '-');
        }
    }

    return result;
}

TEST(UuidStringTest, Works)
{
    std::array<std::byte, 16> id = { { std::byte{ 0 } } };
    EXPECT_EQ("00000000-0000-0000-0000-000000000000", uuidToString(id));
    for (size_t i = 0; i < id.size(); ++i)
    {
        id[i] = static_cast<std::byte>(i);
    }
    EXPECT_EQ("00010203-0405-0607-0809-0A0B0C0D0E0F", uuidToString(id));
}

// We can't actually test UUID detection because the value returned
// will be different on each piece of hardware. This test case is a
// service to GROMACS developers to permit them to check manually that
// UUID detection is working on a platform. When the UUID can be
// detected, this test will still fail, but in doing so it will print
// out the UUID, which can then be compared with the output of another
// tool.
//
// Add HARDWARE_DETECTION to the call to gmx_add_unit_test(HardwareUnitTests)
// Then build and run it with
// "hardware-test --gtest_also_run_disabled_tests --gtest_filter=DevicesManagerTest.DISABLED_DetectsUuid"
// and compare the resulting strings. They should match.
TEST(DevicesManagerTest, DISABLED_DetectsUuid)
{
    if (canPerformDeviceDetection(nullptr))
    {
        std::vector<std::unique_ptr<DeviceInformation>> deviceInfoList = findDevices();
        for (int deviceId = 0; deviceId < static_cast<int>(deviceInfoList.size()); deviceId++)
        {
            SCOPED_TRACE(gmx::formatString("Testing device with ID %d", deviceId));
            ASSERT_TRUE(deviceInfoList[deviceId].get()) << "Invalid handle to DeviceInfo";
            const auto uuid = uuidForDevice(*deviceInfoList[deviceId].get());
            ASSERT_TRUE(uuid.has_value());
            SCOPED_TRACE("Device had UUID " + uuidToString(uuid.value()));
            // Force GoogleTest to print out the messages
            ADD_FAILURE();
        }
    }
}

} // namespace
} // namespace test
} // namespace gmx
