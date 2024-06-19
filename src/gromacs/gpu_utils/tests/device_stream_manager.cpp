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
 * \brief Tests GPU stream manager
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_gpu_utils
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/device_stream_manager.h"

#include "config.h"

#include <initializer_list>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/utility/enumerationhelpers.h"

#include "testutils/test_device.h"
#include "testutils/test_hardware_environment.h"

struct DeviceInformation;

namespace gmx
{

namespace test
{

namespace
{

//! GPU device stream names for outputs.
const EnumerationArray<DeviceStreamType, std::string> c_deviceStreamNames = {
    { "non-bonded local", "non-bonded non-local", "PME", "PME-PP transfer", "update" }
};

/*! \brief Non-GPU builds return nullptr instead of streams,
 * so we have to expect that in such build configurations. */
constexpr bool c_canExpectValidStreams = (GMX_GPU != 0);

//! Helper function to implement readable testing
void expectValidStreams(DeviceStreamManager* manager, std::initializer_list<DeviceStreamType> types)
{
    if (c_canExpectValidStreams)
    {
        for (const DeviceStreamType type : types)
        {
            SCOPED_TRACE("Testing " + c_deviceStreamNames[type] + " stream.");
            EXPECT_TRUE(manager->streamIsValid(type));
        }
    }
}
//! Helper function to implement readable testing
void expectInvalidStreams(DeviceStreamManager* manager, std::initializer_list<DeviceStreamType> types)
{
    for (const DeviceStreamType type : types)
    {
        SCOPED_TRACE("Testing " + c_deviceStreamNames[type] + " stream.");
        EXPECT_FALSE(manager->streamIsValid(type));
    }
}

//! Test fixture
class DeviceStreamManagerTest : public ::testing::Test
{
public:
};

TEST_F(DeviceStreamManagerTest, CorrectStreamsAreReturnedOnNonbondedDevice)
{
    // It would be nice to test that the priority is high when it can
    // be, but that requires calling the same API calls we're testing
    // that we've called, so it is not very useful.
    const bool useTiming = false;

    const auto& testDeviceList = getTestHardwareEnvironment()->getTestDeviceList();
    for (const auto& testDevice : testDeviceList)
    {
        const DeviceInformation& deviceInfo = testDevice->deviceInfo();

        {
            SCOPED_TRACE("No DD, no PME rank, no GPU update");
            SimulationWorkload simulationWork;
            simulationWork.useGpuPme                 = false;
            simulationWork.useGpuPmePpCommunication  = false;
            simulationWork.useGpuUpdate              = false;
            simulationWork.havePpDomainDecomposition = false;
            DeviceStreamManager manager(deviceInfo, simulationWork, useTiming);

            expectValidStreams(&manager, { DeviceStreamType::NonBondedLocal });
            expectInvalidStreams(&manager,
                                 { DeviceStreamType::NonBondedNonLocal,
                                   DeviceStreamType::Pme,
                                   DeviceStreamType::PmePpTransfer,
                                   DeviceStreamType::UpdateAndConstraints });
        }

        {
            SCOPED_TRACE("With DD, no PME rank, no GPU update");
            SimulationWorkload simulationWork;
            simulationWork.useGpuPme                 = false;
            simulationWork.useGpuPmePpCommunication  = false;
            simulationWork.useGpuUpdate              = false;
            simulationWork.havePpDomainDecomposition = true;
            DeviceStreamManager manager(deviceInfo, simulationWork, useTiming);

            expectValidStreams(
                    &manager, { DeviceStreamType::NonBondedLocal, DeviceStreamType::NonBondedNonLocal });
            expectInvalidStreams(&manager,
                                 { DeviceStreamType::Pme,
                                   DeviceStreamType::PmePpTransfer,
                                   DeviceStreamType::UpdateAndConstraints });
        }

        {
            SCOPED_TRACE("No DD, with PME rank, no GPU update");
            SimulationWorkload simulationWork;
            simulationWork.useGpuPme                 = true;
            simulationWork.useGpuPmePpCommunication  = true;
            simulationWork.useGpuUpdate              = false;
            simulationWork.havePpDomainDecomposition = false;
            DeviceStreamManager manager(deviceInfo, simulationWork, useTiming);

            expectValidStreams(&manager,
                               { DeviceStreamType::Pme,
                                 DeviceStreamType::NonBondedLocal,
                                 DeviceStreamType::PmePpTransfer,
                                 DeviceStreamType::UpdateAndConstraints });
            expectInvalidStreams(&manager, { DeviceStreamType::NonBondedNonLocal });
        }

        {
            SCOPED_TRACE("With DD, with PME rank, no GPU update");
            SimulationWorkload simulationWork;
            simulationWork.useGpuPme                 = true;
            simulationWork.useGpuPmePpCommunication  = true;
            simulationWork.useGpuUpdate              = false;
            simulationWork.havePpDomainDecomposition = true;
            DeviceStreamManager manager(deviceInfo, simulationWork, useTiming);

            expectValidStreams(&manager,
                               { DeviceStreamType::Pme,
                                 DeviceStreamType::NonBondedLocal,
                                 DeviceStreamType::NonBondedNonLocal,
                                 DeviceStreamType::PmePpTransfer,
                                 DeviceStreamType::UpdateAndConstraints });
        }

        {
            SCOPED_TRACE("No DD, no PME rank, with GPU update");
            SimulationWorkload simulationWork;
            simulationWork.useGpuPme                 = false;
            simulationWork.useGpuPmePpCommunication  = false;
            simulationWork.useGpuUpdate              = true;
            simulationWork.havePpDomainDecomposition = false;
            DeviceStreamManager manager(deviceInfo, simulationWork, useTiming);

            expectValidStreams(
                    &manager, { DeviceStreamType::NonBondedLocal, DeviceStreamType::UpdateAndConstraints });
            expectInvalidStreams(&manager,
                                 { DeviceStreamType::NonBondedNonLocal,
                                   DeviceStreamType::Pme,
                                   DeviceStreamType::PmePpTransfer });
        }

        {
            SCOPED_TRACE("With DD, no PME rank, with GPU update");
            SimulationWorkload simulationWork;
            simulationWork.useGpuPme                 = false;
            simulationWork.useGpuPmePpCommunication  = false;
            simulationWork.useGpuUpdate              = true;
            simulationWork.havePpDomainDecomposition = true;
            DeviceStreamManager manager(deviceInfo, simulationWork, useTiming);

            expectValidStreams(&manager,
                               { DeviceStreamType::NonBondedLocal,
                                 DeviceStreamType::NonBondedNonLocal,
                                 DeviceStreamType::UpdateAndConstraints });
            expectInvalidStreams(&manager, { DeviceStreamType::Pme, DeviceStreamType::PmePpTransfer });
        }

        {
            SCOPED_TRACE("No DD, with PME rank, with GPU update");
            SimulationWorkload simulationWork;
            simulationWork.useGpuPme                 = true;
            simulationWork.useGpuPmePpCommunication  = true;
            simulationWork.useGpuUpdate              = true;
            simulationWork.havePpDomainDecomposition = false;
            DeviceStreamManager manager(deviceInfo, simulationWork, useTiming);

            expectValidStreams(&manager,
                               { DeviceStreamType::Pme,
                                 DeviceStreamType::NonBondedLocal,
                                 DeviceStreamType::PmePpTransfer,
                                 DeviceStreamType::UpdateAndConstraints });
            expectInvalidStreams(&manager, { DeviceStreamType::NonBondedNonLocal });
        }

        {
            SCOPED_TRACE("With DD, with PME rank, with GPU update");
            SimulationWorkload simulationWork;
            simulationWork.useGpuPme                 = true;
            simulationWork.useGpuPmePpCommunication  = true;
            simulationWork.useGpuUpdate              = true;
            simulationWork.havePpDomainDecomposition = true;
            DeviceStreamManager manager(deviceInfo, simulationWork, useTiming);

            expectValidStreams(&manager,
                               { DeviceStreamType::Pme,
                                 DeviceStreamType::NonBondedLocal,
                                 DeviceStreamType::NonBondedNonLocal,
                                 DeviceStreamType::PmePpTransfer,
                                 DeviceStreamType::UpdateAndConstraints });
        }
    }
}

} // namespace
} // namespace test
} // namespace gmx
