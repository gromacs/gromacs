/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Tests for NVSHMEM put operation.
 *
 * \author Mahesh Doijade <mdoijade@nvidia.com>
 */
#include "gmxpre.h"

#include "config.h"

#if GMX_NVSHMEM

#    include <nvshmem.h>
#    include <nvshmemx.h>

#    include "gromacs/gpu_utils/gputraits.h"
#    include "gromacs/hardware/device_management.h"
#    include "gromacs/utility/gmxmpi.h"

#    include "testutils/mpitest.h"
#    include "testutils/test_device.h"
#    include "testutils/test_hardware_environment.h"
#    include "testutils/testasserts.h"
#    include "testutils/testmatchers.h"

#    include "nvshmem_simple_put.cuh"

namespace gmx
{

namespace test
{

TEST(NvshmemTests, NvshmemSimplePut)
{
    const int requiredNumOfDevices = 2;
    const int minComputeCapability = 7;
    GMX_MPI_TEST(RequireRankCount<requiredNumOfDevices>);
    int      rank;
    MPI_Comm mpiComm = MPI_COMM_WORLD;
    MPI_Comm_rank(mpiComm, &rank);

    if (!getTestHardwareEnvironment()->hasCompatibleDevices())
    {
        FAIL() << "No compatible GPUs to test on.";
    }

    const auto& testDeviceList = getTestHardwareEnvironment()->getTestDeviceList();

    if (testDeviceList.size() < requiredNumOfDevices)
    {
        FAIL() << "Requested " << requiredNumOfDevices << " device(s), "
               << "but only " << testDeviceList.size() << " detected in the test environment.";
    }

    const DeviceInformation& deviceInfo = testDeviceList[rank].get()->deviceInfo();
    if (deviceInfo.prop.major < minComputeCapability)
    {
        FAIL() << "NVSHMEM requires minimum Volta or higher GPUs";
    }

    nvshmemRunSimplePut(testDeviceList[rank].get());
}

} // namespace test
} // namespace gmx

#endif // GMX_NVSHMEM
