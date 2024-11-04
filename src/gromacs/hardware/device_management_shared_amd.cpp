/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 *  \brief Defines shared AMD specific methods for device management
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 * \ingroup module_hardware
 */

#include "gmxpre.h"

#include "config.h"

#if GMX_GPU_HIP
#    include "gromacs/gpu_utils/hiputils.h"
#else
#    error AMD device specific header included in unsupported build
#endif

#include "device_management_shared_amd.h"

std::optional<std::array<std::byte, 16>> getAmdDeviceUuid(int deviceId)
{
    auto uuid = std::make_optional<std::array<std::byte, 16>>();

    hipUUID     hipUuid;
    hipDevice_t hipDevice;
    hipError_t  stat = hipDeviceGet(&hipDevice, deviceId);
    gmx::checkDeviceError(stat, "Error getting device information");
    stat = hipDeviceGetUuid(&hipUuid, hipDevice);
    gmx::checkDeviceError(stat, "Error getting device UUID");
    std::memcpy(uuid.value().begin(), hipUuid.bytes, 16);
    // Check the last 6 bytes for sanity, if they are all 0,
    // revert to the fallback handling
    if (std::all_of(uuid.value().begin() + 10,
                    uuid.value().end(),
                    [](const std::byte c) { return c == static_cast<std::byte>(0); }))
    {
        uuid = std::nullopt;
    }

    return uuid;
}

void doubleCheckGpuAwareMpiWillWork(const DeviceInformation& /* deviceInfo */) {}
