/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 *  \brief Defines the implementations of device management functions that
 *         are common for CPU, CUDA and OpenCL.
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_hardware
 */
#include "gmxpre.h"

#include <algorithm>

#include "gromacs/hardware/device_management.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

#include "device_information.h"

bool canPerformDeviceDetection(std::string* errorMessage)
{
    return isDeviceDetectionEnabled() && isDeviceDetectionFunctional(errorMessage);
}

bool isDeviceDetectionEnabled()
{
    if (c_binarySupportsGpus)
    {
        return getenv("GMX_DISABLE_GPU_DETECTION") == nullptr;
    }
    else
    {
        return false;
    }
}

DeviceVendor getDeviceVendor(const char* vendorName)
{
    if (vendorName)
    {
        if (strstr(vendorName, "NVIDIA"))
        {
            return DeviceVendor::Nvidia;
        }
        else if (strstr(vendorName, "AMD") || strstr(vendorName, "Advanced Micro Devices"))
        {
            return DeviceVendor::Amd;
        }
        else if (strstr(vendorName, "Intel"))
        {
            return DeviceVendor::Intel;
        }
        else if (strstr(vendorName, "Apple"))
        {
            return DeviceVendor::Apple;
        }
    }
    return DeviceVendor::Unknown;
}

int getDeviceComputeUnitFactor(const DeviceInformation& deviceInfo)
{
    switch (deviceInfo.deviceVendor)
    {
        case DeviceVendor::Amd:
            // TODO: Handle RDNA, where it should be 2
            return 1;
        case DeviceVendor::Intel:
            // There are 16 XVEs per XC on Gen9-Gen12 and Xe
            return 16;
        case DeviceVendor::Nvidia:
        case DeviceVendor::Apple: return 1;
        default:
            // Unknown vendor, we don't know any better.
            GMX_ASSERT(false, "Vendor not supported");
            return 1;
    }
}


std::vector<std::reference_wrapper<DeviceInformation>>
getCompatibleDevices(const std::vector<std::unique_ptr<DeviceInformation>>& deviceInfoList)
{
    // Possible minor over-allocation here, but not important for anything
    std::vector<std::reference_wrapper<DeviceInformation>> compatibleDeviceInfoList;
    compatibleDeviceInfoList.reserve(deviceInfoList.size());
    for (const auto& deviceInfo : deviceInfoList)
    {
        if (deviceInfo->status == DeviceStatus::Compatible)
        {
            compatibleDeviceInfoList.emplace_back(*deviceInfo);
        }
    }
    return compatibleDeviceInfoList;
}

std::vector<int> getCompatibleDeviceIds(gmx::ArrayRef<const std::unique_ptr<DeviceInformation>> deviceInfoList)
{
    // Possible minor over-allocation here, but not important for anything
    std::vector<int> compatibleDeviceIds;
    compatibleDeviceIds.reserve(deviceInfoList.size());
    for (const auto& deviceInfo : deviceInfoList)
    {
        if (deviceInfo->status == DeviceStatus::Compatible)
        {
            compatibleDeviceIds.emplace_back(deviceInfo->id);
        }
    }
    return compatibleDeviceIds;
}

bool deviceIdIsCompatible(gmx::ArrayRef<const std::unique_ptr<DeviceInformation>> deviceInfoList,
                          const int                                               deviceId)
{
    auto foundIt = std::find_if(deviceInfoList.begin(),
                                deviceInfoList.end(),
                                [deviceId](auto& deviceInfo) { return deviceInfo->id == deviceId; });
    if (foundIt == deviceInfoList.end())
    {
        GMX_THROW(gmx::RangeError(gmx::formatString(
                "Device ID %d did not correspond to any of the %zu detected device(s)",
                deviceId,
                deviceInfoList.size())));
    }
    return (*foundIt)->status == DeviceStatus::Compatible;
}

gmx::GpuAwareMpiStatus getMinimalSupportedGpuAwareMpiStatus(
        gmx::ArrayRef<const std::unique_ptr<DeviceInformation>> deviceInfoList)
{
    gmx::GpuAwareMpiStatus minVal                 = gmx::GpuAwareMpiStatus::Supported;
    bool                   foundACompatibleDevice = false;
    for (const auto& deviceInfo : deviceInfoList)
    {
        // In all cases, the level of GPU-aware support that is
        // relevant for further GROMACS decision making is the minimal
        // level found over the compatible devices.
        //
        // There can be incompatible devices detected for several
        // reasons, including a mix of devices from a single vendor,
        // of which some are too old.
        //
        // Also, by default dpcpp makes the same Intel GPU visible
        // from two backends. In such cases, GROMACS chooses a
        // preferred backend and leaves only the matching devices with
        // DeviceStatus::Compatible.
        if (deviceInfo->status == DeviceStatus::Compatible)
        {
            minVal                 = std::min(minVal, deviceInfo->gpuAwareMpiStatus);
            foundACompatibleDevice = true;
        }
    }
    // If there were no compatible devices, then there is no
    // support for GPU-aware MPI.
    if (!foundACompatibleDevice)
    {
        minVal = gmx::GpuAwareMpiStatus::NotSupported;
    }
    return minVal;
}

std::string getDeviceCompatibilityDescription(const gmx::ArrayRef<const std::unique_ptr<DeviceInformation>> deviceInfoList,
                                              int deviceId)
{
    return (deviceId >= static_cast<int>(deviceInfoList.size())
                    ? c_deviceStateString[DeviceStatus::Nonexistent]
                    : c_deviceStateString[deviceInfoList[deviceId]->status]);
}

void serializeDeviceInformations(const std::vector<std::unique_ptr<DeviceInformation>>& deviceInfoList,
                                 gmx::ISerializer*                                      serializer)
{
    GMX_RELEASE_ASSERT(c_canSerializeDeviceInformation,
                       "DeviceInformation for OpenCL/SYCL can not be serialized");
    int numDevices = deviceInfoList.size();
    serializer->doInt(&numDevices);
    for (const auto& deviceInfo : deviceInfoList)
    {
        serializer->doOpaque(reinterpret_cast<char*>(deviceInfo.get()), sizeof(DeviceInformation));
    }
}

std::vector<std::unique_ptr<DeviceInformation>> deserializeDeviceInformations(gmx::ISerializer* serializer)
{
    GMX_RELEASE_ASSERT(c_canSerializeDeviceInformation,
                       "DeviceInformation for OpenCL/SYCL can not be deserialized");
    int numDevices = 0;
    serializer->doInt(&numDevices);
    std::vector<std::unique_ptr<DeviceInformation>> deviceInfoList(numDevices);
    for (int i = 0; i < numDevices; i++)
    {
        deviceInfoList[i] = std::make_unique<DeviceInformation>();
        serializer->doOpaque(reinterpret_cast<char*>(deviceInfoList[i].get()), sizeof(DeviceInformation));
    }
    return deviceInfoList;
}
