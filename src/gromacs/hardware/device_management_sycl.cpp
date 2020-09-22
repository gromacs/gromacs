/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 *  \brief Defines the SYCL implementations of the device management.
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *  \author Erik Lindahl <erik.lindahl@gmail.com>
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \ingroup module_hardware
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

#include "device_information.h"


bool isDeviceDetectionFunctional(std::string* errorMessage)
{
    try
    {
        const std::vector<cl::sycl::platform> platforms = cl::sycl::platform::get_platforms();
        // SYCL should always have the "host" platform, but just in case:
        if (platforms.empty() && errorMessage != nullptr)
        {
            errorMessage->assign("No SYCL platforms found.");
        }
        return !platforms.empty();
    }
    catch (const std::exception& e)
    {
        if (errorMessage != nullptr)
        {
            errorMessage->assign(
                    gmx::formatString("Unable to get the list of SYCL platforms: %s", e.what()));
        }
        return false;
    }
}


/*!
 * \brief Checks that device \c deviceInfo is compatible with GROMACS.
 *
 *  For now, only checks that the vendor is Intel and it is a GPU.
 *
 * \param[in]  syclDevice  The SYCL device pointer.
 * \returns                The status enumeration value for the checked device:
 */
static DeviceStatus isDeviceCompatible(const cl::sycl::device& syclDevice)
{
    if (getenv("GMX_GPU_DISABLE_COMPATIBILITY_CHECK") != nullptr)
    {
        // Assume the device is compatible because checking has been disabled.
        return DeviceStatus::Compatible;
    }

    if (syclDevice.is_accelerator()) // FPGAs and FPGA emulators
    {
        return DeviceStatus::Incompatible;
    }
    else
    {
        return DeviceStatus::Compatible;
    }
}


/*!
 * \brief Checks that device \c deviceInfo is sane (ie can run a kernel).
 *
 * Compiles and runs a dummy kernel to determine whether the given
 * SYCL device functions properly.
 *
 *
 * \param[in]  syclDevice      The device info pointer.
 * \param[out] errorMessage    An error message related to a SYCL error.
 * \throws     std::bad_alloc  When out of memory.
 * \returns                    Whether the device passed sanity checks
 */
static bool isDeviceFunctional(const cl::sycl::device& syclDevice, std::string* errorMessage)
{
    static const int numThreads = 8;
    cl::sycl::queue  queue;
    try
    {
        queue = cl::sycl::queue(syclDevice);
        cl::sycl::buffer<int, 1> buffer(numThreads);
        queue.submit([&](cl::sycl::handler& cgh) {
                 auto d_buffer = buffer.get_access<cl::sycl::access::mode::discard_write>(cgh);
                 cgh.parallel_for<class DummyKernel>(numThreads, [=](cl::sycl::id<1> threadId) {
                     d_buffer[threadId] = threadId.get(0);
                 });
             })
                .wait_and_throw();
        const auto h_Buffer = buffer.get_access<cl::sycl::access::mode::read>();
        for (int i = 0; i < numThreads; i++)
        {
            if (h_Buffer[i] != i)
            {
                if (errorMessage != nullptr)
                {
                    errorMessage->assign("Dummy kernel produced invalid values");
                }
                return false;
            }
        }
    }
    catch (const std::exception& e)
    {
        if (errorMessage != nullptr)
        {
            errorMessage->assign(gmx::formatString(
                    "Unable to run dummy kernel on device %s: %s",
                    syclDevice.get_info<cl::sycl::info::device::name>().c_str(), e.what()));
        }
        return false;
    }

    return true;
}

/*!
 * \brief Checks that device \c deviceInfo is compatible and functioning.
 *
 * Checks the given SYCL device for compatibility and runs a dummy kernel on it to determine
 * whether the device functions properly.
 *
 *
 * \param[in] deviceId         Device number (internal to GROMACS).
 * \param[in] deviceInfo       The device info pointer.
 * \returns                    The status of device.
 */
static DeviceStatus checkDevice(size_t deviceId, const DeviceInformation& deviceInfo)
{

    DeviceStatus supportStatus = isDeviceCompatible(deviceInfo.syclDevice);
    if (supportStatus != DeviceStatus::Compatible)
    {
        return supportStatus;
    }

    std::string errorMessage;
    if (!isDeviceFunctional(deviceInfo.syclDevice, &errorMessage))
    {
        gmx_warning("While sanity checking device #%zu, %s", deviceId, errorMessage.c_str());
        return DeviceStatus::NonFunctional;
    }

    return DeviceStatus::Compatible;
}

std::vector<std::unique_ptr<DeviceInformation>> findDevices()
{
    std::vector<std::unique_ptr<DeviceInformation>> deviceInfos(0);
    const std::vector<cl::sycl::device>             devices = cl::sycl::device::get_devices();
    deviceInfos.reserve(devices.size());
    for (const auto& syclDevice : devices)
    {
        deviceInfos.emplace_back(std::make_unique<DeviceInformation>());

        size_t i = deviceInfos.size() - 1;

        deviceInfos[i]->id         = i;
        deviceInfos[i]->syclDevice = syclDevice;
        deviceInfos[i]->status     = checkDevice(i, *deviceInfos[i]);
    }
    return deviceInfos;
}

void setActiveDevice(const DeviceInformation& /*deviceInfo*/) {}

void releaseDevice(DeviceInformation* /* deviceInfo */) {}

std::string getDeviceInformationString(const DeviceInformation& deviceInfo)
{
    bool deviceExists = (deviceInfo.status != DeviceStatus::Nonexistent
                         && deviceInfo.status != DeviceStatus::NonFunctional);

    if (!deviceExists)
    {
        return gmx::formatString("#%d: %s, status: %s", deviceInfo.id, "N/A",
                                 c_deviceStateString[deviceInfo.status]);
    }
    else
    {
        return gmx::formatString(
                "#%d: name: %s, vendor: %s, device version: %s, status: %s", deviceInfo.id,
                deviceInfo.syclDevice.get_info<cl::sycl::info::device::name>().c_str(),
                deviceInfo.syclDevice.get_info<cl::sycl::info::device::vendor>().c_str(),
                deviceInfo.syclDevice.get_info<cl::sycl::info::device::version>().c_str(),
                c_deviceStateString[deviceInfo.status]);
    }
}
