/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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

#include <map>
#include <optional>

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

#include "device_information.h"


void warnWhenDeviceNotTargeted(const gmx::MDLogger& /* mdlog */, const DeviceInformation& /* deviceInfo */)
{
}

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

    if (syclDevice.get_info<cl::sycl::info::device::local_mem_type>() == cl::sycl::info::local_mem_type::none)
    {
        // While some kernels (leapfrog) can run without shared/local memory, this is a bad sign
        return DeviceStatus::Incompatible;
    }

    if (syclDevice.is_host())
    {
        // Host device does not support subgroups or even querying for sub_group_sizes
        return DeviceStatus::Incompatible;
    }

    const std::vector<size_t> supportedSubGroupSizes =
            syclDevice.get_info<cl::sycl::info::device::sub_group_sizes>();

    // Ensure any changes stay in sync with subGroupSize in src/gromacs/nbnxm/sycl/nbnxm_sycl_kernel.cpp
    constexpr size_t requiredSubGroupSizeForNbnxm =
#if defined(HIPSYCL_PLATFORM_ROCM)
            GMX_GPU_NB_CLUSTER_SIZE * GMX_GPU_NB_CLUSTER_SIZE;
#else
            GMX_GPU_NB_CLUSTER_SIZE * GMX_GPU_NB_CLUSTER_SIZE / 2;
#endif

    if (std::find(supportedSubGroupSizes.begin(), supportedSubGroupSizes.end(), requiredSubGroupSizeForNbnxm)
        == supportedSubGroupSizes.end())
    {
        return DeviceStatus::IncompatibleClusterSize;
    }

    /* Host device can not be used, because NBNXM requires sub-groups, which are not supported.
     * Accelerators (FPGAs and their emulators) are not supported.
     * So, the only viable options are CPUs and GPUs. */
    const bool forceCpu = (getenv("GMX_SYCL_FORCE_CPU") != nullptr);

    if (forceCpu && syclDevice.is_cpu())
    {
        return DeviceStatus::Compatible;
    }
    else if (!forceCpu && syclDevice.is_gpu())
    {
        return DeviceStatus::Compatible;
    }
    else
    {
        return DeviceStatus::Incompatible;
    }
}

// Declaring the class here to avoid long unreadable name in the profiler report
//! \brief Class name for test kernel
class DummyKernel;

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
    try
    {
        cl::sycl::queue          queue(syclDevice);
        cl::sycl::buffer<int, 1> buffer(numThreads);
        queue.submit([&](cl::sycl::handler& cgh) {
                 auto d_buffer = buffer.get_access<cl::sycl::access::mode::discard_write>(cgh);
                 cl::sycl::range<1> range{ numThreads };
                 cgh.parallel_for<DummyKernel>(range, [=](cl::sycl::id<1> threadId) {
                     d_buffer[threadId] = threadId.get(0);
                 });
             }).wait_and_throw();
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
            errorMessage->assign(
                    gmx::formatString("Unable to run dummy kernel on device %s: %s",
                                      syclDevice.get_info<cl::sycl::info::device::name>().c_str(),
                                      e.what()));
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

/* In DPCPP, the same physical device can appear as different virtual devices provided
 * by different backends (e.g., the same GPU can be accessible via both OpenCL and L0).
 * Thus, using devices from two backends is more likely to be a user error than the
 * desired behavior. In this function, we choose the backend with the most compatible
 * devices. In case of a tie, we choose OpenCL (if present), or some arbitrary backend
 * among those with the most devices.
 *
 * In hipSYCL, this problem is unlikely to manifest. It has (as of 2021-03-03) another
 * issues: D2D copy between different backends is not allowed. We don't use D2D in
 * SYCL yet. Additionally, hipSYCL does not implement the `sycl::platform::get_backend()`
 * function.
 * Thus, we only do the backend filtering with DPCPP.
 * */
#if GMX_SYCL_DPCPP
static std::optional<cl::sycl::backend>
chooseBestBackend(const std::vector<std::unique_ptr<DeviceInformation>>& deviceInfos)
{
    // Count the number of compatible devices per backend
    std::map<cl::sycl::backend, int> countDevicesByBackend; // Default initialized with zeros
    for (const auto& deviceInfo : deviceInfos)
    {
        if (deviceInfo->status == DeviceStatus::Compatible)
        {
            const cl::sycl::backend backend = deviceInfo->syclDevice.get_platform().get_backend();
            ++countDevicesByBackend[backend];
        }
    }
    // If we have devices from more than one backend...
    if (countDevicesByBackend.size() > 1)
    {
        // Find backend with most devices
        const auto backendWithMostDevices = std::max_element(
                countDevicesByBackend.cbegin(),
                countDevicesByBackend.cend(),
                [](const auto& kv1, const auto& kv2) { return kv1.second < kv2.second; });
        // Count devices provided by OpenCL. Will be zero if no OpenCL devices found.
        const int devicesInOpenCL = countDevicesByBackend[cl::sycl::backend::opencl];
        if (devicesInOpenCL == backendWithMostDevices->second)
        {
            // Prefer OpenCL backend as more stable, if it has as many devices as others
            return cl::sycl::backend::opencl;
        }
        else
        {
            // Otherwise, just return max
            return backendWithMostDevices->first;
        }
    }
    else if (countDevicesByBackend.size() == 1)
    {
        return countDevicesByBackend.cbegin()->first;
    }
    else // No devices found
    {
        return std::nullopt;
    }
}
#endif

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
        deviceInfos[i]->deviceVendor =
                getDeviceVendor(syclDevice.get_info<cl::sycl::info::device::vendor>().c_str());
    }
#if GMX_SYCL_DPCPP
    // Now, filter by the backend if we did not disable compatibility check
    if (getenv("GMX_GPU_DISABLE_COMPATIBILITY_CHECK") == nullptr)
    {
        std::optional<cl::sycl::backend> preferredBackend = chooseBestBackend(deviceInfos);
        if (preferredBackend.has_value())
        {
            for (auto& deviceInfo : deviceInfos)
            {
                if (deviceInfo->syclDevice.get_platform().get_backend() != *preferredBackend)
                {
                    deviceInfo->status = DeviceStatus::NotPreferredBackend;
                }
            }
        }
    }
#endif
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
        return gmx::formatString(
                "#%d: %s, status: %s", deviceInfo.id, "N/A", c_deviceStateString[deviceInfo.status]);
    }
    else
    {
        return gmx::formatString(
                "#%d: name: %s, vendor: %s, device version: %s, status: %s",
                deviceInfo.id,
                deviceInfo.syclDevice.get_info<cl::sycl::info::device::name>().c_str(),
                deviceInfo.syclDevice.get_info<cl::sycl::info::device::vendor>().c_str(),
                deviceInfo.syclDevice.get_info<cl::sycl::info::device::version>().c_str(),
                c_deviceStateString[deviceInfo.status]);
    }
}
