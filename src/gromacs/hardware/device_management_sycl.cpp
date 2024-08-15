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

#include "config.h"

#include <algorithm>
#include <map>
#include <optional>
#include <tuple>
#include <vector>

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/device_management_sycl_intel_device_ids.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/mpiinfo.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "device_information.h"


static std::optional<std::tuple<int, int>> parseHardwareVersionNvidia(const std::string& archName)
{
    // archName could be either '8.6' (DPC++) or 'sm_86' (AdaptiveCpp/hipSYCL)
    try
    {
        if (gmx::startsWith(archName, "sm_"))
        {
            int numericCCBegin = 3; // first character after sm_
            int numericCCEnd   = archName.length();
            if (numericCCEnd - numericCCBegin >= 2)
            {
                std::string majorStr = archName.substr(numericCCBegin, numericCCEnd - numericCCBegin - 1);
                std::string minorStr = archName.substr(numericCCEnd - 1, 1);
                const int   major    = gmx::fromStdString<int>(majorStr);
                const int   minor    = gmx::fromStdString<int>(minorStr);
                return std::make_tuple(major, minor);
            }
        }
        else if (gmx::contains(archName, "."))
        {
            const std::vector<std::string> ccTokens = gmx::splitDelimitedString(archName, '.');
            if (ccTokens.size() == 2)
            {
                const int major = gmx::fromStdString<int>(ccTokens[0]);
                const int minor = gmx::fromStdString<int>(ccTokens[1]);
                return std::make_tuple(major, minor);
            }
        }
    }
    catch (gmx::InvalidInputError)
    {
        return std::nullopt;
    }
    return std::nullopt;
}


static std::optional<std::tuple<int, int>> getHardwareVersionNvidia(const sycl::device& device)
{
    /* First, check device::info::version:
     * - AdaptiveCpp/hipSYCL supports that since AdaptiveCpp 2023.10.0 (merged in July 2023),
     * - Intel DPC++ supports that since 2023.2.0 (merged in July 2023).
     * If device::info::version cannot be parsed, fall back on backend-specific solutions.
     * Fallbacks can be removed once we no longer support older versions. */
    const std::string deviceVersion = device.get_info<sycl::info::device::version>();
    if (auto result = parseHardwareVersionNvidia(deviceVersion); result.has_value())
    {
        return result;
    }
#if (GMX_SYCL_HIPSYCL && GMX_HIPSYCL_HAVE_CUDA_TARGET) // hipSYCL uses CUDA Runtime API
    const int             nativeDeviceId = sycl::get_native<sycl::backend::cuda>(device);
    struct cudaDeviceProp prop;
    cudaError_t           status = cudaGetDeviceProperties(&prop, nativeDeviceId);
    if (status == cudaSuccess)
    {
        return std::make_tuple(prop.major, prop.minor);
    }
    else
    {
        return std::nullopt;
    }
#elif (GMX_SYCL_DPCPP && defined(SYCL_EXT_ONEAPI_BACKEND_CUDA))
    // oneAPI uses CUDA Driver API, but does not link the application to it
    // Instead, we have to use info::device::backend_version, and parse it
    const std::string ccStr = device.get_info<sycl::info::device::backend_version>();
    return parseHardwareVersionNvidia(ccStr);
#else
    GMX_UNUSED_VALUE(device);
    return std::nullopt;
#endif
}

static std::optional<std::tuple<int, int, int>> parseHardwareVersionAmd(const std::string& archName)
{
    // archName is something like 'gfx90a:sramecc+:xnack-' or 'gfx1034'
    std::vector<std::string> archNameTokens = gmx::splitDelimitedString(archName, ':');
    if (!gmx::startsWith(archNameTokens[0], "gfx"))
    {
        return std::nullopt;
    }
    const std::string archNumber = archNameTokens[0].substr(3);
    // Now we have the main part, e.g., "906", "90a", or "1033"
    // Last two chars are minor version and patch, the first one or two chars are major
    const int split = archNumber.size() - 2;
    if (split != 1 && split != 2)
    {
        return std::nullopt;
    }
    try
    {
        const int  major     = gmx::fromStdString<int>(archNumber.substr(0, split));
        const int  minor     = gmx::fromStdString<int>(archNumber.substr(split, 1));
        const char patchChar = archNumber[split + 1];
        int        patch;
        if (patchChar >= '0' && patchChar <= '9')
        {
            patch = patchChar - '0';
        }
        else if (patchChar >= 'a' && patchChar <= 'z')
        {
            patch = 10 + patchChar - 'a';
        }
        // If we have not thrown any exceptions, save the data
        return std::make_tuple(major, minor, patch);
    }
    catch (gmx::InvalidInputError)
    {
        return std::nullopt;
    }
}


static std::optional<std::tuple<int, int, int>> getHardwareVersionAmd(const sycl::device& device)
{
    /* First, check device::info::version:
     * - AdaptiveCpp/hipSYCL supports that since AdaptiveCpp 2023.10.0 (merged in July 2023),
     * - Intel DPC++ supports that since 2023.2.0 (merged in July 2023).
     * If device::info::version cannot be parsed, fall back to backend-specific solutions.
     * Fallbacks can be removed once we no longer support older versions. */
    const std::string deviceVersion = device.get_info<sycl::info::device::version>();
    if (auto result = parseHardwareVersionAmd(deviceVersion); result.has_value())
    {
        return result;
    }
#if (GMX_SYCL_HIPSYCL && GMX_HIPSYCL_HAVE_HIP_TARGET)
    // Fall back on the native device query
    const int              nativeDeviceId = sycl::get_native<sycl::backend::hip>(device);
    struct hipDeviceProp_t prop;
    hipError_t             status = hipGetDeviceProperties(&prop, nativeDeviceId);
    if (status != hipSuccess)
    {
        return std::nullopt;
    }
    // prop.major and prop.minor indicate the closest CUDA CC
    // gcnArch is deprecated, so we have to parse gcnArchName
    return parseHardwareVersionAmd(prop.gcnArchName);
#elif (GMX_SYCL_DPCPP)
    // Device name might contain the desired string, but it depends on the ROCm version
    const std::string deviceName = device.get_info<sycl::info::device::version>();
    return parseHardwareVersionAmd(deviceName);
#else
    GMX_UNUSED_VALUE(device);
    return std::nullopt;
#endif
}

static std::optional<std::tuple<int, int, int>> getHardwareVersionIntel(const sycl::device& device)
{
    // ext::intel::info::device::device_id is supported since oneAPI 2023.1
#if defined(SYCL_EXT_INTEL_DEVICE_INFO) && SYCL_EXT_INTEL_DEVICE_INFO >= 5
    if (device.has(sycl::aspect::ext_intel_device_id))
    {
        auto pciId = device.get_info<sycl::ext::intel::info::device::device_id>();
        return getIntelHardwareVersionFromPciExpressID(pciId);
    }
#endif
    // Otherwise, try to find PCI ID in the device string. That only works for old drivers
    const std::string deviceName = device.get_info<sycl::info::device::name>();
    // Find " [0xABCD]" and extract 0xABCD:
    const size_t prefixStart = deviceName.find(" [0x");
    const size_t pciIdStart  = prefixStart + 2;
    const size_t pciIdEnd    = pciIdStart + 6;
    const size_t suffix      = pciIdEnd; // Should be ']'
    if (prefixStart == std::string::npos || deviceName.length() < suffix + 1 || deviceName[suffix] != ']')
    {
        return std::nullopt;
    }
    const auto    pciIdStr = deviceName.substr(pciIdStart, 6); // 0xABCD
    unsigned long pciId;
    try
    {
        pciId = std::stoul(pciIdStr, nullptr, 16);
    }
    catch (std::invalid_argument)
    {
        return std::nullopt;
    }

    // Now, find the matching device
    return getIntelHardwareVersionFromPciExpressID(pciId);
}

static gmx::GpuAwareMpiStatus getDeviceGpuAwareMpiStatus(const sycl::backend backend)
{
    if (GMX_LIB_MPI == 0)
    {
        return gmx::GpuAwareMpiStatus::NotSupported;
    }
    switch (backend)
    {
#if GMX_SYCL_DPCPP
        case sycl::backend::opencl: return gmx::GpuAwareMpiStatus::NotSupported;
        case sycl::backend::ext_oneapi_level_zero: return gmx::checkMpiZEAwareSupport();
        case sycl::backend::ext_oneapi_cuda: return gmx::checkMpiCudaAwareSupport();
        case sycl::backend::ext_oneapi_hip: return gmx::checkMpiHipAwareSupport();
#elif GMX_SYCL_ACPP
        case sycl::backend::cuda: return gmx::checkMpiCudaAwareSupport();
        case sycl::backend::hip: return gmx::checkMpiHipAwareSupport();
        case sycl::backend::level_zero: return gmx::checkMpiZEAwareSupport();
#endif
        default: return gmx::GpuAwareMpiStatus::NotSupported;
    }
}

void warnWhenDeviceNotTargeted(const gmx::MDLogger& /* mdlog */, const DeviceInformation& /* deviceInfo */)
{
}

bool isDeviceDetectionFunctional(std::string* errorMessage)
{
    try
    {
        const std::vector<sycl::platform> platforms = sycl::platform::get_platforms();
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
 * \param[in]  syclDevice              SYCL device handle.
 * \param[in]  deviceVendor            Device vendor.
 * \param[in]  supportedSubGroupSizes  List of supported sub-group sizes as reported by the device.
 * \returns                            The status enumeration value for the checked device.
 */
static DeviceStatus isDeviceCompatible(const sycl::device&           syclDevice,
                                       const DeviceVendor gmx_unused deviceVendor,
                                       gmx::ArrayRef<const int>      supportedSubGroupSizes)
{
    try
    {
        if (getenv("GMX_GPU_DISABLE_COMPATIBILITY_CHECK") != nullptr)
        {
            // Assume the device is compatible because checking has been disabled.
            return DeviceStatus::Compatible;
        }

        if (syclDevice.get_info<sycl::info::device::local_mem_type>() == sycl::info::local_mem_type::none)
        {
            // While some kernels (leapfrog) can run without shared/local memory, this is a bad sign
            return DeviceStatus::Incompatible;
        }

// Ensure any changes are in sync with nbnxm_sycl_kernel.h
#if GMX_GPU_NB_CLUSTER_SIZE == 4
        const std::vector<int> compiledNbnxmSubGroupSizes{ 8 };
#elif GMX_GPU_NB_CLUSTER_SIZE == 8
#    if GMX_SYCL_HIPSYCL && !(GMX_HIPSYCL_HAVE_HIP_TARGET)
        const std::vector<int> compiledNbnxmSubGroupSizes{ 32 }; // Only NVIDIA
#    elif GMX_SYCL_HIPSYCL && (GMX_HIPSYCL_HAVE_HIP_TARGET && !GMX_HIPSYCL_ENABLE_AMD_RDNA_SUPPORT)
        const std::vector<int> compiledNbnxmSubGroupSizes{ 64 }; // Only AMD GCN and CDNA
#    else
        const std::vector<int> compiledNbnxmSubGroupSizes{ 32, 64 };
#    endif
#endif

        const auto subGroupSizeSupportedByDevice = [&supportedSubGroupSizes](const int sgSize) -> bool {
            return std::find(supportedSubGroupSizes.begin(), supportedSubGroupSizes.end(), sgSize)
                   != supportedSubGroupSizes.end();
        };
        if (std::none_of(compiledNbnxmSubGroupSizes.begin(),
                         compiledNbnxmSubGroupSizes.end(),
                         subGroupSizeSupportedByDevice))
        {
#if GMX_SYCL_HIPSYCL && GMX_HIPSYCL_HAVE_HIP_TARGET && !GMX_HIPSYCL_ENABLE_AMD_RDNA_SUPPORT
            if (supportedSubGroupSizes.size() == 1 && supportedSubGroupSizes[0] == 32
                && deviceVendor == DeviceVendor::Amd)
            {
                return DeviceStatus::IncompatibleAmdRdnaNotTargeted;
            }
#endif
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
    catch (sycl::exception const&) // in case a driver bug causes get_info to throw
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
static bool isDeviceFunctional(const sycl::device& syclDevice, std::string* errorMessage)
{
    static const int numThreads = 8;
    try
    {
        sycl::queue queue(syclDevice, { sycl::property::queue::in_order() });
        std::unique_ptr<int, std::function<void(int*)>> buffer = {
            sycl::malloc_device<int>(numThreads, queue), [=](int* ptr) { sycl::free(ptr, queue); }
        };
        int* d_buffer = buffer.get();
        queue.submit([&](sycl::handler& cgh) {
                 sycl::range<1> range{ numThreads };
                 cgh.parallel_for<DummyKernel>(
                         range, [=](sycl::id<1> threadId) { d_buffer[threadId] = threadId.get(0); });
             }).wait_and_throw();

        std::vector<int> h_buffer(numThreads);
        queue.copy<int>(d_buffer, h_buffer.data(), numThreads).wait_and_throw();
        for (int i = 0; i < numThreads; i++)
        {
            if (h_buffer[i] != i)
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
            errorMessage->assign(gmx::formatString("Unable to run dummy kernel on device %s: %s",
                                                   syclDevice.get_info<sycl::info::device::name>().c_str(),
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

    DeviceStatus supportStatus = isDeviceCompatible(
            deviceInfo.syclDevice, deviceInfo.deviceVendor, deviceInfo.supportedSubGroupSizes);
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

/* In DPC++, the same physical device can appear as different virtual devices provided
 * by different backends (e.g., the same GPU can be accessible via both OpenCL and L0).
 * Thus, using devices from two backends is more likely to be a user error than the
 * desired behavior. In this function, we choose L0 if the MPI library is Intel MPI and
 * the user has opted into using its GPU-aware functionality. Otherwise, we choose the
 * backend with the most compatible devices. In case of a tie, we choose OpenCL (if
 * present), or some arbitrary backend among those with the most devices.
 *
 * In hipSYCL, this problem is unlikely to manifest. It has (as of 2021-03-03) another
 * issues: D2D copy between different backends is not allowed. We don't use D2D in
 * SYCL yet. Additionally, hipSYCL does not implement the `sycl::platform::get_backend()`
 * function.
 * Thus, we only do the backend filtering with DPC++.
 * */
#if GMX_SYCL_DPCPP
static std::optional<sycl::backend> chooseBestBackend(const std::vector<std::unique_ptr<DeviceInformation>>& deviceInfos)
{
    // Count the number of compatible devices per backend
    std::map<sycl::backend, int> countDevicesByBackend; // Default initialized with zeros
    for (const auto& deviceInfo : deviceInfos)
    {
        if (deviceInfo->status == DeviceStatus::Compatible)
        {
            const sycl::backend backend = deviceInfo->syclDevice.get_platform().get_backend();
            ++countDevicesByBackend[backend];
        }
    }
    // Prefer L0 backend if GROMACS might be using GPU-aware Intel MPI.
    if (const gmx::GpuAwareMpiStatus status = gmx::checkMpiZEAwareSupport();
        gmx::usingIntelMpi()
        && ((status == gmx::GpuAwareMpiStatus::Supported) || (status == gmx::GpuAwareMpiStatus::Forced)))
    {
        if (countDevicesByBackend[sycl::backend::ext_oneapi_level_zero] > 0)
        {
            return sycl::backend::ext_oneapi_level_zero;
        }
        // If we get here, then no devices with L0 backend were
        // detected. This can be fine if the user is not running a
        // simulation that will use GPU-aware MPI.
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
        const int devicesInOpenCL = countDevicesByBackend[sycl::backend::opencl];
        if (devicesInOpenCL == backendWithMostDevices->second)
        {
            // Prefer OpenCL backend as more stable, if it has as many devices as others
            return sycl::backend::opencl;
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

static std::vector<sycl::device> partitionDevices(const std::vector<sycl::device>&& devices)
{
#if GMX_SYCL_DPCPP
    std::vector<sycl::device> retVal;
    for (const auto& device : devices)
    {
        using sycl::info::partition_property, sycl::info::partition_affinity_domain;
        // Get the partition affinity domains supported by this device for
        // partitioning into at least two sub-devices along that affinity domain.
        const std::vector<partition_affinity_domain> supportedPartitionAffinityDomains =
                device.get_info<sycl::info::device::partition_affinity_domains>();
        // In principle, the above information is enough to decide
        // whether to attempt NUMA-based device partitioning. But the
        // driver changed in late 2023 to show PVC tiles as devices by
        // default and other infrastructure is still catching up. So
        // we make a further check that some sub-devices exist.
        const unsigned int subDeviceCount =
                device.get_info<sycl::info::device::partition_max_sub_devices>();
        // Can this device be partitioned into two or more sub-devices along NUMA domains?
        if ((subDeviceCount > 0)
            && (std::find(supportedPartitionAffinityDomains.begin(),
                          supportedPartitionAffinityDomains.end(),
                          partition_affinity_domain::numa)
                != supportedPartitionAffinityDomains.end()))
        {
            /* Split the device along NUMA domains into sub-devices.
             * For multi-tile Intel GPUs, this corresponds to
             * individual tiles.  If create_sub_devices() is called
             * for a device that does not support NUMA-based partitioning
             * it would throw, but the above logic avoids that. */
            const std::vector<sycl::device> subDevices =
                    device.create_sub_devices<partition_property::partition_by_affinity_domain>(
                            partition_affinity_domain::numa);
            retVal.insert(retVal.end(), subDevices.begin(), subDevices.end());
        }
        else
        {
            retVal.push_back(device);
        }
    }
    return retVal;
#else
    // For hipSYCL, we don't even bother splitting the devices
    return devices;
#endif
}

std::vector<std::unique_ptr<DeviceInformation>> findDevices()
{
    std::vector<std::unique_ptr<DeviceInformation>> deviceInfos(0);
    const std::vector<sycl::device> allDevices = sycl::device::get_devices(sycl::info::device_type::gpu);
    const std::vector<sycl::device> devices = partitionDevices(std::move(allDevices));
    deviceInfos.reserve(devices.size());
    for (const auto& syclDevice : devices)
    {
        deviceInfos.emplace_back(std::make_unique<DeviceInformation>());

        size_t i = deviceInfos.size() - 1;

        deviceInfos[i]->id         = i;
        deviceInfos[i]->syclDevice = syclDevice;
        deviceInfos[i]->deviceVendor =
                getDeviceVendor(syclDevice.get_info<sycl::info::device::vendor>().c_str());

        deviceInfos[i]->gpuAwareMpiStatus = getDeviceGpuAwareMpiStatus(syclDevice.get_backend());

        deviceInfos[i]->supportedSubGroupSizes.clear();
        try
        {
            const auto sgSizes = syclDevice.get_info<sycl::info::device::sub_group_sizes>();
            GMX_RELEASE_ASSERT(sgSizes.size() <= deviceInfos[i]->supportedSubGroupSizes.capacity(),
                               "Device supports too many subgroup sizes");
            for (int sgSize : sgSizes)
            {
                deviceInfos[i]->supportedSubGroupSizes.push_back(sgSize);
            }
        }
        catch (std::exception)
        {
            deviceInfos[i]->supportedSubGroupSizes.clear();
            // The device will be marked incompatible by checkDevice below
        }

        deviceInfos[i]->status = checkDevice(i, *deviceInfos[i]);

        deviceInfos[i]->hardwareVersionMajor = std::nullopt;
        deviceInfos[i]->hardwareVersionMinor = std::nullopt;
        deviceInfos[i]->hardwareVersionPatch = std::nullopt;
        if (deviceInfos[i]->deviceVendor == DeviceVendor::Nvidia)
        {
            if (const auto computeCapability = getHardwareVersionNvidia(syclDevice);
                computeCapability.has_value())
            {
                deviceInfos[i]->hardwareVersionMajor = std::get<0>(*computeCapability);
                deviceInfos[i]->hardwareVersionMinor = std::get<1>(*computeCapability);
                deviceInfos[i]->hardwareVersionPatch = std::nullopt;
            }
        }
        if (deviceInfos[i]->deviceVendor == DeviceVendor::Amd)
        {
            if (const auto gfxVersion = getHardwareVersionAmd(syclDevice); gfxVersion.has_value())
            {
                deviceInfos[i]->hardwareVersionMajor = std::get<0>(*gfxVersion);
                deviceInfos[i]->hardwareVersionMinor = std::get<1>(*gfxVersion);
                deviceInfos[i]->hardwareVersionPatch = std::get<2>(*gfxVersion);
            }
        }
        if (deviceInfos[i]->deviceVendor == DeviceVendor::Intel)
        {
            if (const auto hwVersion = getHardwareVersionIntel(syclDevice); hwVersion.has_value())
            {
                deviceInfos[i]->hardwareVersionMajor = std::get<0>(*hwVersion);
                deviceInfos[i]->hardwareVersionMinor = std::get<1>(*hwVersion);
                deviceInfos[i]->hardwareVersionPatch = std::get<2>(*hwVersion);
            }
        }
    }
#if GMX_SYCL_DPCPP
    // Now, filter by the backend if we did not disable compatibility check
    if (getenv("GMX_GPU_DISABLE_COMPATIBILITY_CHECK") == nullptr)
    {
        std::optional<sycl::backend> preferredBackend = chooseBestBackend(deviceInfos);
        if (preferredBackend.has_value())
        {
            for (auto& deviceInfo : deviceInfos)
            {
                if (deviceInfo->syclDevice.get_platform().get_backend() != *preferredBackend
                    && deviceInfo->status == DeviceStatus::Compatible)
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

void releaseDevice() {}

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
        std::string deviceArchitecture = "";
        if (deviceInfo.hardwareVersionMajor.has_value())
        {
            deviceArchitecture += gmx::formatString(", architecture %d", *deviceInfo.hardwareVersionMajor);
            if (deviceInfo.hardwareVersionMinor.has_value())
            {
                deviceArchitecture += gmx::formatString(".%d", *deviceInfo.hardwareVersionMinor);
                if (deviceInfo.hardwareVersionPatch.has_value())
                {
                    deviceArchitecture += gmx::formatString(".%d", *deviceInfo.hardwareVersionPatch);
                }
            }
        }
        return gmx::formatString(
                "#%d: name: %s%s, vendor: %s, device version: %s, driver version: %s, status: %s",
                deviceInfo.id,
                deviceInfo.syclDevice.get_info<sycl::info::device::name>().c_str(),
                deviceArchitecture.c_str(),
                deviceInfo.syclDevice.get_info<sycl::info::device::vendor>().c_str(),
                deviceInfo.syclDevice.get_info<sycl::info::device::version>().c_str(),
                deviceInfo.syclDevice.get_info<sycl::info::device::driver_version>().c_str(),
                c_deviceStateString[deviceInfo.status]);
    }
}
