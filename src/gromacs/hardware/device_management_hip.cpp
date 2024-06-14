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
 *  \brief Defines the HIP implementations of the device management.
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *  \author Julio Maia <julio.maia@amd.com>
 *
 * \ingroup module_hardware
 */
#include "gmxpre.h"

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/hiputils.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/mpiinfo.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

#include "device_information.h"
#include "device_management.h"

/** Dummy kernel used for sanity checking. */
static __global__ void dummy_kernel() {}

void warnWhenDeviceNotTargeted(const gmx::MDLogger& mdlog, const DeviceInformation& deviceInfo)
{
    if (deviceInfo.status != DeviceStatus::DeviceNotTargeted)
    {
        return;
    }
    gmx::TextLineWrapper wrapper;
    wrapper.settings().setLineLength(80);
    GMX_LOG(mdlog.warning)
            .asParagraph()
            .appendText(wrapper.wrapToString(gmx::formatString(
                    "WARNING: The %s binary does not include support for the HIP architecture of "
                    "the GPU ID #%d (gcn: %s) detected during detection. "
                    "By default, GROMACS supports a number of architectures defined by "
                    "GMX_HIP_TARGET_ARCH, which can be modified during configuration. "
                    "Either your GPU is not covered by this, so you'll need to add the "
                    "architecture to "
                    "the build configuration, or you are trying to use an unsupported consumer "
                    "device with "
                    "rocFFT, "
                    "in which case you should use VkFFT instead as the GPU FFT library. ",
                    gmx::getProgramContext().displayName(),
                    deviceInfo.id,
                    deviceInfo.prop.gcnArchName)));
}

/*! \brief Runs GPU compatibility and sanity checks on the indicated device.
 *
 * Runs a series of checks to determine that the given GPU and underlying HIP
 * driver/runtime functions properly.
 *
 *  As the error handling only permits returning the state of the GPU, this function
 *  does not clear the HIP runtime API status allowing the caller to inspect the error
 *  upon return. Note that this also means it is the caller's responsibility to
 *  reset the HIP runtime state.
 *
 * \todo Currently we do not make a distinction between the type of errors
 *       that can appear during functionality checks. This needs to be improved,
 *       e.g if the dummy test kernel fails to execute with a "device busy message"
 *       we should appropriately report that the device is busy instead of NonFunctional.
 *
 * \todo Introduce errors codes and handle errors more smoothly.
 *
 *
 * \param[in]  deviceInfo  Device information on the device to check.
 * \returns                The status enumeration value for the checked device:
 */
static DeviceStatus checkDeviceStatus(const DeviceInformation& deviceInfo)
{
    hipError_t hip_err;

    // Is the generation of the device supported?
    if (deviceInfo.prop.major < 7)
    {
        return DeviceStatus::Incompatible;
    }

    hip_err = hipSetDevice(deviceInfo.id);
    if (hip_err != hipSuccess)
    {
        fprintf(stderr,
                "Error while switching to device #%d. %s\n",
                deviceInfo.id,
                gmx::getDeviceErrorString(hipGetLastError()).c_str());
        return DeviceStatus::NonFunctional;
    }

    /* try to execute a dummy kernel */
    dim3 gridSize(1, 1, 1);
    dim3 blockSize(512, 1, 1);

    hip_err = hipLaunchKernel(reinterpret_cast<void*>(dummy_kernel), gridSize, blockSize, nullptr);

    if (hip_err == hipErrorSharedObjectInitFailed || hip_err == hipErrorNoBinaryForGpu)
    {
        // Clear the error from attempting to compile and launch the kernel
        fprintf(stderr,
                "Error while trying to compile kernel, architecture %s of device #%d is not in the "
                "list of supported archs. Error: %s\n",
                deviceInfo.prop.gcnArchName,
                deviceInfo.id,
                gmx::getDeviceErrorString(hipGetLastError()).c_str());
        return DeviceStatus::DeviceNotTargeted;
    }
    else if (hip_err != hipSuccess)
    {
        // launchGpuKernel error is not fatal and should continue with marking the device bad
        fprintf(stderr,
                "Error occurred while running dummy kernel sanity check on device #%d:\n, hip "
                "error %s\n",
                deviceInfo.id,
                gmx::getDeviceErrorString(hipGetLastError()).c_str());
        return DeviceStatus::NonFunctional;
    }

    if (hipDeviceSynchronize() != hipSuccess)
    {
        fprintf(stderr,
                "Error while trying to synchronize for device #%d. %s\n",
                deviceInfo.id,
                gmx::getDeviceErrorString(hipGetLastError()).c_str());
        return DeviceStatus::NonFunctional;
    }

    // Skip context teardown when using HIP-aware MPI because this can lead to
    // corruption and a crash in MPI when when mdrunner is invoked multiple times
    // in the same process in gmxapi or mdrun integration tests. Ref #3952
    const bool haveDetectedOrForcedHipAwareMpi =
            (gmx::checkMpiHipAwareSupport() == gmx::GpuAwareMpiStatus::Supported
             || gmx::checkMpiHipAwareSupport() == gmx::GpuAwareMpiStatus::Forced);
    if (!haveDetectedOrForcedHipAwareMpi)
    {
        hip_err = hipDeviceReset();
        gmx::checkDeviceError(hipGetLastError(), "hipDeviceReset failed");
    }

    return DeviceStatus::Compatible;
}

bool isDeviceDetectionFunctional(std::string* errorMessage)
{
    int        numDevices;
    hipError_t stat = hipGetDeviceCount(&numDevices);
    if (stat != hipSuccess)
    {
        if (errorMessage != nullptr)
        {
            /* hipGetDeviceCount failed which means that there is
             * something wrong with the machine: driver-runtime
             * mismatch, all GPUs being busy in exclusive mode,
             * invalid HIP_VISIBLE_DEVICES, or some other condition
             * which should result in GROMACS issuing at least a
             * warning. */
            errorMessage->assign(hipGetErrorString(stat));
        }

        // Consume the error now that we have prepared to handle
        // it. This stops it reappearing next time we check for
        // errors.
        stat = hipGetLastError();
        // Can't detect GPUs
        return false;
    }

    // We don't actually use numDevices here, that's not the job of
    // this function.
    return true;
}

std::vector<std::unique_ptr<DeviceInformation>> findDevices()
{
    int        numDevices;
    hipError_t stat = hipGetDeviceCount(&numDevices);
    gmx::checkDeviceError(stat,
                          "Invalid call of findDevices() when HIP API returned an error, perhaps "
                          "canPerformDeviceDetection() was not called appropriately beforehand.");

    // We expect to start device support/sanity checks with a clean runtime error state
    gmx::ensureNoPendingDeviceError("Trying to find available HIP devices.");

    std::vector<std::unique_ptr<DeviceInformation>> deviceInfoList(numDevices);
    for (int i = 0; i < numDevices; i++)
    {
        hipDeviceProp_t prop;
        memset(&prop, 0, sizeof(hipDeviceProp_t));
        stat = hipGetDeviceProperties(&prop, i);

        deviceInfoList[i]               = std::make_unique<DeviceInformation>();
        deviceInfoList[i]->id           = i;
        deviceInfoList[i]->prop         = prop;
        deviceInfoList[i]->deviceVendor = DeviceVendor::Amd;

        deviceInfoList[i]->supportedSubGroupSizes.push_back(deviceInfoList[i]->prop.warpSize);

        const DeviceStatus checkResult = (stat != hipSuccess) ? DeviceStatus::NonFunctional
                                                              : checkDeviceStatus(*deviceInfoList[i]);

        deviceInfoList[i]->status = checkResult;

        if (checkResult != DeviceStatus::Compatible)
        {
            // TODO:
            //  - we inspect the HIP API state to retrieve and record any
            //    errors that occurred during is_gmx_supported_gpu_id() here,
            //    but this would be more elegant done within is_gmx_supported_gpu_id()
            //    and only return a string with the error if one was encountered.
            //  - we'll be reporting without rank information which is not ideal.
            //  - we'll end up warning also in cases where users would already
            //    get an error before mdrun aborts.
            //
            // Here we also clear the HIP API error state so potential
            // errors during sanity checks don't propagate.
            const std::string errorMessage = gmx::formatString(
                    "An error occurred while sanity checking device #%d.", deviceInfoList[i]->id);
            gmx::ensureNoPendingDeviceError(errorMessage);
        }
    }

    stat = hipPeekAtLastError();
    GMX_RELEASE_ASSERT(
            stat == hipSuccess,
            ("We promise to return with clean HIP state, but non-success state encountered. "
             + gmx::getDeviceErrorString(stat))
                    .c_str());

    return deviceInfoList;
}

void setActiveDevice(const DeviceInformation& deviceInfo)
{
    int        deviceId = deviceInfo.id;
    hipError_t stat;

    stat = hipSetDevice(deviceId);
    if (stat != hipSuccess)
    {
        auto message = gmx::formatString("Failed to initialize GPU #%d", deviceId);
        gmx::checkDeviceError(stat, message);
    }

    if (debug)
    {
        fprintf(stderr, "Initialized GPU ID #%d: %s\n", deviceId, deviceInfo.prop.name);
    }
}

void releaseDevice()
{
    hipError_t stat;

    int gpuid;
    stat = hipGetDevice(&gpuid);
    if (stat == hipSuccess)
    {
        if (debug)
        {
            fprintf(stderr, "Cleaning up context on GPU ID #%d.\n", gpuid);
        }

        stat = hipDeviceReset();
        if (stat != hipSuccess)
        {
            gmx_warning("Failed to free GPU #%d. %s", gpuid, gmx::getDeviceErrorString(stat).c_str());
        }
    }
}

std::string getDeviceInformationString(const DeviceInformation& deviceInfo)
{
    bool gpuExists = (deviceInfo.status != DeviceStatus::Nonexistent
                      && deviceInfo.status != DeviceStatus::NonFunctional);

    if (!gpuExists)
    {
        return gmx::formatString(
                "#%d: %s, stat: %s", deviceInfo.id, "N/A", c_deviceStateString[deviceInfo.status]);
    }
    else
    {
        return gmx::formatString("#%d: AMD %s, gcn: %s, ECC: %3s, stat: %s",
                                 deviceInfo.id,
                                 deviceInfo.prop.name,
                                 deviceInfo.prop.gcnArchName,
                                 deviceInfo.prop.ECCEnabled ? "yes" : " no",
                                 c_deviceStateString[deviceInfo.status]);
    }
}
