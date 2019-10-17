/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
 *  \brief Define functions for detection and initialization for OpenCL devices.
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 */

#include "gmxpre.h"

#include "config.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cstdio>
#ifdef __APPLE__
#    include <sys/sysctl.h>
#endif

#include <memory.h>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/ocl_compiler.h"
#include "gromacs/gpu_utils/oclraii.h"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

/*! \brief Return true if executing on compatible OS for AMD OpenCL.
 *
 * This is assumed to be true for OS X version of at least 10.10.4 and
 * all other OS flavors.
 *
 * Uses the BSD sysctl() interfaces to extract the kernel version.
 *
 * \return true if version is 14.4 or later (= OS X version 10.10.4),
 *         or OS is not Darwin.
 */
static bool runningOnCompatibleOSForAmd()
{
#ifdef __APPLE__
    int    mib[2];
    char   kernelVersion[256];
    size_t len = sizeof(kernelVersion);

    mib[0] = CTL_KERN;
    mib[1] = KERN_OSRELEASE;

    sysctl(mib, sizeof(mib) / sizeof(mib[0]), kernelVersion, &len, NULL, 0);

    int major = strtod(kernelVersion, NULL);
    int minor = strtod(strchr(kernelVersion, '.') + 1, NULL);

    // Kernel 14.4 corresponds to OS X 10.10.4
    return (major > 14 || (major == 14 && minor >= 4));
#else
    return true;
#endif
}

namespace gmx
{

/*! \brief Make an error string following an OpenCL API call.
 *
 *  It is meant to be called with \p status != CL_SUCCESS, but it will
 *  work correctly even if it is called with no OpenCL failure.
 *
 * \param[in]  message  Supplies context, e.g. the name of the API call that returned the error.
 * \param[in]  status   OpenCL API status code
 * \returns             A string describing the OpenCL error.
 */
static std::string makeOpenClInternalErrorString(const char* message, cl_int status)
{
    if (message != nullptr)
    {
        return formatString("%s did %ssucceed %d: %s", message, ((status != CL_SUCCESS) ? "not " : ""),
                            status, ocl_get_error_string(status).c_str());
    }
    else
    {
        return formatString("%sOpenCL error encountered %d: %s", ((status != CL_SUCCESS) ? "" : "No "),
                            status, ocl_get_error_string(status).c_str());
    }
}

/*!
 * \brief Checks that device \c devInfo is sane (ie can run a kernel).
 *
 * Compiles and runs a dummy kernel to determine whether the given
 * OpenCL device functions properly.
 *
 *
 * \param[in]  devInfo         The device info pointer.
 * \param[out] errorMessage    An error message related to a failing OpenCL API call.
 * \throws     std::bad_alloc  When out of memory.
 * \returns                    Whether the device passed sanity checks
 */
static bool isDeviceSane(const gmx_device_info_t* devInfo, std::string* errorMessage)
{
    cl_context_properties properties[] = {
        CL_CONTEXT_PLATFORM, reinterpret_cast<cl_context_properties>(devInfo->ocl_gpu_id.ocl_platform_id), 0
    };
    // uncrustify spacing

    cl_int    status;
    auto      deviceId = devInfo->ocl_gpu_id.ocl_device_id;
    ClContext context(clCreateContext(properties, 1, &deviceId, nullptr, nullptr, &status));
    if (status != CL_SUCCESS)
    {
        errorMessage->assign(makeOpenClInternalErrorString("clCreateContext", status));
        return false;
    }
    ClCommandQueue commandQueue(clCreateCommandQueue(context, deviceId, 0, &status));
    if (status != CL_SUCCESS)
    {
        errorMessage->assign(makeOpenClInternalErrorString("clCreateCommandQueue", status));
        return false;
    }

    // Some compilers such as Apple's require kernel functions to have at least one argument
    const char* lines[] = { "__kernel void dummyKernel(__global void* input){}" };
    ClProgram   program(clCreateProgramWithSource(context, 1, lines, nullptr, &status));
    if (status != CL_SUCCESS)
    {
        errorMessage->assign(makeOpenClInternalErrorString("clCreateProgramWithSource", status));
        return false;
    }

    if ((status = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr)) != CL_SUCCESS)
    {
        errorMessage->assign(makeOpenClInternalErrorString("clBuildProgram", status));
        return false;
    }

    ClKernel kernel(clCreateKernel(program, "dummyKernel", &status));
    if (status != CL_SUCCESS)
    {
        errorMessage->assign(makeOpenClInternalErrorString("clCreateKernel", status));
        return false;
    }

    clSetKernelArg(kernel, 0, sizeof(void*), nullptr);

    const size_t localWorkSize = 1, globalWorkSize = 1;
    if ((status = clEnqueueNDRangeKernel(commandQueue, kernel, 1, nullptr, &globalWorkSize,
                                         &localWorkSize, 0, nullptr, nullptr))
        != CL_SUCCESS)
    {
        errorMessage->assign(makeOpenClInternalErrorString("clEnqueueNDRangeKernel", status));
        return false;
    }
    return true;
}

/*!
 * \brief Checks that device \c devInfo is compatible with GROMACS.
 *
 *  Vendor and OpenCL version support checks are executed an the result
 *  of these returned.
 *
 * \param[in]  devInfo         The device info pointer.
 * \returns                    The result of the compatibility checks.
 */
static int isDeviceSupported(const gmx_device_info_t* devInfo)
{
    if (getenv("GMX_OCL_DISABLE_COMPATIBILITY_CHECK") != nullptr)
    {
        // Assume the device is compatible because checking has been disabled.
        return egpuCompatible;
    }

    // OpenCL device version check, ensure >= REQUIRED_OPENCL_MIN_VERSION
    constexpr unsigned int minVersionMajor = REQUIRED_OPENCL_MIN_VERSION_MAJOR;
    constexpr unsigned int minVersionMinor = REQUIRED_OPENCL_MIN_VERSION_MINOR;

    // Based on the OpenCL spec we're checking the version supported by
    // the device which has the following format:
    //      OpenCL<space><major_version.minor_version><space><vendor-specific information>
    unsigned int deviceVersionMinor, deviceVersionMajor;
    const int    valuesScanned = std::sscanf(devInfo->device_version, "OpenCL %u.%u",
                                          &deviceVersionMajor, &deviceVersionMinor);
    const bool   versionLargeEnough =
            ((valuesScanned == 2)
             && ((deviceVersionMajor > minVersionMajor)
                 || (deviceVersionMajor == minVersionMajor && deviceVersionMinor >= minVersionMinor)));
    if (!versionLargeEnough)
    {
        return egpuIncompatible;
    }

    /* Only AMD, Intel, and NVIDIA GPUs are supported for now */
    switch (devInfo->vendor_e)
    {
        case OCL_VENDOR_NVIDIA: return egpuCompatible;
        case OCL_VENDOR_AMD:
            return runningOnCompatibleOSForAmd() ? egpuCompatible : egpuIncompatible;
        case OCL_VENDOR_INTEL:
            return GMX_OPENCL_NB_CLUSTER_SIZE == 4 ? egpuCompatible : egpuIncompatibleClusterSize;
        default: return egpuIncompatible;
    }
}


/*! \brief Check whether the \c ocl_gpu_device is suitable for use by mdrun
 *
 * Runs sanity checks: checking that the runtime can compile a dummy kernel
 * and this can be executed;
 * Runs compatibility checks verifying the device OpenCL version requirement
 * and vendor/OS support.
 *
 * \param[in]  deviceId      The runtime-reported numeric ID of the device.
 * \param[in]  deviceInfo    The device info pointer.
 * \returns  An e_gpu_detect_res_t to indicate how the GPU coped with
 *           the sanity and compatibility check.
 */
static int checkGpu(size_t deviceId, const gmx_device_info_t* deviceInfo)
{

    int supportStatus = isDeviceSupported(deviceInfo);
    if (supportStatus != egpuCompatible)
    {
        return supportStatus;
    }

    std::string errorMessage;
    if (!isDeviceSane(deviceInfo, &errorMessage))
    {
        gmx_warning("While sanity checking device #%zu, %s", deviceId, errorMessage.c_str());
        return egpuInsane;
    }

    return egpuCompatible;
}

} // namespace gmx

/*! \brief Returns an ocl_vendor_id_t value corresponding to the input OpenCL vendor name.
 *
 *  \param[in] vendor_name String with OpenCL vendor name.
 *  \returns               ocl_vendor_id_t value for the input vendor_name
 */
static ocl_vendor_id_t get_vendor_id(char* vendor_name)
{
    if (vendor_name)
    {
        if (strstr(vendor_name, "NVIDIA"))
        {
            return OCL_VENDOR_NVIDIA;
        }
        else if (strstr(vendor_name, "AMD") || strstr(vendor_name, "Advanced Micro Devices"))
        {
            return OCL_VENDOR_AMD;
        }
        else if (strstr(vendor_name, "Intel"))
        {
            return OCL_VENDOR_INTEL;
        }
    }
    return OCL_VENDOR_UNKNOWN;
}

bool isGpuDetectionFunctional(std::string* errorMessage)
{
    cl_uint numPlatforms;
    cl_int  status = clGetPlatformIDs(0, nullptr, &numPlatforms);
    GMX_ASSERT(status != CL_INVALID_VALUE, "Incorrect call of clGetPlatformIDs detected");
#ifdef cl_khr_icd
    if (status == CL_PLATFORM_NOT_FOUND_KHR)
    {
        // No valid ICDs found
        if (errorMessage != nullptr)
        {
            errorMessage->assign("No valid OpenCL driver found");
        }
        return false;
    }
#endif
    GMX_RELEASE_ASSERT(
            status == CL_SUCCESS,
            gmx::formatString("An unexpected value was returned from clGetPlatformIDs %d: %s",
                              status, ocl_get_error_string(status).c_str())
                    .c_str());
    bool foundPlatform = (numPlatforms > 0);
    if (!foundPlatform && errorMessage != nullptr)
    {
        errorMessage->assign("No OpenCL platforms found even though the driver was valid");
    }
    return foundPlatform;
}

void findGpus(gmx_gpu_info_t* gpu_info)
{
    cl_uint         ocl_platform_count;
    cl_platform_id* ocl_platform_ids;
    cl_device_type  req_dev_type = CL_DEVICE_TYPE_GPU;

    ocl_platform_ids = nullptr;

    if (getenv("GMX_OCL_FORCE_CPU") != nullptr)
    {
        req_dev_type = CL_DEVICE_TYPE_CPU;
    }

    while (true)
    {
        cl_int status = clGetPlatformIDs(0, nullptr, &ocl_platform_count);
        if (CL_SUCCESS != status)
        {
            GMX_THROW(gmx::InternalError(
                    gmx::formatString("An unexpected value %d was returned from clGetPlatformIDs: ", status)
                    + ocl_get_error_string(status)));
        }

        if (1 > ocl_platform_count)
        {
            // TODO this should have a descriptive error message that we only support one OpenCL platform
            break;
        }

        snew(ocl_platform_ids, ocl_platform_count);

        status = clGetPlatformIDs(ocl_platform_count, ocl_platform_ids, nullptr);
        if (CL_SUCCESS != status)
        {
            GMX_THROW(gmx::InternalError(
                    gmx::formatString("An unexpected value %d was returned from clGetPlatformIDs: ", status)
                    + ocl_get_error_string(status)));
        }

        for (unsigned int i = 0; i < ocl_platform_count; i++)
        {
            cl_uint ocl_device_count;

            /* If requesting req_dev_type devices fails, just go to the next platform */
            if (CL_SUCCESS != clGetDeviceIDs(ocl_platform_ids[i], req_dev_type, 0, nullptr, &ocl_device_count))
            {
                continue;
            }

            if (1 <= ocl_device_count)
            {
                gpu_info->n_dev += ocl_device_count;
            }
        }

        if (1 > gpu_info->n_dev)
        {
            break;
        }

        snew(gpu_info->gpu_dev, gpu_info->n_dev);

        {
            int           device_index;
            cl_device_id* ocl_device_ids;

            snew(ocl_device_ids, gpu_info->n_dev);
            device_index = 0;

            for (unsigned int i = 0; i < ocl_platform_count; i++)
            {
                cl_uint ocl_device_count;

                /* If requesting req_dev_type devices fails, just go to the next platform */
                if (CL_SUCCESS
                    != clGetDeviceIDs(ocl_platform_ids[i], req_dev_type, gpu_info->n_dev,
                                      ocl_device_ids, &ocl_device_count))
                {
                    continue;
                }

                if (1 > ocl_device_count)
                {
                    break;
                }

                for (unsigned int j = 0; j < ocl_device_count; j++)
                {
                    gpu_info->gpu_dev[device_index].ocl_gpu_id.ocl_platform_id = ocl_platform_ids[i];
                    gpu_info->gpu_dev[device_index].ocl_gpu_id.ocl_device_id   = ocl_device_ids[j];

                    gpu_info->gpu_dev[device_index].device_name[0] = 0;
                    clGetDeviceInfo(ocl_device_ids[j], CL_DEVICE_NAME,
                                    sizeof(gpu_info->gpu_dev[device_index].device_name),
                                    gpu_info->gpu_dev[device_index].device_name, nullptr);

                    gpu_info->gpu_dev[device_index].device_version[0] = 0;
                    clGetDeviceInfo(ocl_device_ids[j], CL_DEVICE_VERSION,
                                    sizeof(gpu_info->gpu_dev[device_index].device_version),
                                    gpu_info->gpu_dev[device_index].device_version, nullptr);

                    gpu_info->gpu_dev[device_index].device_vendor[0] = 0;
                    clGetDeviceInfo(ocl_device_ids[j], CL_DEVICE_VENDOR,
                                    sizeof(gpu_info->gpu_dev[device_index].device_vendor),
                                    gpu_info->gpu_dev[device_index].device_vendor, nullptr);

                    gpu_info->gpu_dev[device_index].compute_units = 0;
                    clGetDeviceInfo(ocl_device_ids[j], CL_DEVICE_MAX_COMPUTE_UNITS,
                                    sizeof(gpu_info->gpu_dev[device_index].compute_units),
                                    &(gpu_info->gpu_dev[device_index].compute_units), nullptr);

                    gpu_info->gpu_dev[device_index].adress_bits = 0;
                    clGetDeviceInfo(ocl_device_ids[j], CL_DEVICE_ADDRESS_BITS,
                                    sizeof(gpu_info->gpu_dev[device_index].adress_bits),
                                    &(gpu_info->gpu_dev[device_index].adress_bits), nullptr);

                    gpu_info->gpu_dev[device_index].vendor_e =
                            get_vendor_id(gpu_info->gpu_dev[device_index].device_vendor);

                    clGetDeviceInfo(ocl_device_ids[j], CL_DEVICE_MAX_WORK_ITEM_SIZES, 3 * sizeof(size_t),
                                    &gpu_info->gpu_dev[device_index].maxWorkItemSizes, nullptr);

                    clGetDeviceInfo(ocl_device_ids[j], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t),
                                    &gpu_info->gpu_dev[device_index].maxWorkGroupSize, nullptr);

                    gpu_info->gpu_dev[device_index].stat =
                            gmx::checkGpu(device_index, gpu_info->gpu_dev + device_index);

                    if (egpuCompatible == gpu_info->gpu_dev[device_index].stat)
                    {
                        gpu_info->n_dev_compatible++;
                    }

                    device_index++;
                }
            }

            gpu_info->n_dev = device_index;

            /* Dummy sort of devices -  AMD first, then NVIDIA, then Intel */
            // TODO: Sort devices based on performance.
            if (0 < gpu_info->n_dev)
            {
                int last = -1;
                for (int i = 0; i < gpu_info->n_dev; i++)
                {
                    if (OCL_VENDOR_AMD == gpu_info->gpu_dev[i].vendor_e)
                    {
                        last++;

                        if (last < i)
                        {
                            gmx_device_info_t ocl_gpu_info;
                            ocl_gpu_info            = gpu_info->gpu_dev[i];
                            gpu_info->gpu_dev[i]    = gpu_info->gpu_dev[last];
                            gpu_info->gpu_dev[last] = ocl_gpu_info;
                        }
                    }
                }

                /* if more than 1 device left to be sorted */
                if ((gpu_info->n_dev - 1 - last) > 1)
                {
                    for (int i = 0; i < gpu_info->n_dev; i++)
                    {
                        if (OCL_VENDOR_NVIDIA == gpu_info->gpu_dev[i].vendor_e)
                        {
                            last++;

                            if (last < i)
                            {
                                gmx_device_info_t ocl_gpu_info;
                                ocl_gpu_info            = gpu_info->gpu_dev[i];
                                gpu_info->gpu_dev[i]    = gpu_info->gpu_dev[last];
                                gpu_info->gpu_dev[last] = ocl_gpu_info;
                            }
                        }
                    }
                }
            }

            sfree(ocl_device_ids);
        }

        break;
    }

    sfree(ocl_platform_ids);
}

void get_gpu_device_info_string(char* s, const gmx_gpu_info_t& gpu_info, int index)
{
    assert(s);

    if (index < 0 && index >= gpu_info.n_dev)
    {
        return;
    }

    gmx_device_info_t* dinfo = &gpu_info.gpu_dev[index];

    bool bGpuExists = (dinfo->stat != egpuNonexistent && dinfo->stat != egpuInsane);

    if (!bGpuExists)
    {
        sprintf(s, "#%d: %s, stat: %s", index, "N/A", gpu_detect_res_str[dinfo->stat]);
    }
    else
    {
        sprintf(s, "#%d: name: %s, vendor: %s, device version: %s, stat: %s", index, dinfo->device_name,
                dinfo->device_vendor, dinfo->device_version, gpu_detect_res_str[dinfo->stat]);
    }
}


void init_gpu(const gmx_device_info_t* deviceInfo)
{
    assert(deviceInfo);

    // If the device is NVIDIA, for safety reasons we disable the JIT
    // caching as this is known to be broken at least until driver 364.19;
    // the cache does not always get regenerated when the source code changes,
    // e.g. if the path to the kernel sources remains the same

    if (deviceInfo->vendor_e == OCL_VENDOR_NVIDIA)
    {
        // Ignore return values, failing to set the variable does not mean
        // that something will go wrong later.
#ifdef _MSC_VER
        _putenv("CUDA_CACHE_DISABLE=1");
#else
        // Don't override, maybe a dev is testing.
        setenv("CUDA_CACHE_DISABLE", "1", 0);
#endif
    }
}

gmx_device_info_t* getDeviceInfo(const gmx_gpu_info_t& gpu_info, int deviceId)
{
    if (deviceId < 0 || deviceId >= gpu_info.n_dev)
    {
        gmx_incons("Invalid GPU deviceId requested");
    }
    return &gpu_info.gpu_dev[deviceId];
}

size_t sizeof_gpu_dev_info()
{
    return sizeof(gmx_device_info_t);
}

int gpu_info_get_stat(const gmx_gpu_info_t& info, int index)
{
    return info.gpu_dev[index].stat;
}
