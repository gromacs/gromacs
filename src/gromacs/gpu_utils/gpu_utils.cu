/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015,2016, The GROMACS development team.
 * Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
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
/*! \file
 *  \brief Define functions for detection and initialization for CUDA devices.
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 */

#include "gmxpre.h"

#include "gpu_utils.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <cuda_profiler_api.h>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/pmalloc_cuda.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/stringutil.h"

/*! \internal \brief
 * Max number of devices supported by CUDA (for consistency checking).
 *
 * In reality it is 16 with CUDA <=v5.0, but let's stay on the safe side.
 */
static int cuda_max_device_count = 32;

static bool cudaProfilerRun = ((getenv("NVPROF_ID") != nullptr));

/** Dummy kernel used for sanity checking. */
static __global__ void k_dummy_test(void) {}

static cudaError_t checkCompiledTargetCompatibility(int deviceId, const cudaDeviceProp& deviceProp)
{
    cudaFuncAttributes attributes;
    cudaError_t        stat = cudaFuncGetAttributes(&attributes, k_dummy_test);

    if (cudaErrorInvalidDeviceFunction == stat)
    {
        fprintf(stderr,
                "\nWARNING: The %s binary does not include support for the CUDA architecture of "
                "the GPU ID #%d (compute capability %d.%d) detected during detection. "
                "By default, GROMACS supports all architectures of compute "
                "capability >= 3.0, so your GPU "
                "might be rare, or some architectures were disabled in the build. \n"
                "Consult the install guide for how to use the GMX_CUDA_TARGET_SM and "
                "GMX_CUDA_TARGET_COMPUTE CMake variables to add this architecture. \n",
                gmx::getProgramContext().displayName(), deviceId, deviceProp.major, deviceProp.minor);
    }

    return stat;
}

bool isHostMemoryPinned(const void* h_ptr)
{
    cudaPointerAttributes memoryAttributes;
    cudaError_t           stat = cudaPointerGetAttributes(&memoryAttributes, h_ptr);

    bool result = false;
    switch (stat)
    {
        case cudaSuccess: result = true; break;

        case cudaErrorInvalidValue:
            // If the buffer was not pinned, then it will not be recognized by CUDA at all
            result = false;
            // Reset the last error status
            cudaGetLastError();
            break;

        default: CU_RET_ERR(stat, "Unexpected CUDA error");
    }
    return result;
}

/*!
 * \brief Runs GPU sanity checks.
 *
 * Runs a series of checks to determine that the given GPU and underlying CUDA
 * driver/runtime functions properly.
 *
 * \param[in]  dev_id      the device ID of the GPU or -1 if the device has already been initialized
 * \param[in]  dev_prop    The device properties structure
 * \returns                0 if the device looks OK, -1 if it sanity checks failed, and -2 if the device is busy
 *
 * TODO: introduce errors codes and handle errors more smoothly.
 */
static int do_sanity_checks(int dev_id, const cudaDeviceProp& dev_prop)
{
    cudaError_t cu_err;
    int         dev_count, id;

    cu_err = cudaGetDeviceCount(&dev_count);
    if (cu_err != cudaSuccess)
    {
        fprintf(stderr, "Error %d while querying device count: %s\n", cu_err, cudaGetErrorString(cu_err));
        return -1;
    }

    /* no CUDA compatible device at all */
    if (dev_count == 0)
    {
        return -1;
    }

    /* things might go horribly wrong if cudart is not compatible with the driver */
    if (dev_count < 0 || dev_count > cuda_max_device_count)
    {
        return -1;
    }

    if (dev_id == -1) /* device already selected let's not destroy the context */
    {
        cu_err = cudaGetDevice(&id);
        if (cu_err != cudaSuccess)
        {
            fprintf(stderr, "Error %d while querying device id: %s\n", cu_err, cudaGetErrorString(cu_err));
            return -1;
        }
    }
    else
    {
        id = dev_id;
        if (id > dev_count - 1) /* pfff there's no such device */
        {
            fprintf(stderr,
                    "The requested device with id %d does not seem to exist (device count=%d)\n",
                    dev_id, dev_count);
            return -1;
        }
    }

    /* both major & minor is 9999 if no CUDA capable devices are present */
    if (dev_prop.major == 9999 && dev_prop.minor == 9999)
    {
        return -1;
    }
    /* we don't care about emulation mode */
    if (dev_prop.major == 0)
    {
        return -1;
    }

    if (id != -1)
    {
        cu_err = cudaSetDevice(id);
        if (cu_err != cudaSuccess)
        {
            fprintf(stderr, "Error %d while switching to device #%d: %s\n", cu_err, id,
                    cudaGetErrorString(cu_err));
            return -1;
        }
    }

    cu_err = checkCompiledTargetCompatibility(dev_id, dev_prop);
    // Avoid triggering an error if GPU devices are in exclusive or prohibited mode;
    // it is enough to check for cudaErrorDevicesUnavailable only here because
    // if we encounter it that will happen in cudaFuncGetAttributes in the above function.
    if (cu_err == cudaErrorDevicesUnavailable)
    {
        return -2;
    }
    else if (cu_err != cudaSuccess)
    {
        return -1;
    }

    /* try to execute a dummy kernel */
    try
    {
        KernelLaunchConfig config;
        config.blockSize[0]                = 512;
        const auto          dummyArguments = prepareGpuKernelArguments(k_dummy_test, config);
        DeviceInformation   deviceInfo;
        const DeviceContext deviceContext(deviceInfo);
        const DeviceStream  deviceStream(deviceContext, DeviceStreamPriority::Normal, false);
        launchGpuKernel(k_dummy_test, config, deviceStream, nullptr, "Dummy kernel", dummyArguments);
    }
    catch (gmx::GromacsException& ex)
    {
        // launchGpuKernel error is not fatal and should continue with marking the device bad
        fprintf(stderr,
                "Error occurred while running dummy kernel sanity check on device #%d:\n %s\n", id,
                formatExceptionMessageToString(ex).c_str());
        return -1;
    }

    if (cudaDeviceSynchronize() != cudaSuccess)
    {
        return -1;
    }

    /* destroy context if we created one */
    if (id != -1)
    {
        cu_err = cudaDeviceReset();
        CU_RET_ERR(cu_err, "cudaDeviceReset failed");
    }

    return 0;
}

void init_gpu(const DeviceInformation* deviceInfo)
{
    cudaError_t stat;

    assert(deviceInfo);

    stat = cudaSetDevice(deviceInfo->id);
    if (stat != cudaSuccess)
    {
        auto message = gmx::formatString("Failed to initialize GPU #%d", deviceInfo->id);
        CU_RET_ERR(stat, message.c_str());
    }

    if (debug)
    {
        fprintf(stderr, "Initialized GPU ID #%d: %s\n", deviceInfo->id, deviceInfo->prop.name);
    }
}

void free_gpu(const DeviceInformation* deviceInfo)
{
    // One should only attempt to clear the device context when
    // it has been used, but currently the only way to know that a GPU
    // device was used is that deviceInfo will be non-null.
    if (deviceInfo == nullptr)
    {
        return;
    }

    cudaError_t stat;

    if (debug)
    {
        int gpuid;
        stat = cudaGetDevice(&gpuid);
        CU_RET_ERR(stat, "cudaGetDevice failed");
        fprintf(stderr, "Cleaning up context on GPU ID #%d\n", gpuid);
    }

    stat = cudaDeviceReset();
    if (stat != cudaSuccess)
    {
        gmx_warning("Failed to free GPU #%d: %s", deviceInfo->id, cudaGetErrorString(stat));
    }
}

DeviceInformation* getDeviceInfo(const gmx_gpu_info_t& gpu_info, int deviceId)
{
    if (deviceId < 0 || deviceId >= gpu_info.n_dev)
    {
        gmx_incons("Invalid GPU deviceId requested");
    }
    return &gpu_info.deviceInfo[deviceId];
}

/*! \brief Returns true if the gpu characterized by the device properties is
 *  supported by the native gpu acceleration.
 *
 * \param[in] dev_prop  the CUDA device properties of the gpus to test.
 * \returns             true if the GPU properties passed indicate a compatible
 *                      GPU, otherwise false.
 */
static bool is_gmx_supported_gpu(const cudaDeviceProp& dev_prop)
{
    return (dev_prop.major >= 3);
}

/*! \brief Checks if a GPU with a given ID is supported by the native GROMACS acceleration.
 *
 *  Returns a status value which indicates compatibility or one of the following
 *  errors: incompatibility or insanity (=unexpected behavior).
 *
 *  As the error handling only permits returning the state of the GPU, this function
 *  does not clear the CUDA runtime API status allowing the caller to inspect the error
 *  upon return. Note that this also means it is the caller's responsibility to
 *  reset the CUDA runtime state.
 *
 *  \param[in]  deviceId   the ID of the GPU to check.
 *  \param[in]  deviceProp the CUDA device properties of the device checked.
 *  \returns               the status of the requested device
 */
static int is_gmx_supported_gpu_id(int deviceId, const cudaDeviceProp& deviceProp)
{
    if (!is_gmx_supported_gpu(deviceProp))
    {
        return egpuIncompatible;
    }

    /* TODO: currently we do not make a distinction between the type of errors
     * that can appear during sanity checks. This needs to be improved, e.g if
     * the dummy test kernel fails to execute with a "device busy message" we
     * should appropriately report that the device is busy instead of insane.
     */
    const int checkResult = do_sanity_checks(deviceId, deviceProp);
    switch (checkResult)
    {
        case 0: return egpuCompatible;
        case -1: return egpuInsane;
        case -2: return egpuUnavailable;
        default:
            GMX_RELEASE_ASSERT(false, "Invalid do_sanity_checks() return value");
            return egpuCompatible;
    }
}

bool isGpuDetectionFunctional(std::string* errorMessage)
{
    cudaError_t stat;
    int         driverVersion = -1;
    stat                      = cudaDriverGetVersion(&driverVersion);
    GMX_ASSERT(stat != cudaErrorInvalidValue,
               "An impossible null pointer was passed to cudaDriverGetVersion");
    GMX_RELEASE_ASSERT(
            stat == cudaSuccess,
            gmx::formatString("An unexpected value was returned from cudaDriverGetVersion %s: %s",
                              cudaGetErrorName(stat), cudaGetErrorString(stat))
                    .c_str());
    bool foundDriver = (driverVersion > 0);
    if (!foundDriver)
    {
        // Can't detect GPUs if there is no driver
        if (errorMessage != nullptr)
        {
            errorMessage->assign("No valid CUDA driver found");
        }
        return false;
    }

    int numDevices;
    stat = cudaGetDeviceCount(&numDevices);
    if (stat != cudaSuccess)
    {
        if (errorMessage != nullptr)
        {
            /* cudaGetDeviceCount failed which means that there is
             * something wrong with the machine: driver-runtime
             * mismatch, all GPUs being busy in exclusive mode,
             * invalid CUDA_VISIBLE_DEVICES, or some other condition
             * which should result in GROMACS issuing at least a
             * warning. */
            errorMessage->assign(cudaGetErrorString(stat));
        }

        // Consume the error now that we have prepared to handle
        // it. This stops it reappearing next time we check for
        // errors. Note that if CUDA_VISIBLE_DEVICES does not contain
        // valid devices, then cudaGetLastError returns the
        // (undocumented) cudaErrorNoDevice, but this should not be a
        // problem as there should be no future CUDA API calls.
        // NVIDIA bug report #2038718 has been filed.
        cudaGetLastError();
        // Can't detect GPUs
        return false;
    }

    // We don't actually use numDevices here, that's not the job of
    // this function.
    return true;
}

void findGpus(gmx_gpu_info_t* gpu_info)
{
    assert(gpu_info);

    gpu_info->n_dev_compatible = 0;

    int         ndev;
    cudaError_t stat = cudaGetDeviceCount(&ndev);
    if (stat != cudaSuccess)
    {
        GMX_THROW(gmx::InternalError(
                "Invalid call of findGpus() when CUDA API returned an error, perhaps "
                "canDetectGpus() was not called appropriately beforehand."));
    }

    // We expect to start device support/sanity checks with a clean runtime error state
    gmx::ensureNoPendingCudaError("");

    DeviceInformation* devs;
    snew(devs, ndev);
    for (int i = 0; i < ndev; i++)
    {
        cudaDeviceProp prop;
        memset(&prop, 0, sizeof(cudaDeviceProp));
        stat = cudaGetDeviceProperties(&prop, i);
        int checkResult;
        if (stat != cudaSuccess)
        {
            // Will handle the error reporting below
            checkResult = egpuInsane;
        }
        else
        {
            checkResult = is_gmx_supported_gpu_id(i, prop);
        }

        devs[i].id   = i;
        devs[i].prop = prop;
        devs[i].stat = checkResult;

        if (checkResult == egpuCompatible)
        {
            gpu_info->n_dev_compatible++;
        }
        else
        {
            // TODO:
            //  - we inspect the CUDA API state to retrieve and record any
            //    errors that occurred during is_gmx_supported_gpu_id() here,
            //    but this would be more elegant done within is_gmx_supported_gpu_id()
            //    and only return a string with the error if one was encountered.
            //  - we'll be reporting without rank information which is not ideal.
            //  - we'll end up warning also in cases where users would already
            //    get an error before mdrun aborts.
            //
            // Here we also clear the CUDA API error state so potential
            // errors during sanity checks don't propagate.
            if ((stat = cudaGetLastError()) != cudaSuccess)
            {
                gmx_warning("An error occurred while sanity checking device #%d; %s: %s",
                            devs[i].id, cudaGetErrorName(stat), cudaGetErrorString(stat));
            }
        }
    }

    stat = cudaPeekAtLastError();
    GMX_RELEASE_ASSERT(stat == cudaSuccess,
                       gmx::formatString("We promise to return with clean CUDA state, but "
                                         "non-success state encountered: %s: %s",
                                         cudaGetErrorName(stat), cudaGetErrorString(stat))
                               .c_str());

    gpu_info->n_dev      = ndev;
    gpu_info->deviceInfo = devs;
}

void get_gpu_device_info_string(char* s, const gmx_gpu_info_t& gpu_info, int index)
{
    assert(s);

    if (index < 0 && index >= gpu_info.n_dev)
    {
        return;
    }

    DeviceInformation* dinfo = &gpu_info.deviceInfo[index];

    bool bGpuExists = (dinfo->stat != egpuNonexistent && dinfo->stat != egpuInsane);

    if (!bGpuExists)
    {
        sprintf(s, "#%d: %s, stat: %s", dinfo->id, "N/A", gpu_detect_res_str[dinfo->stat]);
    }
    else
    {
        sprintf(s, "#%d: NVIDIA %s, compute cap.: %d.%d, ECC: %3s, stat: %s", dinfo->id,
                dinfo->prop.name, dinfo->prop.major, dinfo->prop.minor,
                dinfo->prop.ECCEnabled ? "yes" : " no", gpu_detect_res_str[dinfo->stat]);
    }
}

int get_current_cuda_gpu_device_id(void)
{
    int gpuid;
    CU_RET_ERR(cudaGetDevice(&gpuid), "cudaGetDevice failed");

    return gpuid;
}

size_t sizeof_gpu_dev_info(void)
{
    return sizeof(DeviceInformation);
}

void startGpuProfiler(void)
{
    /* The NVPROF_ID environment variable is set by nvprof and indicates that
       mdrun is executed in the CUDA profiler.
       If nvprof was run is with "--profile-from-start off", the profiler will
       be started here. This way we can avoid tracing the CUDA events from the
       first part of the run. Starting the profiler again does nothing.
     */
    if (cudaProfilerRun)
    {
        cudaError_t stat;
        stat = cudaProfilerStart();
        CU_RET_ERR(stat, "cudaProfilerStart failed");
    }
}

void stopGpuProfiler(void)
{
    /* Stopping the nvidia here allows us to eliminate the subsequent
       API calls from the trace, e.g. uninitialization and cleanup. */
    if (cudaProfilerRun)
    {
        cudaError_t stat;
        stat = cudaProfilerStop();
        CU_RET_ERR(stat, "cudaProfilerStop failed");
    }
}

void resetGpuProfiler(void)
{
    /* With CUDA <=7.5 the profiler can't be properly reset; we can only start
     *  the profiling here (can't stop it) which will achieve the desired effect if
     *  the run was started with the profiling disabled.
     *
     * TODO: add a stop (or replace it with reset) when this will work correctly in CUDA.
     * stopGpuProfiler();
     */
    if (cudaProfilerRun)
    {
        startGpuProfiler();
    }
}

int gpu_info_get_stat(const gmx_gpu_info_t& info, int index)
{
    return info.deviceInfo[index].stat;
}

/*! \brief Check status returned from peer access CUDA call, and error out or warn appropriately
 * \param[in] stat           CUDA call return status
 * \param[in] gpuA           ID for GPU initiating peer access call
 * \param[in] gpuB           ID for remote GPU
 * \param[in] mdlog          Logger object
 * \param[in] cudaCallName   name of CUDA peer access call
 */
static void peerAccessCheckStat(const cudaError_t    stat,
                                const int            gpuA,
                                const int            gpuB,
                                const gmx::MDLogger& mdlog,
                                const char*          cudaCallName)
{
    if ((stat == cudaErrorInvalidDevice) || (stat == cudaErrorInvalidValue))
    {
        std::string errorString =
                gmx::formatString("%s from GPU %d to GPU %d failed", cudaCallName, gpuA, gpuB);
        CU_RET_ERR(stat, errorString.c_str());
    }
    if (stat != cudaSuccess)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendTextFormatted(
                        "GPU peer access not enabled between GPUs %d and %d due to unexpected "
                        "return value from %s: %s",
                        gpuA, gpuB, cudaCallName, cudaGetErrorString(stat));
    }
}

void setupGpuDevicePeerAccess(const std::vector<int>& gpuIdsToUse, const gmx::MDLogger& mdlog)
{
    cudaError_t stat;

    // take a note of currently-set GPU
    int currentGpu;
    stat = cudaGetDevice(&currentGpu);
    CU_RET_ERR(stat, "cudaGetDevice in setupGpuDevicePeerAccess failed");

    std::string message = gmx::formatString(
            "Note: Peer access enabled between the following GPU pairs in the node:\n ");
    bool peerAccessEnabled = false;

    for (unsigned int i = 0; i < gpuIdsToUse.size(); i++)
    {
        int gpuA = gpuIdsToUse[i];
        stat     = cudaSetDevice(gpuA);
        if (stat != cudaSuccess)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "GPU peer access not enabled due to unexpected return value from "
                            "cudaSetDevice(%d): %s",
                            gpuA, cudaGetErrorString(stat));
            return;
        }
        for (unsigned int j = 0; j < gpuIdsToUse.size(); j++)
        {
            if (j != i)
            {
                int gpuB          = gpuIdsToUse[j];
                int canAccessPeer = 0;
                stat              = cudaDeviceCanAccessPeer(&canAccessPeer, gpuA, gpuB);
                peerAccessCheckStat(stat, gpuA, gpuB, mdlog, "cudaDeviceCanAccessPeer");

                if (canAccessPeer)
                {
                    stat = cudaDeviceEnablePeerAccess(gpuB, 0);
                    peerAccessCheckStat(stat, gpuA, gpuB, mdlog, "cudaDeviceEnablePeerAccess");

                    message           = gmx::formatString("%s%d->%d ", message.c_str(), gpuA, gpuB);
                    peerAccessEnabled = true;
                }
            }
        }
    }

    // re-set GPU to that originally set
    stat = cudaSetDevice(currentGpu);
    if (stat != cudaSuccess)
    {
        CU_RET_ERR(stat, "cudaSetDevice in setupGpuDevicePeerAccess failed");
        return;
    }

    if (peerAccessEnabled)
    {
        GMX_LOG(mdlog.info).asParagraph().appendTextFormatted("%s", message.c_str());
    }
}
