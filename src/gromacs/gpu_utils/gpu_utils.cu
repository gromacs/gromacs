/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "config.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <cuda_profiler_api.h>

#include "gromacs/gpu_utils/cudautils.cuh"
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

#if HAVE_NVML
#include <nvml.h>
#define HAVE_NVML_APPLICATION_CLOCKS (NVML_API_VERSION >= 6)
#else  /* HAVE_NVML */
#define HAVE_NVML_APPLICATION_CLOCKS 0
#endif /* HAVE_NVML */

#if defined(CHECK_CUDA_ERRORS) && HAVE_NVML_APPLICATION_CLOCKS
/*! Check for NVML error on the return status of a NVML API call. */
#  define HANDLE_NVML_RET_ERR(status, msg) \
    do { \
        if (status != NVML_SUCCESS) \
        { \
            gmx_warning("%s: %s\n", msg, nvmlErrorString(status)); \
        } \
    } while (0)
#else  /* defined(CHECK_CUDA_ERRORS) && HAVE_NVML_APPLICATION_CLOCKS */
#  define HANDLE_NVML_RET_ERR(status, msg) do { } while (0)
#endif /* defined(CHECK_CUDA_ERRORS) && HAVE_NVML_APPLICATION_CLOCKS */

#if HAVE_NVML_APPLICATION_CLOCKS
static const gmx_bool            bCompiledWithApplicationClockSupport = true;
#else
static const gmx_bool gmx_unused bCompiledWithApplicationClockSupport = false;
#endif

/*! \internal \brief
 * Max number of devices supported by CUDA (for consistency checking).
 *
 * In reality it is 16 with CUDA <=v5.0, but let's stay on the safe side.
 */
static int  cuda_max_device_count = 32;

static bool cudaProfilerRun      = ((getenv("NVPROF_ID") != NULL));

/** Dummy kernel used for sanity checking. */
static __global__ void k_dummy_test(void)
{
}

static void checkCompiledTargetCompatibility(const gmx_device_info_t *devInfo)
{
    assert(devInfo);

    cudaFuncAttributes attributes;
    cudaError_t        stat = cudaFuncGetAttributes(&attributes, k_dummy_test);

    if (cudaErrorInvalidDeviceFunction == stat)
    {
        gmx_fatal(FARGS,
                  "The %s binary does not include support for the CUDA architecture "
                  "of the selected GPU (device ID #%d, compute capability %d.%d). "
                  "By default, GROMACS supports all common architectures, so your GPU "
                  "might be rare, or some architectures were disabled in the build. ",
                  "Consult the install guide for how to use the GMX_CUDA_TARGET_SM and ",
                  "GMX_CUDA_TARGET_COMPUTE CMake variables to add this architecture.",
                  gmx::getProgramContext().displayName(), devInfo->id,
                  devInfo->prop.major, devInfo->prop.minor);
    }

    CU_RET_ERR(stat, "cudaFuncGetAttributes failed");

    if (devInfo->prop.major >= 3 && attributes.ptxVersion < 30)
    {
        gmx_fatal(FARGS,
                  "The GPU device code was compiled at runtime from 2.0 source which is "
                  "not compatible with the selected GPU (device ID #%d, compute capability %d.%d). "
                  "Pass the appropriate target in GMX_CUDA_TARGET_SM or a >=30 value to GMX_CUDA_TARGET_COMPUTE.",
                  devInfo->id,
                  devInfo->prop.major, devInfo->prop.minor);
    }
}

bool isHostMemoryPinned(void *h_ptr)
{
    cudaPointerAttributes memoryAttributes;
    cudaError_t           stat = cudaPointerGetAttributes(&memoryAttributes, h_ptr);

    bool                  result = false;
    switch (stat)
    {
        case cudaSuccess:
            result = true;
            break;

        case cudaErrorInvalidValue:
            // If the buffer was not pinned, then it will not be recognized by CUDA at all
            result = false;
            // Reset the last error status
            cudaGetLastError();
            break;

        default:
            CU_RET_ERR(stat, "Unexpected CUDA error");
    }
    return result;
}

/*!
 * \brief Runs GPU sanity checks.
 *
 * Runs a series of checks to determine that the given GPU and underlying CUDA
 * driver/runtime functions properly.
 * Returns properties of a device with given ID or the one that has
 * already been initialized earlier in the case if of \dev_id == -1.
 *
 * \param[in]  dev_id      the device ID of the GPU or -1 if the device has already been initialized
 * \param[out] dev_prop    pointer to the structure in which the device properties will be returned
 * \returns                0 if the device looks OK
 *
 * TODO: introduce errors codes and handle errors more smoothly.
 */
static int do_sanity_checks(int dev_id, cudaDeviceProp *dev_prop)
{
    cudaError_t cu_err;
    int         dev_count, id;

    cu_err = cudaGetDeviceCount(&dev_count);
    if (cu_err != cudaSuccess)
    {
        fprintf(stderr, "Error %d while querying device count: %s\n", cu_err,
                cudaGetErrorString(cu_err));
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
            fprintf(stderr, "Error %d while querying device id: %s\n", cu_err,
                    cudaGetErrorString(cu_err));
            return -1;
        }
    }
    else
    {
        id = dev_id;
        if (id > dev_count - 1) /* pfff there's no such device */
        {
            fprintf(stderr, "The requested device with id %d does not seem to exist (device count=%d)\n",
                    dev_id, dev_count);
            return -1;
        }
    }

    memset(dev_prop, 0, sizeof(cudaDeviceProp));
    cu_err = cudaGetDeviceProperties(dev_prop, id);
    if (cu_err != cudaSuccess)
    {
        fprintf(stderr, "Error %d while querying device properties: %s\n", cu_err,
                cudaGetErrorString(cu_err));
        return -1;
    }

    /* both major & minor is 9999 if no CUDA capable devices are present */
    if (dev_prop->major == 9999 && dev_prop->minor == 9999)
    {
        return -1;
    }
    /* we don't care about emulation mode */
    if (dev_prop->major == 0)
    {
        return -1;
    }

    if (id != -1)
    {
        cu_err = cudaSetDevice(id);
        if (cu_err != cudaSuccess)
        {
            fprintf(stderr, "Error %d while switching to device #%d: %s\n",
                    cu_err, id, cudaGetErrorString(cu_err));
            return -1;
        }
    }

    /* try to execute a dummy kernel */
    k_dummy_test<<< 1, 512>>> ();
    if (cudaThreadSynchronize() != cudaSuccess)
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

#if HAVE_NVML_APPLICATION_CLOCKS
/*! \brief Determines and adds the NVML device ID to the passed \cuda_dev.
 *
 * Determines and adds the NVML device ID to the passed \cuda_dev. This is done by
 * matching PCI-E information from \cuda_dev with the available NVML devices.
 *
 * \param[in,out] cuda_dev  CUDA device information to enrich with NVML device info
 * \returns                 true if \cuda_dev could be enriched with matching NVML device information.
 */
static bool addNVMLDeviceId(gmx_device_info_t* cuda_dev)
{
    nvmlDevice_t nvml_device_id;
    unsigned int nvml_device_count  = 0;
    nvmlReturn_t nvml_stat          = nvmlDeviceGetCount ( &nvml_device_count );
    bool         nvmlWasInitialized = false;
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetCount failed" );
    for (unsigned int nvml_device_idx = 0; nvml_stat == NVML_SUCCESS && nvml_device_idx < nvml_device_count; ++nvml_device_idx)
    {
        nvml_stat = nvmlDeviceGetHandleByIndex ( nvml_device_idx, &nvml_device_id );
        HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetHandleByIndex failed" );
        if (nvml_stat != NVML_SUCCESS)
        {
            break;
        }

        nvmlPciInfo_t nvml_pci_info;
        nvml_stat = nvmlDeviceGetPciInfo ( nvml_device_id, &nvml_pci_info );
        HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetPciInfo failed" );
        if (nvml_stat != NVML_SUCCESS)
        {
            break;
        }
        if (static_cast<unsigned int>(cuda_dev->prop.pciBusID) == nvml_pci_info.bus &&
            static_cast<unsigned int>(cuda_dev->prop.pciDeviceID) == nvml_pci_info.device &&
            static_cast<unsigned int>(cuda_dev->prop.pciDomainID) == nvml_pci_info.domain)
        {
            nvmlWasInitialized         = true;
            cuda_dev->nvml_device_id   = nvml_device_id;
            break;
        }
    }
    return nvmlWasInitialized;
}

/*! \brief Reads and returns the application clocks for device.
 *
 * \param[in]  device        The GPU device
 * \param[out] app_sm_clock  The current application SM clock
 * \param[out] app_mem_clock The current application memory clock
 * \returns if applacation clocks are supported
 */
static bool getApplicationClocks(const gmx_device_info_t *cuda_dev,
                                 unsigned int            *app_sm_clock,
                                 unsigned int            *app_mem_clock)
{
    nvmlReturn_t nvml_stat;

    nvml_stat = nvmlDeviceGetApplicationsClock(cuda_dev->nvml_device_id, NVML_CLOCK_SM, app_sm_clock);
    if (NVML_ERROR_NOT_SUPPORTED == nvml_stat)
    {
        return false;
    }
    HANDLE_NVML_RET_ERR(nvml_stat, "nvmlDeviceGetApplicationsClock failed for NVIDIA_CLOCK_SM");
    nvml_stat = nvmlDeviceGetApplicationsClock(cuda_dev->nvml_device_id, NVML_CLOCK_MEM, app_mem_clock);
    HANDLE_NVML_RET_ERR(nvml_stat, "nvmlDeviceGetApplicationsClock failed for NVIDIA_CLOCK_MEM");

    return true;
}
#endif /* HAVE_NVML_APPLICATION_CLOCKS */

/*! \brief Tries to set application clocks for the GPU with the given index.
 *
 * Application clocks are set to the max supported value to increase
 * performance if application clock permissions allow this. For future
 * GPU architectures a more sophisticated scheme might be required.
 *
 * \todo Refactor this into a detection phase and a work phase. Also
 * refactor to remove compile-time dependence on logging header.
 *
 * \param     mdlog         log file to write to
 * \param[in] cuda_dev      GPU device info for the GPU in use
 * \returns                 true if no error occurs during application clocks handling.
 */
static gmx_bool init_gpu_application_clocks(
        const gmx::MDLogger &mdlog,
        gmx_device_info_t   *cuda_dev)
{
    const cudaDeviceProp *prop                        = &cuda_dev->prop;
    int                   cuda_compute_capability     = prop->major * 10 + prop->minor;
    gmx_bool              bGpuCanUseApplicationClocks =
        ((0 == gmx_wcmatch("*Tesla*", prop->name) && cuda_compute_capability >= 35 ) ||
         (0 == gmx_wcmatch("*Quadro*", prop->name) && cuda_compute_capability >= 52 ));
    if (!bGpuCanUseApplicationClocks)
    {
        return true;
    }
#if !HAVE_NVML
    GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
            "NOTE: GROMACS was configured without NVML support hence it can not exploit\n"
            "      application clocks of the detected %s GPU to improve performance.\n"
            "      Recompile with the NVML library (compatible with the driver used) or set application clocks manually.",
            prop->name);
    return true;
#else
    if (!bCompiledWithApplicationClockSupport)
    {
        GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                "NOTE: GROMACS was compiled with an old NVML library which does not support\n"
                "      managing application clocks of the detected %s GPU to improve performance.\n"
                "      If your GPU supports application clocks, upgrade NVML (and driver) and recompile or set the clocks manually.",
                prop->name );
        return true;
    }

    /* We've compiled with NVML application clocks support, and have a GPU that can use it */
    nvmlReturn_t nvml_stat = NVML_SUCCESS;
    char        *env;
    //TODO: GMX_GPU_APPLICATION_CLOCKS is currently only used to enable/disable setting of application clocks
    //      this variable can be later used to give a user more fine grained control.
    env = getenv("GMX_GPU_APPLICATION_CLOCKS");
    if (env != NULL && ( strcmp( env, "0") == 0 ||
                         gmx_strcasecmp( env, "OFF") == 0 ||
                         gmx_strcasecmp( env, "DISABLE") == 0 ))
    {
        return true;
    }
    nvml_stat = nvmlInit();
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlInit failed." );
    if (nvml_stat != NVML_SUCCESS)
    {
        return false;
    }

    if (!addNVMLDeviceId(cuda_dev))
    {
        return false;
    }
    //get current application clocks setting
    if (!getApplicationClocks(cuda_dev,
                              &cuda_dev->nvml_orig_app_sm_clock,
                              &cuda_dev->nvml_orig_app_mem_clock))
    {
        return false;
    }
    //get max application clocks
    unsigned int max_sm_clock  = 0;
    unsigned int max_mem_clock = 0;
    nvml_stat = nvmlDeviceGetMaxClockInfo(cuda_dev->nvml_device_id, NVML_CLOCK_SM, &max_sm_clock);
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetMaxClockInfo failed" );
    nvml_stat = nvmlDeviceGetMaxClockInfo(cuda_dev->nvml_device_id, NVML_CLOCK_MEM, &max_mem_clock);
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetMaxClockInfo failed" );

    cuda_dev->nvml_is_restricted      = NVML_FEATURE_ENABLED;
    cuda_dev->nvml_app_clocks_changed = false;

    if (cuda_dev->nvml_orig_app_sm_clock >= max_sm_clock)
    {
        //TODO: This should probably be integrated into the GPU Properties table.
        GMX_LOG(mdlog.info).appendTextFormatted(
                "Application clocks (GPU clocks) for %s are (%d,%d)",
                cuda_dev->prop.name, cuda_dev->nvml_orig_app_mem_clock, cuda_dev->nvml_orig_app_sm_clock);
        return true;
    }

    if (cuda_compute_capability >= 60)
    {
        // Only warn about not being able to change clocks if they are not already at the max values
        if (max_mem_clock > cuda_dev->nvml_orig_app_mem_clock || max_sm_clock > cuda_dev->nvml_orig_app_sm_clock)
        {
            GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                    "Cannot change application clocks for %s to optimal values due to insufficient permissions. Current values are (%d,%d), max values are (%d,%d).\nPlease contact your admin to change application clocks.\n",
                    cuda_dev->prop.name, cuda_dev->nvml_orig_app_mem_clock, cuda_dev->nvml_orig_app_sm_clock, max_mem_clock, max_sm_clock);
        }
        return true;
    }

    nvml_stat = nvmlDeviceGetAPIRestriction(cuda_dev->nvml_device_id, NVML_RESTRICTED_API_SET_APPLICATION_CLOCKS, &(cuda_dev->nvml_is_restricted));
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetAPIRestriction failed" );

    if (nvml_stat != NVML_SUCCESS)
    {
        GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                "Cannot change GPU application clocks to optimal values due to NVML error (%d): %s.",
                nvml_stat, nvmlErrorString(nvml_stat));
        return false;
    }

    if (cuda_dev->nvml_is_restricted != NVML_FEATURE_DISABLED)
    {
        // Only warn about not being able to change clocks if they are not already at the max values
        if (max_mem_clock > cuda_dev->nvml_orig_app_mem_clock || max_sm_clock > cuda_dev->nvml_orig_app_sm_clock)
        {
            GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                    "Cannot change application clocks for %s to optimal values due to insufficient permissions. Current values are (%d,%d), max values are (%d,%d).\nUse sudo nvidia-smi -acp UNRESTRICTED or contact your admin to change application clocks.",
                    cuda_dev->prop.name, cuda_dev->nvml_orig_app_mem_clock, cuda_dev->nvml_orig_app_sm_clock, max_mem_clock, max_sm_clock);
        }
        return true;
    }

    /* Note: Distinguishing between different types of GPUs here might be necessary in the future,
       e.g. if max application clocks should not be used for certain GPUs. */
    GMX_LOG(mdlog.warning).appendTextFormatted(
            "Changing GPU application clocks for %s to (%d,%d)",
            cuda_dev->prop.name, max_mem_clock, max_sm_clock);
    nvml_stat = nvmlDeviceSetApplicationsClocks(cuda_dev->nvml_device_id, max_mem_clock, max_sm_clock);
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetApplicationsClock failed" );
    cuda_dev->nvml_app_clocks_changed = true;
    cuda_dev->nvml_set_app_sm_clock   = max_sm_clock;
    cuda_dev->nvml_set_app_mem_clock  = max_mem_clock;

    return true;
#endif /* HAVE_NVML */
}

/*! \brief Resets application clocks if changed and cleans up NVML for the passed \gpu_dev.
 *
 * \param[in] gpu_dev  CUDA device information
 */
static gmx_bool reset_gpu_application_clocks(const gmx_device_info_t gmx_unused * cuda_dev)
{
#if !HAVE_NVML_APPLICATION_CLOCKS
    GMX_UNUSED_VALUE(cuda_dev);
    return true;
#else /* HAVE_NVML_APPLICATION_CLOCKS */
    nvmlReturn_t nvml_stat = NVML_SUCCESS;
    if (cuda_dev &&
        cuda_dev->nvml_is_restricted == NVML_FEATURE_DISABLED &&
        cuda_dev->nvml_app_clocks_changed)
    {
        /* Check if the clocks are still what we set them to.
         * If so, set them back to the state we originally found them in.
         * If not, don't touch them, because something else set them later.
         */
        unsigned int app_sm_clock, app_mem_clock;
        getApplicationClocks(cuda_dev, &app_sm_clock, &app_mem_clock);
        if (app_sm_clock  == cuda_dev->nvml_set_app_sm_clock &&
            app_mem_clock == cuda_dev->nvml_set_app_mem_clock)
        {
            nvml_stat = nvmlDeviceSetApplicationsClocks(cuda_dev->nvml_device_id, cuda_dev->nvml_orig_app_mem_clock, cuda_dev->nvml_orig_app_sm_clock);
            HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceSetApplicationsClock failed" );
        }
    }
    nvml_stat = nvmlShutdown();
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlShutdown failed" );
    return (nvml_stat == NVML_SUCCESS);
#endif /* HAVE_NVML_APPLICATION_CLOCKS */
}

void init_gpu(const gmx::MDLogger &mdlog,
              gmx_device_info_t   *deviceInfo)
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

    checkCompiledTargetCompatibility(deviceInfo);

    //Ignoring return value as NVML errors should be treated not critical.
    init_gpu_application_clocks(mdlog, deviceInfo);
}

void free_gpu(const gmx_device_info_t *deviceInfo)
{
    // One should only attempt to clear the device context when
    // it has been used, but currently the only way to know that a GPU
    // device was used is that deviceInfo will be non-null.
    if (deviceInfo == nullptr)
    {
        return;
    }

    cudaError_t  stat;

    if (debug)
    {
        int gpuid;
        stat = cudaGetDevice(&gpuid);
        CU_RET_ERR(stat, "cudaGetDevice failed");
        fprintf(stderr, "Cleaning up context on GPU ID #%d\n", gpuid);
    }

    if (!reset_gpu_application_clocks(deviceInfo))
    {
        gmx_warning("Failed to reset GPU application clocks on GPU #%d", deviceInfo->id);
    }

    stat = cudaDeviceReset();
    if (stat != cudaSuccess)
    {
        gmx_warning("Failed to free GPU #%d: %s", deviceInfo->id, cudaGetErrorString(stat));
    }
}

gmx_device_info_t *getDeviceInfo(const gmx_gpu_info_t &gpu_info,
                                 int                   deviceId)
{
    if (deviceId < 0 || deviceId >= gpu_info.n_dev)
    {
        gmx_incons("Invalid GPU deviceId requested");
    }
    return &gpu_info.gpu_dev[deviceId];
}

/*! \brief Returns true if the gpu characterized by the device properties is
 *  supported by the native gpu acceleration.
 *
 * \param[in] dev_prop  the CUDA device properties of the gpus to test.
 * \returns             true if the GPU properties passed indicate a compatible
 *                      GPU, otherwise false.
 */
static bool is_gmx_supported_gpu(const cudaDeviceProp *dev_prop)
{
    return (dev_prop->major >= 2);
}

/*! \brief Checks if a GPU with a given ID is supported by the native GROMACS acceleration.
 *
 *  Returns a status value which indicates compatibility or one of the following
 *  errors: incompatibility, insistence, or insanity (=unexpected behavior).
 *  It also returns the respective device's properties in \dev_prop (if applicable).
 *
 *  As the error handling only permits returning the state of the GPU, this function
 *  does not clear the CUDA runtime API status allowing the caller to inspect the error
 *  upon return. Note that this also means it is the caller's responsibility to
 *  reset the CUDA runtime state.
 *
 *  \param[in]  dev_id   the ID of the GPU to check.
 *  \param[out] dev_prop the CUDA device properties of the device checked.
 *  \returns             the status of the requested device
 */
static int is_gmx_supported_gpu_id(int dev_id, cudaDeviceProp *dev_prop)
{
    cudaError_t stat;
    int         ndev;

    stat = cudaGetDeviceCount(&ndev);
    if (stat != cudaSuccess)
    {
        return egpuInsane;
    }

    if (dev_id > ndev - 1)
    {
        return egpuNonexistent;
    }

    /* TODO: currently we do not make a distinction between the type of errors
     * that can appear during sanity checks. This needs to be improved, e.g if
     * the dummy test kernel fails to execute with a "device busy message" we
     * should appropriately report that the device is busy instead of insane.
     */
    if (do_sanity_checks(dev_id, dev_prop) == 0)
    {
        if (is_gmx_supported_gpu(dev_prop))
        {
            return egpuCompatible;
        }
        else
        {
            return egpuIncompatible;
        }
    }
    else
    {
        return egpuInsane;
    }
}

bool canDetectGpus(std::string *errorMessage)
{
    cudaError_t        stat;
    int                driverVersion = -1;
    stat = cudaDriverGetVersion(&driverVersion);
    GMX_ASSERT(stat != cudaErrorInvalidValue, "An impossible null pointer was passed to cudaDriverGetVersion");
    GMX_RELEASE_ASSERT(stat == cudaSuccess,
                       gmx::formatString("An unexpected value was returned from cudaDriverGetVersion %s: %s",
                                         cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());
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

void findGpus(gmx_gpu_info_t *gpu_info)
{
    int                i, ndev, checkres;
    cudaError_t        stat;
    cudaDeviceProp     prop;
    gmx_device_info_t *devs;

    assert(gpu_info);

    gpu_info->n_dev_compatible = 0;

    ndev    = 0;
    devs    = NULL;

    stat = cudaGetDeviceCount(&ndev);
    if (stat != cudaSuccess)
    {
        GMX_THROW(gmx::InternalError("Invalid call of findGpus() when CUDA API returned an error, perhaps "
                                     "canDetectGpus() was not called appropriately beforehand."));
    }

    // We expect to start device support/sanity checks with a clean runtime error state
    gmx::ensureNoPendingCudaError("");

    snew(devs, ndev);
    for (i = 0; i < ndev; i++)
    {
        checkres = is_gmx_supported_gpu_id(i, &prop);

        devs[i].id   = i;
        devs[i].prop = prop;
        devs[i].stat = checkres;

        if (checkres == egpuCompatible)
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
                gmx_warning(gmx::formatString("An error occurred while sanity checking device #%d; %s: %s",
                                              devs[i].id, cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());
            }
        }
    }

    stat = cudaPeekAtLastError();
    GMX_RELEASE_ASSERT(stat == cudaSuccess,
                       gmx::formatString("We promise to return with clean CUDA state, but non-success state encountered: %s: %s",
                                         cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());

    gpu_info->n_dev   = ndev;
    gpu_info->gpu_dev = devs;
}

std::vector<int> getCompatibleGpus(const gmx_gpu_info_t &gpu_info)
{
    // Possible minor over-allocation here, but not important for anything
    std::vector<int> compatibleGpus;
    compatibleGpus.reserve(gpu_info.n_dev);
    for (int i = 0; i < gpu_info.n_dev; i++)
    {
        assert(gpu_info.gpu_dev);
        if (gpu_info.gpu_dev[i].stat == egpuCompatible)
        {
            compatibleGpus.push_back(i);
        }
    }
    return compatibleGpus;
}

const char *getGpuCompatibilityDescription(const gmx_gpu_info_t &gpu_info,
                                           int                   index)
{
    return (index >= gpu_info.n_dev ?
            gpu_detect_res_str[egpuNonexistent] :
            gpu_detect_res_str[gpu_info.gpu_dev[index].stat]);
}

void free_gpu_info(const gmx_gpu_info_t *gpu_info)
{
    if (gpu_info == NULL)
    {
        return;
    }

    sfree(gpu_info->gpu_dev);
}

void get_gpu_device_info_string(char *s, const gmx_gpu_info_t &gpu_info, int index)
{
    assert(s);

    if (index < 0 && index >= gpu_info.n_dev)
    {
        return;
    }

    gmx_device_info_t *dinfo = &gpu_info.gpu_dev[index];

    bool               bGpuExists =
        dinfo->stat == egpuCompatible ||
        dinfo->stat == egpuIncompatible;

    if (!bGpuExists)
    {
        sprintf(s, "#%d: %s, stat: %s",
                dinfo->id, "N/A",
                gpu_detect_res_str[dinfo->stat]);
    }
    else
    {
        sprintf(s, "#%d: NVIDIA %s, compute cap.: %d.%d, ECC: %3s, stat: %s",
                dinfo->id, dinfo->prop.name,
                dinfo->prop.major, dinfo->prop.minor,
                dinfo->prop.ECCEnabled ? "yes" : " no",
                gpu_detect_res_str[dinfo->stat]);
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
    return sizeof(gmx_device_info_t);
}

void gpu_set_host_malloc_and_free(bool               bUseGpuKernels,
                                  gmx_host_alloc_t **nb_alloc,
                                  gmx_host_free_t  **nb_free)
{
    if (bUseGpuKernels)
    {
        *nb_alloc = &pmalloc;
        *nb_free  = &pfree;
    }
    else
    {
        *nb_alloc = NULL;
        *nb_free  = NULL;
    }
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
