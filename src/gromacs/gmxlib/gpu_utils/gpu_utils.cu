/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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

#include "gromacs/gmxlib/cuda_tools/cudautils.cuh"
#include "gromacs/gmxlib/cuda_tools/pmalloc_cuda.h"
#include "gromacs/legacyheaders/types/hw_info.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

/*! \internal \brief
 * Max number of devices supported by CUDA (for consistency checking).
 *
 * In reality it is 16 with CUDA <=v5.0, but let's stay on the safe side.
 */
static int cuda_max_device_count = 32;

/** Dummy kernel used for sanity checking. */
__global__ void k_dummy_test()
{
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

#ifdef HAVE_NVML
/* TODO: We should actually be using md_print_warn in md_logging.c,
 * but we can't include mpi.h in CUDA code.
 */
static void md_print_info(FILE       *fplog,
                          const char *fmt, ...)
{
    va_list ap;

    if (fplog != NULL)
    {
        /* We should only print to stderr on the master node,
         * in most cases fplog is only set on the master node, so this works.
         */
        va_start(ap, fmt);
        vfprintf(stderr, fmt, ap);
        va_end(ap);

        va_start(ap, fmt);
        vfprintf(fplog, fmt, ap);
        va_end(ap);
    }
}
#endif /*HAVE_NVML*/

/* TODO: We should actually be using md_print_warn in md_logging.c,
 * but we can't include mpi.h in CUDA code.
 * This is replicated from nbnxn_cuda_data_mgmt.cu.
 */
static void md_print_warn(FILE       *fplog,
                          const char *fmt, ...)
{
    va_list ap;

    if (fplog != NULL)
    {
        /* We should only print to stderr on the master node,
         * in most cases fplog is only set on the master node, so this works.
         */
        va_start(ap, fmt);
        fprintf(stderr, "\n");
        vfprintf(stderr, fmt, ap);
        fprintf(stderr, "\n");
        va_end(ap);

        va_start(ap, fmt);
        fprintf(fplog, "\n");
        vfprintf(fplog, fmt, ap);
        fprintf(fplog, "\n");
        va_end(ap);
    }
}

#ifdef HAVE_NVML
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
    nvmlReturn_t nvml_stat = NVML_SUCCESS;
    nvmlDevice_t nvml_device_id;
    unsigned int nvml_device_count = 0;
    cuda_dev->nvml_initialized = false;
    nvml_stat                  = nvmlDeviceGetCount ( &nvml_device_count );
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
            cuda_dev->nvml_initialized = true;
            cuda_dev->nvml_device_id   = nvml_device_id;
            break;
        }
    }
    return cuda_dev->nvml_initialized;
}
#endif /*HAVE_NVML*/

/*! \brief Tries to set application clocks for the GPU with the given index.
 *
 * The variable \gpuid is the index of the GPU in the gpu_info.cuda_dev array
 * to handle the application clocks for. Application clocks are set to the
 * max supported value to increase performance if application clock permissions
 * allow this. For future GPU architectures a more sophisticated scheme might be
 * required.
 *
 * \param[out] fplog        log file to write to
 * \param[in] gpuid         index of the GPU to set application clocks for
 * \param[in] gpu_info      GPU info of all detected devices in the system.
 * \returns                 true if no error occurs during application clocks handling.
 */
static gmx_bool init_gpu_application_clocks(FILE gmx_unused *fplog, int gmx_unused gpuid, const gmx_gpu_info_t gmx_unused *gpu_info)
{
    const cudaDeviceProp *prop                        = &gpu_info->gpu_dev[gpuid].prop;
    int                   cuda_version_number         = prop->major * 10 + prop->minor;
    gmx_bool              bGpuCanUseApplicationClocks =
        ((0 == gmx_wcmatch("*Tesla*", prop->name) && cuda_version_number >= 35 ) ||
         (0 == gmx_wcmatch("*Quadro*", prop->name) && cuda_version_number >= 52 ));
#ifndef HAVE_NVML
    if (bGpuCanUseApplicationClocks)
    {
        int cuda_driver  = 0;
        int cuda_runtime = 0;
        cudaDriverGetVersion(&cuda_driver);
        cudaRuntimeGetVersion(&cuda_runtime);
        md_print_warn( fplog, "Note: NVML support was not found (CUDA runtime %d.%d, driver %d.%d), so your\n"
                       "      %s GPU cannot use application clock support to improve performance.\n",
                       cuda_runtime/1000, cuda_runtime%100,
                       cuda_driver/1000, cuda_driver%100,
                       prop->name );
    }
    return true;
#else /* HAVE_NVML defined */
    nvmlReturn_t nvml_stat = NVML_SUCCESS;
    char        *env;
    if (!bGpuCanUseApplicationClocks)
    {
        return true;
    }
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
    if (!addNVMLDeviceId( &(gpu_info->gpu_dev[gpuid])))
    {
        return false;
    }
    //get current application clocks setting
    unsigned int app_sm_clock  = 0;
    unsigned int app_mem_clock = 0;
    nvml_stat = nvmlDeviceGetApplicationsClock ( gpu_info->gpu_dev[gpuid].nvml_device_id, NVML_CLOCK_SM, &app_sm_clock );
    if (NVML_ERROR_NOT_SUPPORTED == nvml_stat)
    {
        return false;
    }
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetApplicationsClock failed" );
    nvml_stat = nvmlDeviceGetApplicationsClock ( gpu_info->gpu_dev[gpuid].nvml_device_id, NVML_CLOCK_MEM, &app_mem_clock );
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetApplicationsClock failed" );
    //get max application clocks
    unsigned int max_sm_clock  = 0;
    unsigned int max_mem_clock = 0;
    nvml_stat = nvmlDeviceGetMaxClockInfo ( gpu_info->gpu_dev[gpuid].nvml_device_id, NVML_CLOCK_SM, &max_sm_clock );
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetMaxClockInfo failed" );
    nvml_stat = nvmlDeviceGetMaxClockInfo ( gpu_info->gpu_dev[gpuid].nvml_device_id, NVML_CLOCK_MEM, &max_mem_clock );
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetMaxClockInfo failed" );

    gpu_info->gpu_dev[gpuid].nvml_is_restricted     = NVML_FEATURE_ENABLED;
    gpu_info->gpu_dev[gpuid].nvml_ap_clocks_changed = false;

    nvml_stat = nvmlDeviceGetAPIRestriction ( gpu_info->gpu_dev[gpuid].nvml_device_id, NVML_RESTRICTED_API_SET_APPLICATION_CLOCKS, &(gpu_info->gpu_dev[gpuid].nvml_is_restricted) );
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetAPIRestriction failed" );

    //TODO: Need to distinguish between different type of GPUs might be necessary in the future, e.g. if max application clocks should not be used
    //      for certain GPUs.
    if (nvml_stat == NVML_SUCCESS && app_sm_clock < max_sm_clock && gpu_info->gpu_dev[gpuid].nvml_is_restricted == NVML_FEATURE_DISABLED)
    {
        //TODO: Maybe need to think about something more user friendly here.
        md_print_info( fplog, "Changing GPU clock rates by setting application clocks for %s to (%d,%d)\n", gpu_info->gpu_dev[gpuid].prop.name, max_mem_clock, max_sm_clock);
        nvml_stat = nvmlDeviceSetApplicationsClocks ( gpu_info->gpu_dev[gpuid].nvml_device_id, max_mem_clock, max_sm_clock );
        HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceGetApplicationsClock failed" );
        gpu_info->gpu_dev[gpuid].nvml_ap_clocks_changed = true;
    }
    else if (nvml_stat == NVML_SUCCESS && app_sm_clock < max_sm_clock)
    {
        //TODO: Maybe need to think about something more user friendly here.
        md_print_warn( fplog, "Not possible to change GPU clocks to optimal value because of insufficient permissions to set application clocks for %s. Current values are (%d,%d). Max values are (%d,%d)\nUse sudo nvidia-smi -acp UNRESTRICTED or contact your admin to change application clock permissions.\n", gpu_info->gpu_dev[gpuid].prop.name, app_mem_clock, app_sm_clock, max_mem_clock, max_sm_clock);
    }
    else if (nvml_stat == NVML_SUCCESS && app_sm_clock == max_sm_clock)
    {
        //TODO: This should probably be integrated into the GPU Properties table.
        md_print_info( fplog, "Application clocks (GPU clocks) for %s are (%d,%d)\n", gpu_info->gpu_dev[gpuid].prop.name, app_mem_clock, app_sm_clock);
    }
    else
    {
        //TODO: Maybe need to think about something more user friendly here.
        md_print_warn( fplog,  "Not possible to change GPU clocks to optimal value because application clocks handling failed with NVML error (%d): %s.\n", nvml_stat, nvmlErrorString(nvml_stat));
    }
    return (nvml_stat == NVML_SUCCESS);
#endif /*HAVE_NVML*/
}

/*! \brief Resets application clocks if changed and cleans up NVML for the passed \gpu_dev.
 *
 * \param[in] gpu_dev  CUDA device information
 */
static gmx_bool reset_gpu_application_clocks(const gmx_device_info_t gmx_unused * cuda_dev)
{
#ifndef HAVE_NVML
    GMX_UNUSED_VALUE(cuda_dev);
    return true;
#else
    nvmlReturn_t nvml_stat = NVML_SUCCESS;
    if (cuda_dev &&
        cuda_dev->nvml_is_restricted == NVML_FEATURE_DISABLED &&
        cuda_dev->nvml_ap_clocks_changed)
    {
        nvml_stat = nvmlDeviceResetApplicationsClocks( cuda_dev->nvml_device_id );
        HANDLE_NVML_RET_ERR( nvml_stat, "nvmlDeviceResetApplicationsClocks failed" );
    }
    nvml_stat = nvmlShutdown();
    HANDLE_NVML_RET_ERR( nvml_stat, "nvmlShutdown failed" );
    return (nvml_stat == NVML_SUCCESS);
#endif /*HAVE_NVML*/
}

gmx_bool init_gpu(FILE gmx_unused *fplog, int mygpu, char *result_str,
                  const gmx_gpu_info_t *gpu_info,
                  const gmx_gpu_opt_t *gpu_opt)
{
    cudaError_t stat;
    char        sbuf[STRLEN];
    int         gpuid;

    assert(gpu_info);
    assert(result_str);

    if (mygpu < 0 || mygpu >= gpu_opt->n_dev_use)
    {
        sprintf(sbuf, "Trying to initialize an inexistent GPU: "
                "there are %d %s-selected GPU(s), but #%d was requested.",
                gpu_opt->n_dev_use, gpu_opt->bUserSet ? "user" : "auto", mygpu);
        gmx_incons(sbuf);
    }

    gpuid = gpu_info->gpu_dev[gpu_opt->dev_use[mygpu]].id;

    stat = cudaSetDevice(gpuid);
    strncpy(result_str, cudaGetErrorString(stat), STRLEN);

    if (debug)
    {
        fprintf(stderr, "Initialized GPU ID #%d: %s\n", gpuid, gpu_info->gpu_dev[gpuid].prop.name);
    }

    //Ignoring return value as NVML errors should be treated not critical.
    if (stat == cudaSuccess)
    {
        init_gpu_application_clocks(fplog, gpuid, gpu_info);
    }
    return (stat == cudaSuccess);
}

gmx_bool free_cuda_gpu(
        int gmx_unused mygpu, char *result_str,
        const gmx_gpu_info_t gmx_unused *gpu_info,
        const gmx_gpu_opt_t gmx_unused *gpu_opt
        )
{
    cudaError_t  stat;
    gmx_bool     reset_gpu_application_clocks_status = true;
    int          gpuid;

    assert(result_str);

    if (debug)
    {
        int gpuid;
        stat = cudaGetDevice(&gpuid);
        CU_RET_ERR(stat, "cudaGetDevice failed");
        fprintf(stderr, "Cleaning up context on GPU ID #%d\n", gpuid);
    }

    gpuid = gpu_opt ? gpu_opt->dev_use[mygpu] : -1;
    if (gpuid != -1)
    {
        reset_gpu_application_clocks_status = reset_gpu_application_clocks( &(gpu_info->gpu_dev[gpuid]) );
    }

    stat = cudaDeviceReset();
    strncpy(result_str, cudaGetErrorString(stat), STRLEN);
    return (stat == cudaSuccess) && reset_gpu_application_clocks_status;
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

/*! \brief Helper function that checks whether a given GPU status indicates compatible GPU.
 *
 * \param[in] stat  GPU status.
 * \returns         true if the provided status is egpuCompatible, otherwise false.
 */
static bool is_compatible_gpu(int stat)
{
    return (stat == egpuCompatible);
}

/*! \brief Checks if a GPU with a given ID is supported by the native GROMACS acceleration.
 *
 *  Returns a status value which indicates compatibility or one of the following
 *  errors: incompatibility, insistence, or insanity (=unexpected behavior).
 *  It also returns the respective device's properties in \dev_prop (if applicable).
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


int detect_gpus(gmx_gpu_info_t *gpu_info, char *err_str)
{
    int                i, ndev, checkres, retval;
    cudaError_t        stat;
    cudaDeviceProp     prop;
    gmx_device_info_t *devs;

    assert(gpu_info);
    assert(err_str);

    gpu_info->n_dev_compatible = 0;

    ndev    = 0;
    devs    = NULL;

    stat = cudaGetDeviceCount(&ndev);
    if (stat != cudaSuccess)
    {
        const char *s;

        /* cudaGetDeviceCount failed which means that there is something
         * wrong with the machine: driver-runtime mismatch, all GPUs being
         * busy in exclusive mode, or some other condition which should
         * result in us issuing a warning a falling back to CPUs. */
        retval = -1;
        s      = cudaGetErrorString(stat);
        strncpy(err_str, s, STRLEN*sizeof(err_str[0]));
    }
    else
    {
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
        }
        retval = 0;
    }

    gpu_info->n_dev   = ndev;
    gpu_info->gpu_dev = devs;

    return retval;
}

void pick_compatible_gpus(const gmx_gpu_info_t *gpu_info,
                          gmx_gpu_opt_t        *gpu_opt)
{
    int  i, ncompat;
    int *compat;

    assert(gpu_info);
    /* gpu_dev/n_dev have to be either NULL/0 or not (NULL/0) */
    assert((gpu_info->n_dev != 0 ? 0 : 1) ^ (gpu_info->gpu_dev == NULL ? 0 : 1));

    snew(compat, gpu_info->n_dev);
    ncompat = 0;
    for (i = 0; i < gpu_info->n_dev; i++)
    {
        if (is_compatible_gpu(gpu_info->gpu_dev[i].stat))
        {
            ncompat++;
            compat[ncompat - 1] = i;
        }
    }

    gpu_opt->n_dev_compatible = ncompat;
    snew(gpu_opt->dev_compatible, ncompat);
    memcpy(gpu_opt->dev_compatible, compat, ncompat*sizeof(*compat));
    sfree(compat);
}

gmx_bool check_selected_gpus(int                  *checkres,
                             const gmx_gpu_info_t *gpu_info,
                             gmx_gpu_opt_t        *gpu_opt)
{
    int  i, id;
    bool bAllOk;

    assert(checkres);
    assert(gpu_info);
    assert(gpu_opt->n_dev_use >= 0);

    if (gpu_opt->n_dev_use == 0)
    {
        return TRUE;
    }

    assert(gpu_opt->dev_use);

    /* we will assume that all GPUs requested are valid IDs,
       otherwise we'll bail anyways */

    bAllOk = true;
    for (i = 0; i < gpu_opt->n_dev_use; i++)
    {
        id = gpu_opt->dev_use[i];

        /* devices are stored in increasing order of IDs in gpu_dev */
        gpu_opt->dev_use[i] = id;

        checkres[i] = (id >= gpu_info->n_dev) ?
            egpuNonexistent : gpu_info->gpu_dev[id].stat;

        bAllOk = bAllOk && is_compatible_gpu(checkres[i]);
    }

    return bAllOk;
}

void free_gpu_info(const gmx_gpu_info_t *gpu_info)
{
    if (gpu_info == NULL)
    {
        return;
    }

    sfree(gpu_info->gpu_dev);
}

void get_gpu_device_info_string(char *s, const gmx_gpu_info_t *gpu_info, int index)
{
    assert(s);
    assert(gpu_info);

    if (index < 0 && index >= gpu_info->n_dev)
    {
        return;
    }

    gmx_device_info_t *dinfo = &gpu_info->gpu_dev[index];

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

int get_gpu_device_id(const gmx_gpu_info_t *gpu_info,
                      const gmx_gpu_opt_t  *gpu_opt,
                      int                   idx)
{
    assert(gpu_info);
    assert(gpu_opt);
    assert(idx >= 0 && idx < gpu_opt->n_dev_use);

    return gpu_info->gpu_dev[gpu_opt->dev_use[idx]].id;
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
