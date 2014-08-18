/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "gromacs/legacyheaders/types/hw_info.h"

#include "gromacs/legacyheaders/gpu_utils.h"
#include "../cuda_tools/cudautils.cuh"

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

/*! \brief
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
#if CUDA_VERSION < 4000
        cu_err = cudaThreadExit();
        CU_RET_ERR(cu_err, "cudaThreadExit failed");
#else
        cu_err = cudaDeviceReset();
        CU_RET_ERR(cu_err, "cudaDeviceReset failed");
#endif
    }

    return 0;
}

/*! \brief Initializes the GPU with the given index.
 *
 * The varible \mygpu is the index of the GPU to initialize in the
 * gpu_info.cuda_dev array.
 *
 * \param[in]  mygpu        index of the GPU to initialize
 * \param[out] result_str   the message related to the error that occurred
 *                          during the initialization (if there was any).
 * \param[in] gpu_info      GPU info of all detected devices in the system.
 * \param[in] gpu_opt       options for using the GPUs in gpu_info
 * \returns                 true if no error occurs during initialization.
 */
gmx_bool init_gpu(int mygpu, char *result_str,
                  const gmx_gpu_info_t *gpu_info,
                  const gmx_gpu_opt_t *gpu_opt)
{
    cudaError_t stat;
    char        sbuf[STRLEN];
    int         gpuid;

    assert(gpu_info);
    assert(result_str);

    if (mygpu < 0 || mygpu >= gpu_opt->ncuda_dev_use)
    {
        sprintf(sbuf, "Trying to initialize an inexistent GPU: "
                "there are %d %s-selected GPU(s), but #%d was requested.",
                gpu_opt->ncuda_dev_use, gpu_opt->bUserSet ? "user" : "auto", mygpu);
        gmx_incons(sbuf);
    }

    gpuid = gpu_info->cuda_dev[gpu_opt->cuda_dev_use[mygpu]].id;

    stat = cudaSetDevice(gpuid);
    strncpy(result_str, cudaGetErrorString(stat), STRLEN);

    if (debug)
    {
        fprintf(stderr, "Initialized GPU ID #%d: %s\n", gpuid, gpu_info->cuda_dev[gpuid].prop.name);
    }

    return (stat == cudaSuccess);
}

/*! \brief Frees up the CUDA GPU used by the active context at the time of calling.
 *
 * The context is explicitly destroyed and therefore all data uploaded to the GPU
 * is lost. This should only be called when none of this data is required anymore.
 *
 * \param[out] result_str   the message related to the error that occurred
 *                          during the initialization (if there was any).
 * \returns                 true if no error occurs during the freeing.
 */
gmx_bool free_gpu(char *result_str)
{
    cudaError_t stat;

    assert(result_str);

    if (debug)
    {
        int gpuid;
        stat = cudaGetDevice(&gpuid);
        CU_RET_ERR(stat, "cudaGetDevice failed");
        fprintf(stderr, "Cleaning up context on GPU ID #%d\n", gpuid);
    }

#if CUDA_VERSION < 4000
    stat = cudaThreadExit();
#else
    stat = cudaDeviceReset();
#endif
    strncpy(result_str, cudaGetErrorString(stat), STRLEN);

    return (stat == cudaSuccess);
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


/*! \brief Detect all NVIDIA GPUs in the system.
 *
 *  Will detect every NVIDIA GPU supported by the device driver in use. Also
 *  check for the compatibility of each and fill the gpu_info->cuda_dev array
 *  with the required information on each the device: ID, device properties,
 *  status.
 *
 *  \param[in] gpu_info    pointer to structure holding GPU information.
 *  \param[out] err_str    The error message of any CUDA API error that caused
 *                         the detection to fail (if there was any). The memory
 *                         the pointer points to should be managed externally.
 *  \returns               non-zero if the detection encountered a failure, zero otherwise.
 */
int detect_cuda_gpus(gmx_gpu_info_t *gpu_info, char *err_str)
{
    int              i, ndev, checkres, retval;
    cudaError_t      stat;
    cudaDeviceProp   prop;
    cuda_dev_info_t *devs;

    assert(gpu_info);
    assert(err_str);

    gpu_info->ncuda_dev_compatible = 0;

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
                gpu_info->ncuda_dev_compatible++;
            }
        }
        retval = 0;
    }

    gpu_info->ncuda_dev = ndev;
    gpu_info->cuda_dev  = devs;

    return retval;
}

/*! \brief Select the GPUs compatible with the native GROMACS acceleration.
 *
 * This function selects the compatible gpus and initializes
 * gpu_info->cuda_dev_use and gpu_info->ncuda_dev_use.
 *
 * Given the list of GPUs available in the system check each device in
 * gpu_info->cuda_dev and place the indices of the compatible GPUs into
 * cuda_dev_use with this marking the respective GPUs as "available for use."
 * Note that \detect_cuda_gpus must have been called before.
 *
 * \param[in]     gpu_info    pointer to structure holding GPU information
 * \param[in,out] gpu_opt     pointer to structure holding GPU options
 */
void pick_compatible_gpus(const gmx_gpu_info_t *gpu_info,
                          gmx_gpu_opt_t        *gpu_opt)
{
    int  i, ncompat;
    int *compat;

    assert(gpu_info);
    /* cuda_dev/ncuda_dev have to be either NULL/0 or not (NULL/0) */
    assert((gpu_info->ncuda_dev != 0 ? 0 : 1) ^ (gpu_info->cuda_dev == NULL ? 0 : 1));

    snew(compat, gpu_info->ncuda_dev);
    ncompat = 0;
    for (i = 0; i < gpu_info->ncuda_dev; i++)
    {
        if (is_compatible_gpu(gpu_info->cuda_dev[i].stat))
        {
            ncompat++;
            compat[ncompat - 1] = i;
        }
    }

    gpu_opt->ncuda_dev_use = ncompat;
    snew(gpu_opt->cuda_dev_use, ncompat);
    memcpy(gpu_opt->cuda_dev_use, compat, ncompat*sizeof(*compat));
    sfree(compat);
}

/*! \brief Check the existence/compatibility of a set of GPUs specified by their device IDs.
 *
 * Given the a list of gpu->ncuda_dev_use GPU device IDs stored in
 * gpu_opt->cuda_dev_use check the existence and compatibility
 * of the respective GPUs. Also provide the caller with an array containing
 * the result of checks in \checkres.
 *
 * \param[out]  checkres    check result for each ID passed in \requested_devs
 * \param[in]   gpu_info    pointer to structure holding GPU information
 * \param[out]  gpu_opt     pointer to structure holding GPU options
 * \returns                 TRUE if every the requested GPUs are compatible
 */
gmx_bool check_selected_cuda_gpus(int                  *checkres,
                                  const gmx_gpu_info_t *gpu_info,
                                  gmx_gpu_opt_t        *gpu_opt)
{
    int  i, id;
    bool bAllOk;

    assert(checkres);
    assert(gpu_info);
    assert(gpu_opt->ncuda_dev_use >= 0);

    if (gpu_opt->ncuda_dev_use == 0)
    {
        return TRUE;
    }

    assert(gpu_opt->cuda_dev_use);

    /* we will assume that all GPUs requested are valid IDs,
       otherwise we'll bail anyways */

    bAllOk = true;
    for (i = 0; i < gpu_opt->ncuda_dev_use; i++)
    {
        id = gpu_opt->cuda_dev_use[i];

        /* devices are stored in increasing order of IDs in cuda_dev */
        gpu_opt->cuda_dev_use[i] = id;

        checkres[i] = (id >= gpu_info->ncuda_dev) ?
            egpuNonexistent : gpu_info->cuda_dev[id].stat;

        bAllOk = bAllOk && is_compatible_gpu(checkres[i]);
    }

    return bAllOk;
}

/*! \brief Frees the cuda_dev and cuda_dev_use array fields of \gpu_info.
 *
 * \param[in]    gpu_info    pointer to structure holding GPU information
 */
void free_gpu_info(const gmx_gpu_info_t *gpu_info)
{
    if (gpu_info == NULL)
    {
        return;
    }

    sfree(gpu_info->cuda_dev);
}

/*! \brief Formats and returns a device information string for a given GPU.
 *
 * Given an index *directly* into the array of available GPUs (cuda_dev)
 * returns a formatted info string for the respective GPU which includes
 * ID, name, compute capability, and detection status.
 *
 * \param[out]  s           pointer to output string (has to be allocated externally)
 * \param[in]   gpu_info    pointer to structure holding GPU information
 * \param[in]   index       an index *directly* into the array of available GPUs
 */
void get_gpu_device_info_string(char *s, const gmx_gpu_info_t *gpu_info, int index)
{
    assert(s);
    assert(gpu_info);

    if (index < 0 && index >= gpu_info->ncuda_dev)
    {
        return;
    }

    cuda_dev_info_t *dinfo = &gpu_info->cuda_dev[index];

    bool             bGpuExists =
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

/*! \brief Returns the device ID of the GPU with a given index into the array of used GPUs.
 *
 * Getter function which, given an index into the array of GPUs in use
 * (cuda_dev_use) -- typically a tMPI/MPI rank --, returns the device ID of the
 * respective CUDA GPU.
 *
 * \param[in]    gpu_info   pointer to structure holding GPU information
 * \param[in]    gpu_opt    pointer to structure holding GPU options
 * \param[in]    idx        index into the array of used GPUs
 * \returns                 device ID of the requested GPU
 */
int get_gpu_device_id(const gmx_gpu_info_t *gpu_info,
                      const gmx_gpu_opt_t  *gpu_opt,
                      int                   idx)
{
    assert(gpu_info);
    assert(gpu_opt);
    assert(idx >= 0 && idx < gpu_opt->ncuda_dev_use);

    return gpu_info->cuda_dev[gpu_opt->cuda_dev_use[idx]].id;
}

/*! \brief Returns the device ID of the GPU currently in use.
 *
 * The GPU used is the one that is active at the time of the call in the active context.
 *
 * \param[in]    gpu_info   pointer to structure holding GPU information
 * \returns                 device ID of the GPU in use at the time of the call
 */
int get_current_gpu_device_id(void)
{
    int gpuid;
    CU_RET_ERR(cudaGetDevice(&gpuid), "cudaGetDevice failed");

    return gpuid;
}

/*! \brief Returns the size of the cuda_dev_info struct.
 *
 * The size of cuda_dev_info can be used for allocation and communication.
 *
 * \returns                 size in bytes of cuda_dev_info
 */
size_t sizeof_cuda_dev_info(void)
{
    return sizeof(cuda_dev_info);
}
