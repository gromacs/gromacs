/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015,2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_GPU_UTILS_CUDAUTILS_CUH
#define GMX_GPU_UTILS_CUDAUTILS_CUH

#include "config.h"

#include <stdio.h>
#if HAVE_NVML
#include <nvml.h>
#endif /* HAVE_NVML */

#include "gromacs/utility/fatalerror.h"

/* TODO error checking needs to be rewritten. We have 2 types of error checks needed
   based on where they occur in the code:
   - non performance-critical: these errors are unsafe to be ignored and must be
     _always_ checked for, e.g. initializations
   - performance critical: handling errors might hurt performance so care need to be taken
     when/if we should check for them at all, e.g. in cu_upload_X. However, we should be
     able to turn the check for these errors on!

   Probably we'll need two sets of the macros below...

 */
#define CHECK_CUDA_ERRORS

#ifdef CHECK_CUDA_ERRORS

/*! Check for CUDA error on the return status of a CUDA RT API call. */
#define CU_RET_ERR(status, msg) \
    do { \
        if (status != cudaSuccess) \
        { \
            gmx_fatal(FARGS, "%s: %s\n", msg, cudaGetErrorString(status)); \
        } \
    } while (0)

/*! Check for any previously occurred uncaught CUDA error. */
#define CU_CHECK_PREV_ERR() \
    do { \
        cudaError_t _CU_CHECK_PREV_ERR_status = cudaGetLastError(); \
        if (_CU_CHECK_PREV_ERR_status != cudaSuccess) { \
            gmx_warning("Just caught a previously occurred CUDA error (%s), will try to continue.", cudaGetErrorString(_CU_CHECK_PREV_ERR_status)); \
        } \
    } while (0)

/*! Check for any previously occurred uncaught CUDA error
   -- aimed at use after kernel calls. */
#define CU_LAUNCH_ERR(msg) \
    do { \
        cudaError_t _CU_LAUNCH_ERR_status = cudaGetLastError(); \
        if (_CU_LAUNCH_ERR_status != cudaSuccess) { \
            gmx_fatal(FARGS, "Error while launching kernel %s: %s\n", msg, cudaGetErrorString(_CU_LAUNCH_ERR_status)); \
        } \
    } while (0)

/*! Synchronize with GPU and check for any previously occurred uncaught CUDA error
   -- aimed at use after kernel calls. */
#define CU_LAUNCH_ERR_SYNC(msg) \
    do { \
        cudaError_t _CU_SYNC_LAUNCH_ERR_status = cudaThreadSynchronize(); \
        if (_CU_SYNC_LAUNCH_ERR_status != cudaSuccess) { \
            gmx_fatal(FARGS, "Error while launching kernel %s: %s\n", msg, cudaGetErrorString(_CU_SYNC_LAUNCH_ERR_status)); \
        } \
    } while (0)

#else /* CHECK_CUDA_ERRORS */

#define CU_RET_ERR(status, msg) do { } while (0)
#define CU_CHECK_PREV_ERR()     do { } while (0)
#define CU_LAUNCH_ERR(msg)      do { } while (0)
#define CU_LAUNCH_ERR_SYNC(msg) do { } while (0)
#define HANDLE_NVML_RET_ERR(status, msg) do { } while (0)

#endif /* CHECK_CUDA_ERRORS */

/*! \brief CUDA device information.
 *
 * The CUDA device information is queried and set at detection and contains
 * both information about the device/hardware returned by the runtime as well
 * as additional data like support status.
 *
 * \todo extract an object to manage NVML details
 */
struct gmx_device_info_t
{
    int                 id;                      /* id of the CUDA device */
    cudaDeviceProp      prop;                    /* CUDA device properties */
    int                 stat;                    /* result of the device check */
    unsigned int        nvml_orig_app_sm_clock;  /* The original SM clock before we changed it */
    unsigned int        nvml_orig_app_mem_clock; /* The original memory clock before we changed it */
    gmx_bool            nvml_app_clocks_changed; /* If application clocks have been changed */
    unsigned int        nvml_set_app_sm_clock;   /* The SM clock we set */
    unsigned int        nvml_set_app_mem_clock;  /* The memory clock we set */
#if HAVE_NVML
    nvmlDevice_t        nvml_device_id;          /* NVML device id */
    // TODO This can become a bool with a more useful name
    nvmlEnableState_t   nvml_is_restricted;      /* Status of application clocks permission */
#endif                                           /* HAVE_NVML */
};


/*! Launches asynchronous host to device memory copy in stream 0. */
int cu_copy_D2H(void * /*h_dest*/, void * /*d_src*/, size_t /*bytes*/);

/*! Launches asynchronous host to device memory copy in stream s. */
int cu_copy_D2H_async(void * /*h_dest*/, void * /*d_src*/, size_t /*bytes*/, cudaStream_t /*s = 0*/);

/*! Launches synchronous host to device memory copy. */
int cu_copy_H2D(void * /*d_dest*/, void * /*h_src*/, size_t /*bytes*/);

/*! Launches asynchronous host to device memory copy in stream s. */
int cu_copy_H2D_async(void * /*d_dest*/, void * /*h_src*/, size_t /*bytes*/, cudaStream_t /*s = 0*/);

/*! Frees device memory and resets the size and allocation size to -1. */
void cu_free_buffered(void *d_ptr, int *n = NULL, int *nalloc = NULL);

/*! Reallocates the device memory and copies data from the host. */
void cu_realloc_buffered(void **d_dest, void *h_src,
                         size_t type_size,
                         int *curr_size, int *curr_alloc_size,
                         int req_size,
                         cudaStream_t s,
                         bool bAsync);

/*! Waits for event e to complete, */
int cu_wait_event(cudaEvent_t /*e*/);

/*! Calculates and returns the time elapsed between event start and end. */
float cu_event_elapsed(cudaEvent_t /*start*/, cudaEvent_t /*end*/);

/*! Waits for event end to complete and calculates the time between start and end. */
int cu_wait_event_time(cudaEvent_t /*end*/, cudaEvent_t /*begin*/, float * /*time*/);

/*! \brief Return whether texture objects are used on this device.
 *
 * \todo This should be static in cudautils.cu, as soon as texture destruction code is moved there as well
 *
 * \param[in]   pointer to the GPU device info structure to inspect for texture objects support
 * \return      true if texture objects are used on this device
 */
bool use_texobj(const gmx_device_info_t *dev_info);

/*! \brief Initialize parameter lookup table.
 *
 * Initializes device memory, copies data from host and binds
 * a texture to allocated device memory to be used for parameter lookup.
 *
 * \tparam[in] T         Raw data type
 * \param[out] d_ptr     device pointer to the memory to be allocated
 * \param[out] texObj    texture object to be initialized
 * \param[out] texRef    texture reference to be initialized
 * \param[in]  h_ptr     pointer to the host memory to be uploaded to the device
 * \param[in]  numElem   number of elements in the h_ptr
 * \param[in]  devInfo   pointer to the info struct of the device in use
 */
template <typename T>
void initParamLookupTable(T                        * &d_ptr,
                          cudaTextureObject_t       &texObj,
                          const struct texture<T, 1, cudaReadModeElementType> *texRef,
                          const T                   *h_ptr,
                          int                        numElem,
                          const gmx_device_info_t   *devInfo);

#endif
