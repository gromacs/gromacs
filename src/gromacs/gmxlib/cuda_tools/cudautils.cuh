/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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

#ifndef CUDAUTILS_CUH
#define CUDAUTILS_CUH

#include "config.h"

#include <stdio.h>
#ifdef HAVE_NVML
#include <nvml.h>
#endif /* HAVE_NVML */

#include "gromacs/utility/fatalerror.h"

/* CUDA library and hardware related defines */
/* TODO list some constants instead that can be used for consistency checks to
   detect future devices with features that make the currect code incompatible
   with them (e.g. expected warp size = 32, check against the dev_info->props.warpsize). */
#define WARP_SIZE           32

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

/*! Check for NVML error on the return status of a NVML API call. */
#ifdef HAVE_NVML
#define HANDLE_NVML_RET_ERR(status, msg) \
    do { \
        if (status != NVML_SUCCESS) \
        { \
            gmx_warning("%s: %s\n", msg, nvmlErrorString(status)); \
        } \
    } while (0)
#endif /* HAVE_NVML */
#else

#define CU_RET_ERR(status, msg) do { } while (0)
#define CU_CHECK_PREV_ERR()     do { } while (0)
#define CU_LAUNCH_ERR(msg)      do { } while (0)
#define CU_LAUNCH_ERR_SYNC(msg) do { } while (0)
#define HANDLE_NVML_RET_ERR(status, msg) do { } while (0)

#endif /* CHECK_CUDA_ERRORS */

#ifdef __cplusplus
extern "C" {
#endif

/*! CUDA device information. */
struct gmx_device_info_t
{
    int                 id;                     /* id of the CUDA device */
    cudaDeviceProp      prop;                   /* CUDA device properties */
    int                 stat;                   /* result of the device check */
    gmx_bool            nvml_initialized;       /* If NVML was initialized */
    gmx_bool            nvml_ap_clocks_changed; /* If application clocks have been changed */
#ifdef HAVE_NVML
    nvmlDevice_t        nvml_device_id;         /* NVML device id */
    nvmlEnableState_t   nvml_is_restricted;     /* Status of application clocks permission */
#endif                                          /* HAVE_NVML */
};


/*! Launches asynchronous host to device memory copy in stream 0. */
int cu_copy_D2H(void * /*h_dest*/, void * /*d_src*/, size_t /*bytes*/);

/*! Launches asynchronous host to device memory copy in stream s. */
int cu_copy_D2H_async(void * /*h_dest*/, void * /*d_src*/, size_t /*bytes*/, cudaStream_t /*s = 0*/);

/*! Allocates host memory and launches synchronous host to device memory copy. */
int cu_copy_D2H_alloc(void ** /*h_dest*/, void * /*d_src*/, size_t /*bytes*/);


/*! Launches synchronous host to device memory copy. */
int cu_copy_H2D(void * /*d_dest*/, void * /*h_src*/, size_t /*bytes*/);

/*! Launches asynchronous host to device memory copy in stream s. */
int cu_copy_H2D_async(void * /*d_dest*/, void * /*h_src*/, size_t /*bytes*/, cudaStream_t /*s = 0*/);

/*! Allocates device memory and launches synchronous host to device memory copy. */
int cu_copy_H2D_alloc(void ** /*d_dest*/, void * /*h_src*/, size_t /*bytes*/);

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

#ifdef __cplusplus
}
#endif

#endif /* CUDAUTILS_CUH */
