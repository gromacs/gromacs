/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef CUDAUTILS_CUH
#define CUDAUTILS_CUH

#include <stdio.h>

#include "gmx_fatal.h"

#include "cuda.h"


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
            gmx_warning("Just caught a previously occured CUDA error (%s), will try to continue.", cudaGetErrorString(_CU_CHECK_PREV_ERR_status)); \
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

#else

#define CU_RET_ERR(status, msg) do { } while (0)
#define CU_CHECK_PREV_ERR()     do { } while (0)
#define CU_LAUNCH_ERR(msg)      do { } while (0)
#define CU_LAUNCH_ERR_SYNC(msg) do { } while (0)

#endif /* CHECK_CUDA_ERRORS */ 

#ifdef __cplusplus
extern "C" {
#endif

/*! Device information: ID and properties of CUDA GPU use in the current process. */
typedef struct cu_dev_info
{
    int dev_id;                 /* id of the CUDA device in use */
    cudaDeviceProp dev_prop;    /* CUDA device properties */
} cu_dev_info_t;

/*! Launches asynchronous host to device memory copy in tstream 0. */
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
                         gmx_bool doAsync);

/*! Waits for event e to complete, */
int cu_wait_event(cudaEvent_t /*e*/);

/*! Caculates and returns the time ellapsed between event start and end. */
float cu_event_elapsed(cudaEvent_t /*start*/, cudaEvent_t /*end*/);

/*! Waits for event end to complete and calculates the time between start and end. */
int cu_wait_event_time(cudaEvent_t /*end*/, cudaEvent_t /*begin*/, float * /*time*/);

/*! Unbinds texture tex_name. */
void cu_unbind_texture(const char * /*tex_name*/);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/*! Binds texture tex_name to the GPU global memory pointed by d_ptr.*/
template <typename T>
size_t cu_bind_texture(const char * /*tex_name*/, const T * /*d_ptr*/, int /*size*/);
#endif


#endif /* CUDAUTILS_CUH */
