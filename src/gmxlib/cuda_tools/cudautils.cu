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

#include <stdlib.h>

#include "gmx_fatal.h"
#include "smalloc.h"

#include "cuda.h"
#include "cudautils.cuh"

/*** Generic CUDA data operation wrappers ***/

/*! Launches synchronous or asynchronous host to device memory copy.
 *
 *  The copy is launched in stream s or if not spefied, in stream 0.
 */
static int cu_copy_D2H_generic(void * h_dest, void * d_src, size_t bytes, 
                               gmx_bool async = FALSE, cudaStream_t s = 0)
{
    cudaError_t stat;
    
    if (h_dest == 0 || d_src == 0 || bytes <= 0)
        return -1;

    if (async)
    {
        stat = cudaMemcpyAsync(h_dest, d_src, bytes, cudaMemcpyDeviceToHost, s);
        CU_RET_ERR(stat, "DtoH cudaMemcpyAsync failed");

    }
    else
    {
        stat = cudaMemcpy(h_dest, d_src, bytes, cudaMemcpyDeviceToHost);
        CU_RET_ERR(stat, "DtoH cudaMemcpy failed");
    }

    return 0;
}

int cu_copy_D2H(void * h_dest, void * d_src, size_t bytes)
{
    return cu_copy_D2H_generic(h_dest, d_src, bytes, FALSE);
}

/*!
 *  The copy is launched in stream s or if not spefied, in stream 0.
 */
int cu_copy_D2H_async(void * h_dest, void * d_src, size_t bytes, cudaStream_t s = 0)
{
    return cu_copy_D2H_generic(h_dest, d_src, bytes, TRUE, s);
}

int cu_copy_D2H_alloc(void ** h_dest, void * d_src, size_t bytes)
{ 
    if (h_dest == 0 || d_src == 0 || bytes <= 0)
        return -1;

    smalloc(*h_dest, bytes);

    return cu_copy_D2H(*h_dest, d_src, bytes);
}

/*! Launches synchronous or asynchronous device to host memory copy.
 *
 *  The copy is launched in stream s or if not spefied, in stream 0.
 */
static int cu_copy_H2D_generic(void * d_dest, void * h_src, size_t bytes, 
                               gmx_bool async = FALSE, cudaStream_t s = 0)
{
    cudaError_t stat;

    if (d_dest == 0 || h_src == 0 || bytes <= 0)
        return -1;

    if (async)
    {
        stat = cudaMemcpyAsync(d_dest, h_src, bytes, cudaMemcpyHostToDevice, s);
        CU_RET_ERR(stat, "HtoD cudaMemcpyAsync failed");
    }
    else
    {
        stat = cudaMemcpy(d_dest, h_src, bytes, cudaMemcpyHostToDevice);
        CU_RET_ERR(stat, "HtoD cudaMemcpy failed");
    }

    return 0;
}

int cu_copy_H2D(void * d_dest, void * h_src, size_t bytes)
{   
    return cu_copy_H2D_generic(d_dest, h_src, bytes, FALSE);
}

/*!
 *  The copy is launched in stream s or if not spefied, in stream 0.
 */
int cu_copy_H2D_async(void * d_dest, void * h_src, size_t bytes, cudaStream_t s = 0)
{   
    return cu_copy_H2D_generic(d_dest, h_src, bytes, TRUE, s);
}

int cu_copy_H2D_alloc(void ** d_dest, void * h_src, size_t bytes)
{
    cudaError_t stat;

    if (d_dest == 0 || h_src == 0 || bytes <= 0)
        return -1;

    stat = cudaMalloc(d_dest, bytes);
    CU_RET_ERR(stat, "cudaMalloc failed in cu_copy_H2D_alloc");

    return cu_copy_H2D(*d_dest, h_src, bytes);
}

float cu_event_elapsed(cudaEvent_t start, cudaEvent_t end)
{
    float t = 0.0;
    cudaError_t stat;

    stat = cudaEventElapsedTime(&t, start, end);
    CU_RET_ERR(stat, "cudaEventElapsedTime failed in cu_event_elapsed");

    return t;
}

int cu_wait_event(cudaEvent_t e)
{
    cudaError_t s;

    s = cudaEventSynchronize(e);
    CU_RET_ERR(s, "cudaEventSynchronize failed in cu_wait_event");

    return 0;
}

/*! 
 *  If time != NULL it also calculates the time ellapsed between start and end and 
 *  return this is milliseconds.
 */ 
int cu_wait_event_time(cudaEvent_t end, cudaEvent_t start, float *time)
{
    cudaError_t s;

    s = cudaEventSynchronize(end);
    CU_RET_ERR(s, "cudaEventSynchronize failed in cu_wait_event");

    if (time)
    {
        *time = cu_event_elapsed(start, end);
    }

    return 0;
}

/*! 
 *  The name of the texture is passed in text_name, while its size in size.
 *  Returns the offset that needs to be used when fetching from the texture.
 */
template <typename T>
size_t cu_bind_texture(const char *tex_name, const T *d_ptr, int size)
{
    cudaError_t             stat;
    cudaChannelFormatDesc   cd;
    const textureReference  *tex;
    char                    str[100];

    size_t offset;

    stat = cudaGetTextureReference(&tex, tex_name);
    sprintf(str, "cudaGetTextureReference on %s failed", tex_name);
    CU_RET_ERR(stat, str);
    cd = cudaCreateChannelDesc<T>();

    stat = cudaBindTexture(&offset, tex, d_ptr, &cd, size*sizeof(*d_ptr));
    sprintf(str, "cudaBindTexture on %s failed ", tex_name);
    CU_RET_ERR(stat, str);

    return offset;
}

/*! Binds a float texture tex_name to the GPU global memory pointed by d_ptr.
 *
 *  Explicit instantiation of cu_bind_texture for float type.
 */
template size_t cu_bind_texture<float>(const char *, const float *, int);

void cu_unbind_texture(const char *tex_name)
{
    cudaError_t             stat;
    const textureReference  *tex;
    char                    str[100];

    stat = cudaGetTextureReference(&tex, tex_name);
    sprintf(str, "cudaGetTextureReference on %s failed", tex_name);
    CU_RET_ERR(stat, str);
    stat = cudaUnbindTexture(tex);
    sprintf(str, "cudaUnbindTexture on %s failed ", tex_name);
    CU_RET_ERR(stat, str);
}


/**** Operation on buffered arrays (arrays with "over-allocation" in gmx wording) *****/

/*!
 * If the pointers to the size variables are NULL no resetting happens.
 */
void cu_free_buffered(void *d_ptr, int *n, int *nalloc)
{
    cudaError_t stat;

    if (d_ptr)
    {
        stat = cudaFree(d_ptr);
        CU_RET_ERR(stat, "cudaFree failed");
    }

    if (n)
    {
        *n = -1;
    }

    if (nalloc)
    {
        *nalloc = -1;
    }
}

/*!
 *  Reallocation of the memory pointed by d_ptr and copying of the data from 
 *  the location pointed by h_src host-side pointer is done. Allocation is 
 *  buffered and therefore freeing is only needed if the previously allocated 
 *  space is not enough.
 *  The H2D copy is launched in stream s and can be done synchronously or 
 *  assynchronously (the default is the latter).
 */
void cu_realloc_buffered(void **d_dest, void *h_src,
                         size_t type_size,
                         int *curr_size, int *curr_alloc_size,
                         int req_size,
                         cudaStream_t s,
                         gmx_bool doAsync = TRUE)
{
    cudaError_t stat;

    if (d_dest == NULL || req_size < 0)
    {
        return;
    }

    /* reallocate only if the data does not fit = allocation size is smaller 
       than the current requested size */
    if (req_size > *curr_alloc_size)
    {
        /* only free if the array has already been initialized */
        if (*curr_alloc_size >= 0)
        {
            cu_free_buffered(*d_dest, curr_size, curr_alloc_size);
        }

        *curr_alloc_size = 1.2 * req_size + 100;

        stat = cudaMalloc(d_dest, *curr_alloc_size * type_size);
        CU_RET_ERR(stat, "cudaMalloc failed in cu_free_buffered");
    }

    /* size could have changed without actual reallocation */
    *curr_size = req_size;

    /* upload to device */
    if (h_src)
    {
        if (doAsync)
        {
            cu_copy_H2D_async(*d_dest, h_src, *curr_size * type_size, s);
        }
        else
        {
            cu_copy_H2D(*d_dest, h_src,  *curr_size * type_size);
        }
    }
}
