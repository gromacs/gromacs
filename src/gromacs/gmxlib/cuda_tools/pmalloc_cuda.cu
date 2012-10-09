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

#include "cudautils.cuh"
#include "pmalloc_cuda.h"

/*! Allocates nbytes of page-locked memory. 
 *  This memory should always be freed using pfree (or with the page-locked 
 *  free functions provied by the CUDA library).
 */
void pmalloc(void **h_ptr, size_t nbytes)
{
    cudaError_t stat;
    char        strbuf[STRLEN];
    int         flag = cudaHostAllocDefault;

    if (nbytes == 0)
    {
        *h_ptr = NULL;
        return;
    }

    CU_CHECK_PREV_ERR();

    stat = cudaMallocHost(h_ptr, nbytes, flag);    
    sprintf(strbuf, "cudaMallocHost of size %d bytes failed", (int)nbytes);
    CU_RET_ERR(stat, strbuf);  
}

/*! Allocates nbytes of page-locked memory with write-combining. 
 *  This memory should always be freed using pfree (or with the page-locked 
 *  free functions provied by the CUDA library).
 */
void pmalloc_wc(void **h_ptr, size_t nbytes)
{
    cudaError_t stat;
    char        strbuf[STRLEN];
    int         flag = cudaHostAllocDefault || cudaHostAllocWriteCombined;

    if (nbytes == 0)
    {
        *h_ptr = NULL;
        return;
    }

    CU_CHECK_PREV_ERR();

    stat = cudaMallocHost(h_ptr, nbytes, flag);    
    sprintf(strbuf, "cudaMallocHost of size %d bytes failed", (int)nbytes);
    CU_RET_ERR(stat, strbuf);  
}

/*! Frees page locked memory allocated with pmalloc.
 *  This function can safely be called also with a pointer to a page-locked 
 *  memory allocated directly with CUDA API calls.
 */
void pfree(void *h_ptr) 
{
    cudaError_t stat; 

    if (h_ptr == NULL)
    {        
        return;
    }

    CU_CHECK_PREV_ERR();

    stat = cudaFreeHost(h_ptr);
    CU_RET_ERR(stat, "cudaFreeHost failed");
}
