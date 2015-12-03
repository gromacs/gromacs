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
/*! \internal \file
 *  \brief Define functions for host-side memory handling when using CUDA devices.
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 */

#include "gmxpre.h"

#include "pmalloc_cuda.h"

#include <stdlib.h>

#include "gromacs/gmxlib/cuda_tools/cudautils.cuh"
#include "gromacs/utility/cstringutil.h"

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
    int         flag = cudaHostAllocDefault | cudaHostAllocWriteCombined;

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
