/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 *  \brief Define utility routines for OpenCL
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 */
#include "gmxpre.h"

#include "oclutils.h"

#include <cassert>
#include <cstdio>
#include <stdlib.h>

#include "gromacs/utility/fatalerror.h"

/*! \brief Launches synchronous or asynchronous host to device memory copy.
 *
 *  If copy_event is not NULL, on return it will contain an event object
 *  identifying this particular host to device operation. The event can further
 *  be used to queue a wait for this operation or to query profiling information.
 */
static int ocl_copy_H2D_generic(cl_mem d_dest, void* h_src,
                                size_t offset, size_t bytes,
                                bool bAsync /* = false*/,
                                cl_command_queue command_queue,
                                cl_event *copy_event)
{
    cl_int gmx_unused cl_error;

    if (d_dest == NULL || h_src == NULL || bytes == 0)
    {
        return -1;
    }

    if (bAsync)
    {
        cl_error = clEnqueueWriteBuffer(command_queue, d_dest, CL_FALSE, offset, bytes, h_src, 0, NULL, copy_event);
        assert(cl_error == CL_SUCCESS);
        // TODO: handle errors
    }
    else
    {
        cl_error = clEnqueueWriteBuffer(command_queue, d_dest, CL_TRUE, offset, bytes, h_src, 0, NULL, copy_event);
        assert(cl_error == CL_SUCCESS);
        // TODO: handle errors
    }

    return 0;
}

/*! \brief Launches asynchronous host to device memory copy.
 *
 *  If copy_event is not NULL, on return it will contain an event object
 *  identifying this particular host to device operation. The event can further
 *  be used to queue a wait for this operation or to query profiling information.
 */
int ocl_copy_H2D_async(cl_mem d_dest, void * h_src,
                       size_t offset, size_t bytes,
                       cl_command_queue command_queue,
                       cl_event *copy_event)
{
    return ocl_copy_H2D_generic(d_dest, h_src, offset, bytes, true, command_queue, copy_event);
}

/*! \brief Launches synchronous host to device memory copy.
 */
int ocl_copy_H2D(cl_mem d_dest, void * h_src,
                 size_t offset, size_t bytes,
                 cl_command_queue command_queue)
{
    return ocl_copy_H2D_generic(d_dest, h_src, offset, bytes, false, command_queue, NULL);
}

/*! \brief Launches synchronous or asynchronous device to host memory copy.
 *
 *  If copy_event is not NULL, on return it will contain an event object
 *  identifying this particular device to host operation. The event can further
 *  be used to queue a wait for this operation or to query profiling information.
 */
int ocl_copy_D2H_generic(void * h_dest, cl_mem d_src,
                         size_t offset, size_t bytes,
                         bool bAsync,
                         cl_command_queue command_queue,
                         cl_event *copy_event)
{
    cl_int gmx_unused cl_error;

    if (h_dest == NULL || d_src == NULL || bytes == 0)
    {
        return -1;
    }

    if (bAsync)
    {
        cl_error = clEnqueueReadBuffer(command_queue, d_src, CL_FALSE, offset, bytes, h_dest, 0, NULL, copy_event);
        assert(cl_error == CL_SUCCESS);
        // TODO: handle errors
    }
    else
    {
        cl_error = clEnqueueReadBuffer(command_queue, d_src, CL_TRUE, offset, bytes, h_dest, 0, NULL, copy_event);
        assert(cl_error == CL_SUCCESS);
        // TODO: handle errors
    }

    return 0;
}

/*! \brief Launches asynchronous device to host memory copy.
 *
 *  If copy_event is not NULL, on return it will contain an event object
 *  identifying this particular host to device operation. The event can further
 *  be used to queue a wait for this operation or to query profiling information.
 */
int ocl_copy_D2H_async(void * h_dest, cl_mem d_src,
                       size_t offset, size_t bytes,
                       cl_command_queue command_queue,
                       cl_event *copy_event)
{
    return ocl_copy_D2H_generic(h_dest, d_src, offset, bytes, true, command_queue, copy_event);
}

/*! \brief \brief Allocates nbytes of host memory. Use ocl_free to free memory allocated with this function.
 *
 *  \todo
 *  This function should allocate page-locked memory to help reduce D2H and H2D
 *  transfer times, similar with pmalloc from pmalloc_cuda.cu.
 *
 * \param[in,out]    h_ptr   Pointer where to store the address of the newly allocated buffer.
 * \param[in]        nbytes  Size in bytes of the buffer to be allocated.
 */
void ocl_pmalloc(void **h_ptr, size_t nbytes)
{
    char        strbuf[1024];

#ifndef NDEBUG
    printf("Warning, pmalloc in OpenCL is doing a normal alloc instead of page-locked alloc\n");
#endif

    if (nbytes == 0)
    {
        *h_ptr = NULL;
        return;
    }

    *h_ptr = calloc(1, nbytes);

    if (!*h_ptr)
    {
        sprintf(strbuf, "calloc of size %d bytes failed", (int)nbytes);
        gmx_fatal(FARGS, "%s: %s\n", __FUNCTION__, strbuf);
    }
}

/*! \brief Frees memory allocated with ocl_pmalloc.
 *
 * \param[in]    h_ptr   Buffer allocated with ocl_pmalloc that needs to be freed.
 */
void ocl_pfree(void *h_ptr)
{

#ifndef NDEBUG
    printf("Warning, pfree in OpenCL is not deallocating page-locked memory\n");
#endif

    if (h_ptr)
    {
        free(h_ptr);
    }
    return;
}
