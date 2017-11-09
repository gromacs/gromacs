/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <stdlib.h>

#include <cassert>
#include <cstdio>

#include <string>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

int ocl_copy_H2D(cl_mem d_dest, void* h_src,
                 size_t offset, size_t bytes,
                 GpuApiCallBehavior transferKind,
                 cl_command_queue command_queue,
                 cl_event *copy_event)
{
    cl_int gmx_unused cl_error;

    if (d_dest == NULL || h_src == NULL || bytes == 0)
    {
        return -1;
    }

    switch (transferKind)
    {
        case GpuApiCallBehavior::Async:
            cl_error = clEnqueueWriteBuffer(command_queue, d_dest, CL_FALSE, offset, bytes, h_src, 0, NULL, copy_event);
            assert(cl_error == CL_SUCCESS);
            // TODO: handle errors
            break;

        case GpuApiCallBehavior::Sync:
            cl_error = clEnqueueWriteBuffer(command_queue, d_dest, CL_TRUE, offset, bytes, h_src, 0, NULL, copy_event);
            assert(cl_error == CL_SUCCESS);
            // TODO: handle errors
            break;

        default:
            throw;
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
    return ocl_copy_H2D(d_dest, h_src, offset, bytes, GpuApiCallBehavior::Async, command_queue, copy_event);
}

/*! \brief Launches synchronous host to device memory copy.
 */
int ocl_copy_H2D_sync(cl_mem d_dest, void * h_src,
                      size_t offset, size_t bytes,
                      cl_command_queue command_queue)
{
    return ocl_copy_H2D(d_dest, h_src, offset, bytes, GpuApiCallBehavior::Sync, command_queue, NULL);
}

int ocl_copy_D2H(void * h_dest, cl_mem d_src,
                 size_t offset, size_t bytes,
                 GpuApiCallBehavior transferKind,
                 cl_command_queue command_queue,
                 cl_event *copy_event)
{
    cl_int gmx_unused cl_error;

    if (h_dest == NULL || d_src == NULL || bytes == 0)
    {
        return -1;
    }

    switch (transferKind)
    {
        case GpuApiCallBehavior::Async:
            cl_error = clEnqueueReadBuffer(command_queue, d_src, CL_FALSE, offset, bytes, h_dest, 0, NULL, copy_event);
            assert(cl_error == CL_SUCCESS);
            // TODO: handle errors
            break;

        case GpuApiCallBehavior::Sync:
            cl_error = clEnqueueReadBuffer(command_queue, d_src, CL_TRUE, offset, bytes, h_dest, 0, NULL, copy_event);
            assert(cl_error == CL_SUCCESS);
            // TODO: handle errors
            break;

        default:
            throw;
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
    return ocl_copy_D2H(h_dest, d_src, offset, bytes, GpuApiCallBehavior::Async, command_queue, copy_event);
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
    /* Need a temporary type whose size is 1 byte, so that the
     * implementation of snew_aligned can cope without issuing
     * warnings. */
    char **temporary = reinterpret_cast<char **>(h_ptr);

    /* 16-byte alignment is required by the neighbour-searching code,
     * because it uses four-wide SIMD for bounding-box calculation.
     * However, when we organize using page-locked memory for
     * device-host transfers, it will probably need to be aligned to a
     * 4kb page, like CUDA does. */
    snew_aligned(*temporary, nbytes, 16);
}

/*! \brief Frees memory allocated with ocl_pmalloc.
 *
 * \param[in]    h_ptr   Buffer allocated with ocl_pmalloc that needs to be freed.
 */
void ocl_pfree(void *h_ptr)
{

    if (h_ptr)
    {
        sfree_aligned(h_ptr);
    }
    return;
}

/*! \brief Convert error code to diagnostic string */
std::string ocl_get_error_string(cl_int error)
{
    switch (error)
    {
        // run-time and JIT compiler errors
        case 0: return "CL_SUCCESS";
        case -1: return "CL_DEVICE_NOT_FOUND";
        case -2: return "CL_DEVICE_NOT_AVAILABLE";
        case -3: return "CL_COMPILER_NOT_AVAILABLE";
        case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case -5: return "CL_OUT_OF_RESOURCES";
        case -6: return "CL_OUT_OF_HOST_MEMORY";
        case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
        case -8: return "CL_MEM_COPY_OVERLAP";
        case -9: return "CL_IMAGE_FORMAT_MISMATCH";
        case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case -11: return "CL_BUILD_PROGRAM_FAILURE";
        case -12: return "CL_MAP_FAILURE";
        case -13: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
        case -14: return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
        case -15: return "CL_COMPILE_PROGRAM_FAILURE";
        case -16: return "CL_LINKER_NOT_AVAILABLE";
        case -17: return "CL_LINK_PROGRAM_FAILURE";
        case -18: return "CL_DEVICE_PARTITION_FAILED";
        case -19: return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

        // compile-time errors
        case -30: return "CL_INVALID_VALUE";
        case -31: return "CL_INVALID_DEVICE_TYPE";
        case -32: return "CL_INVALID_PLATFORM";
        case -33: return "CL_INVALID_DEVICE";
        case -34: return "CL_INVALID_CONTEXT";
        case -35: return "CL_INVALID_QUEUE_PROPERTIES";
        case -36: return "CL_INVALID_COMMAND_QUEUE";
        case -37: return "CL_INVALID_HOST_PTR";
        case -38: return "CL_INVALID_MEM_OBJECT";
        case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        case -40: return "CL_INVALID_IMAGE_SIZE";
        case -41: return "CL_INVALID_SAMPLER";
        case -42: return "CL_INVALID_BINARY";
        case -43: return "CL_INVALID_BUILD_OPTIONS";
        case -44: return "CL_INVALID_PROGRAM";
        case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
        case -46: return "CL_INVALID_KERNEL_NAME";
        case -47: return "CL_INVALID_KERNEL_DEFINITION";
        case -48: return "CL_INVALID_KERNEL";
        case -49: return "CL_INVALID_ARG_INDEX";
        case -50: return "CL_INVALID_ARG_VALUE";
        case -51: return "CL_INVALID_ARG_SIZE";
        case -52: return "CL_INVALID_KERNEL_ARGS";
        case -53: return "CL_INVALID_WORK_DIMENSION";
        case -54: return "CL_INVALID_WORK_GROUP_SIZE";
        case -55: return "CL_INVALID_WORK_ITEM_SIZE";
        case -56: return "CL_INVALID_GLOBAL_OFFSET";
        case -57: return "CL_INVALID_EVENT_WAIT_LIST";
        case -58: return "CL_INVALID_EVENT";
        case -59: return "CL_INVALID_OPERATION";
        case -60: return "CL_INVALID_GL_OBJECT";
        case -61: return "CL_INVALID_BUFFER_SIZE";
        case -62: return "CL_INVALID_MIP_LEVEL";
        case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
        case -64: return "CL_INVALID_PROPERTY";
        case -65: return "CL_INVALID_IMAGE_DESCRIPTOR";
        case -66: return "CL_INVALID_COMPILER_OPTIONS";
        case -67: return "CL_INVALID_LINKER_OPTIONS";
        case -68: return "CL_INVALID_DEVICE_PARTITION_COUNT";

        // extension errors
        case -1000: return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
        case -1001: return "CL_PLATFORM_NOT_FOUND_KHR";
        case -1002: return "CL_INVALID_D3D10_DEVICE_KHR";
        case -1003: return "CL_INVALID_D3D10_RESOURCE_KHR";
        case -1004: return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
        case -1005: return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
        default:    return "Unknown OpenCL error: " +
                   std::to_string(static_cast<int32_t>(error));
    }
}
