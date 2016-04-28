/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
/*! \libinternal \file
 *  \brief Declare utility routines for OpenCL
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_OCLUTILS_H
#define GMX_GPU_UTILS_OCLUTILS_H

/*! \brief Declare to OpenCL SDKs that we intend to use OpenCL API
   features that were deprecated in 2.0, so that they don't warn about
   it. */
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#ifdef __APPLE__
#    include <OpenCL/opencl.h>
#else
#    include <CL/opencl.h>
#endif

#include <string>

/*! \brief OpenCL vendor IDs */
typedef enum {
    OCL_VENDOR_NVIDIA = 0,
    OCL_VENDOR_AMD,
    OCL_VENDOR_INTEL,
    OCL_VENDOR_UNKNOWN
} ocl_vendor_id_t;

/*! \internal
 * \brief OpenCL GPU device identificator
 *
 * An OpenCL device is identified by its ID.
 * The platform ID is also included for caching reasons.
 */
typedef struct
{
    cl_platform_id      ocl_platform_id; /**< Platform ID */
    cl_device_id        ocl_device_id;   /**< Device ID */
} ocl_gpu_id_t;

/*! \internal
 * \brief OpenCL device information.
 *
 * The OpenCL device information is queried and set at detection and contains
 * both information about the device/hardware returned by the runtime as well
 * as additional data like support status.
 */
struct gmx_device_info_t
{
    ocl_gpu_id_t        ocl_gpu_id;          /**< device ID assigned at detection   */
    char                device_name[256];    /**< device name */
    char                device_version[256]; /**< device version */
    char                device_vendor[256];  /**< device vendor */
    int                 compute_units;       /**< number of compute units */
    int                 adress_bits;         /**< number of adress bits the device is capable of */
    int                 stat;                /**< device status takes values of e_gpu_detect_res_t */
    ocl_vendor_id_t     vendor_e;            /**< device vendor as defined by ocl_vendor_id_t */
};

/*! \internal
 * \brief OpenCL GPU runtime data
 *
 * The device runtime data is meant to hold objects associated with a GROMACS rank's
 * (thread or process) use of a single device (multiple devices per rank is not
 * implemented). These objects should be constructed at ther point where a device
 * dets assigned to a rank and released at when this assignment is no longer valid
 * (i.e. at cleanup in the current implementation).
 *
 */
struct gmx_device_runtime_data_t
{
    cl_context context; /**< OpenCL context */
    cl_program program; /**< OpenCL program */
};

#if !defined(NDEBUG)
/* Debugger callable function that prints the name of a kernel function pointer */
cl_int dbg_ocl_kernel_name(const cl_kernel kernel);
cl_int dbg_ocl_kernel_name_address(void* kernel);
#endif


/*! \brief Launches asynchronous host to device memory copy. */
int ocl_copy_H2D_async(cl_mem d_dest, void * h_src,
                       size_t offset, size_t bytes,
                       cl_command_queue command_queue,
                       cl_event *copy_event);

/*! \brief Launches asynchronous device to host memory copy. */
int ocl_copy_D2H_async(void * h_dest, cl_mem d_src,
                       size_t offset, size_t bytes,
                       cl_command_queue command_queue,
                       cl_event *copy_event);

/*! \brief Launches synchronous host to device memory copy. */
int ocl_copy_H2D(cl_mem d_dest, void * h_src,
                 size_t offset, size_t bytes,
                 cl_command_queue command_queue);

/*! \brief Allocate host memory in malloc style */
void ocl_pmalloc(void **h_ptr, size_t nbytes);

/*! \brief Free host memory in malloc style */
void ocl_pfree(void *h_ptr);

/*! \brief Convert error code to diagnostic string */
std::string ocl_get_error_string(cl_int error);

#endif
