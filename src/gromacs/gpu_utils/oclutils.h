/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include <string>

#include "gromacs/gpu_utils/gmxopencl.h"
#include "gromacs/gpu_utils/gputraits_ocl.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

enum class GpuApiCallBehavior;

/*! \brief OpenCL vendor IDs */
typedef enum
{
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
    cl_platform_id ocl_platform_id; /**< Platform ID */
    cl_device_id   ocl_device_id;   /**< Device ID */
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
    ocl_gpu_id_t    ocl_gpu_id;          /**< device ID assigned at detection   */
    char            device_name[256];    /**< device name */
    char            device_version[256]; /**< device version */
    char            device_vendor[256];  /**< device vendor */
    int             compute_units;       /**< number of compute units */
    int             adress_bits;         /**< number of adress bits the device is capable of */
    int             stat;                /**< device status takes values of e_gpu_detect_res_t */
    ocl_vendor_id_t vendor_e;            /**< device vendor as defined by ocl_vendor_id_t */
    size_t maxWorkItemSizes[3]; /**< workgroup size limits (CL_DEVICE_MAX_WORK_ITEM_SIZES) */
    size_t maxWorkGroupSize;    /**< workgroup total size limit (CL_DEVICE_MAX_WORK_GROUP_SIZE) */
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

/*! \brief Launches synchronous or asynchronous device to host memory copy.
 *
 *  If copy_event is not NULL, on return it will contain an event object
 *  identifying this particular device to host operation. The event can further
 *  be used to queue a wait for this operation or to query profiling information.
 */
int ocl_copy_D2H(void*              h_dest,
                 cl_mem             d_src,
                 size_t             offset,
                 size_t             bytes,
                 GpuApiCallBehavior transferKind,
                 cl_command_queue   command_queue,
                 cl_event*          copy_event);


/*! \brief Launches asynchronous device to host memory copy. */
int ocl_copy_D2H_async(void*            h_dest,
                       cl_mem           d_src,
                       size_t           offset,
                       size_t           bytes,
                       cl_command_queue command_queue,
                       cl_event*        copy_event);

/*! \brief Launches synchronous or asynchronous host to device memory copy.
 *
 *  If copy_event is not NULL, on return it will contain an event object
 *  identifying this particular host to device operation. The event can further
 *  be used to queue a wait for this operation or to query profiling information.
 */
int ocl_copy_H2D(cl_mem             d_dest,
                 const void*        h_src,
                 size_t             offset,
                 size_t             bytes,
                 GpuApiCallBehavior transferKind,
                 cl_command_queue   command_queue,
                 cl_event*          copy_event);

/*! \brief Launches asynchronous host to device memory copy. */
int ocl_copy_H2D_async(cl_mem           d_dest,
                       const void*      h_src,
                       size_t           offset,
                       size_t           bytes,
                       cl_command_queue command_queue,
                       cl_event*        copy_event);

/*! \brief Launches synchronous host to device memory copy. */
int ocl_copy_H2D_sync(cl_mem d_dest, const void* h_src, size_t offset, size_t bytes, cl_command_queue command_queue);

/*! \brief Allocate host memory in malloc style */
void pmalloc(void** h_ptr, size_t nbytes);

/*! \brief Free host memory in malloc style */
void pfree(void* h_ptr);

/*! \brief Convert error code to diagnostic string */
std::string ocl_get_error_string(cl_int error);

/*! \brief Calls clFinish() in the stream \p s.
 *
 * \param[in] s stream to synchronize with
 */
static inline void gpuStreamSynchronize(cl_command_queue s)
{
    cl_int cl_error = clFinish(s);
    GMX_RELEASE_ASSERT(CL_SUCCESS == cl_error,
                       ("Error caught during clFinish:" + ocl_get_error_string(cl_error)).c_str());
}

//! A debug checker to track cl_events being released correctly
inline void ensureReferenceCount(const cl_event& event, unsigned int refCount)
{
#ifndef NDEBUG
    cl_int clError = clGetEventInfo(event, CL_EVENT_REFERENCE_COUNT, sizeof(refCount), &refCount, nullptr);
    GMX_ASSERT(CL_SUCCESS == clError, ocl_get_error_string(clError).c_str());
    GMX_ASSERT(refCount == refCount, "Unexpected reference count");
#else
    GMX_UNUSED_VALUE(event);
    GMX_UNUSED_VALUE(refCount);
#endif
}

/*! \brief Pretend to synchronize an OpenCL stream (dummy implementation).
 *
 * \param[in] s queue to check
 *
 *  \returns     True if all tasks enqueued in the stream \p s (at the time of this call) have completed.
 */
static inline bool haveStreamTasksCompleted(cl_command_queue gmx_unused s)
{
    GMX_RELEASE_ASSERT(false, "haveStreamTasksCompleted is not implemented for OpenCL");
    return false;
}

/* Kernel launch helpers */

/*! \brief
 * A function for setting up a single OpenCL kernel argument.
 * This is the tail of the compile-time recursive function below.
 * It has to be seen by the compiler first.
 * As NB kernels might be using dynamic local memory as the last argument,
 * this function also manages that, using sharedMemorySize from \p config.
 *
 * \param[in]     kernel          Kernel function handle
 * \param[in]     config          Kernel configuration for launching
 * \param[in]     argIndex        Index of the current argument
 */
void inline prepareGpuKernelArgument(cl_kernel kernel, const KernelLaunchConfig& config, size_t argIndex)
{
    if (config.sharedMemorySize > 0)
    {
        cl_int gmx_used_in_debug clError =
                clSetKernelArg(kernel, argIndex, config.sharedMemorySize, nullptr);
        GMX_ASSERT(CL_SUCCESS == clError, ocl_get_error_string(clError).c_str());
    }
}

/*! \brief
 * Compile-time recursive function for setting up a single OpenCL kernel argument.
 * This function uses one kernel argument pointer \p argPtr to call clSetKernelArg(),
 * and calls itself on the next argument, eventually calling the tail function above.
 *
 * \tparam        CurrentArg      Type of the current argument
 * \tparam        RemainingArgs   Types of remaining arguments after the current one
 * \param[in]     kernel          Kernel function handle
 * \param[in]     config          Kernel configuration for launching
 * \param[in]     argIndex        Index of the current argument
 * \param[in]     argPtr          Pointer to the current argument
 * \param[in]     otherArgsPtrs   Pack of pointers to arguments remaining to process after the current one
 */
template<typename CurrentArg, typename... RemainingArgs>
void prepareGpuKernelArgument(cl_kernel                 kernel,
                              const KernelLaunchConfig& config,
                              size_t                    argIndex,
                              const CurrentArg*         argPtr,
                              const RemainingArgs*... otherArgsPtrs)
{
    cl_int gmx_used_in_debug clError = clSetKernelArg(kernel, argIndex, sizeof(CurrentArg), argPtr);
    GMX_ASSERT(CL_SUCCESS == clError, ocl_get_error_string(clError).c_str());

    // Assert on types not allowed to be passed to a kernel
    // (as per section 6.9 of the OpenCL spec).
    static_assert(!std::is_same<CurrentArg, bool>::value && !std::is_same<CurrentArg, size_t>::value
                          && !std::is_same<CurrentArg, ptrdiff_t>::value
                          && !std::is_same<CurrentArg, intptr_t>::value
                          && !std::is_same<CurrentArg, uintptr_t>::value,
                  "Invalid type passed to OpenCL kernel functions (see OpenCL spec section 6.9).");

    prepareGpuKernelArgument(kernel, config, argIndex + 1, otherArgsPtrs...);
}

/*! \brief
 * A wrapper function for setting up all the OpenCL kernel arguments.
 * Calls the recursive functions above.
 *
 * \tparam    Args            Types of all the kernel arguments
 * \param[in] kernel          Kernel function handle
 * \param[in] config          Kernel configuration for launching
 * \param[in] argsPtrs        Pointers to all the kernel arguments
 * \returns A handle for the prepared parameter pack to be used with launchGpuKernel() as the last argument
 * - currently always nullptr for OpenCL, as it manages kernel/arguments association by itself.
 */
template<typename... Args>
void* prepareGpuKernelArguments(cl_kernel kernel, const KernelLaunchConfig& config, const Args*... argsPtrs)
{
    prepareGpuKernelArgument(kernel, config, 0, argsPtrs...);
    return nullptr;
}

/*! \brief Launches the OpenCL kernel and handles the errors.
 *
 * \param[in] kernel          Kernel function handle
 * \param[in] config          Kernel configuration for launching
 * \param[in] timingEvent     Timing event, fetched from GpuRegionTimer
 * \param[in] kernelName      Human readable kernel description, for error handling only
 * \throws gmx::InternalError on kernel launch failure
 */
inline void launchGpuKernel(cl_kernel                 kernel,
                            const KernelLaunchConfig& config,
                            CommandEvent*             timingEvent,
                            const char*               kernelName,
                            const void* /*kernelArgs*/)
{
    const int       workDimensions   = 3;
    const size_t*   globalWorkOffset = nullptr;
    const size_t    waitListSize     = 0;
    const cl_event* waitList         = nullptr;
    size_t          globalWorkSize[3];
    for (int i = 0; i < workDimensions; i++)
    {
        globalWorkSize[i] = config.gridSize[i] * config.blockSize[i];
    }
    cl_int clError = clEnqueueNDRangeKernel(config.stream, kernel, workDimensions, globalWorkOffset,
                                            globalWorkSize, config.blockSize, waitListSize,
                                            waitList, timingEvent);
    if (CL_SUCCESS != clError)
    {
        const std::string errorMessage = "GPU kernel (" + std::string(kernelName)
                                         + ") failed to launch: " + ocl_get_error_string(clError);
        GMX_THROW(gmx::InternalError(errorMessage));
    }
}

#endif
