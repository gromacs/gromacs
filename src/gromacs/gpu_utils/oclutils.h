/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gmxopencl.h"
#include "gromacs/gpu_utils/gputraits_ocl.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

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
    //! OpenCL program
    cl_program program;
};

/*! \brief Convert error code to diagnostic string */
std::string ocl_get_error_string(cl_int error);

/*! \brief Pretend to synchronize an OpenCL stream (dummy implementation).
 *
 *  \returns  Not implemented in OpenCL.
 */
static inline bool haveStreamTasksCompleted(const DeviceStream& /* deviceStream */)
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
    static_assert(
            !std::is_same_v<CurrentArg,
                            bool> && !std::is_same_v<CurrentArg, size_t> && !std::is_same_v<CurrentArg, ptrdiff_t> && !std::is_same_v<CurrentArg, intptr_t> && !std::is_same_v<CurrentArg, uintptr_t>,
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
 * \param[in] deviceStream    GPU stream to launch kernel in
 * \param[in] timingEvent     Timing event, fetched from GpuRegionTimer
 * \param[in] kernelName      Human readable kernel description, for error handling only
 * \throws gmx::InternalError on kernel launch failure
 */
inline void launchGpuKernel(cl_kernel                 kernel,
                            const KernelLaunchConfig& config,
                            const DeviceStream&       deviceStream,
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
    cl_int clError = clEnqueueNDRangeKernel(deviceStream.stream(),
                                            kernel,
                                            workDimensions,
                                            globalWorkOffset,
                                            globalWorkSize,
                                            config.blockSize,
                                            waitListSize,
                                            waitList,
                                            timingEvent);
    if (CL_SUCCESS != clError)
    {
        const std::string errorMessage = "GPU kernel (" + std::string(kernelName)
                                         + ") failed to launch: " + ocl_get_error_string(clError);
        GMX_THROW(gmx::InternalError(errorMessage));
    }
}

#endif
