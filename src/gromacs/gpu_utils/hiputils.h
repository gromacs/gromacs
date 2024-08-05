/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
#ifndef GMX_GPU_UTILS_HIPUTILS_H
#define GMX_GPU_UTILS_HIPUTILS_H

#include <cstdio>

#include <array>
#include <string>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gputraits_hip.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace
{

/*! \brief Add the API information on the specific error to the error message.
 *
 * \param[in]  deviceError  The error to assert hipSuccess on.
 *
 * \returns A description of the API error. Returns '(HIP error #0 (hipSuccess): no error)' in case deviceError is hipSuccess.
 */
inline std::string getDeviceErrorString(const hipError_t deviceError)
{
    return formatString("HIP error #%d (%s): %s.",
                        deviceError,
                        hipGetErrorName(deviceError),
                        hipGetErrorString(deviceError));
}

/*! \brief Check if API returned an error and throw an exception with information on it.
 *
 * \param[in]  deviceError  The error to assert hipSuccess on.
 * \param[in]  errorMessage  Undecorated error message.
 *
 *  \throws InternalError if deviceError is not a success.
 */
inline void checkDeviceError(const hipError_t deviceError, const std::string& errorMessage)
{
    if (deviceError != hipSuccess)
    {
        GMX_THROW(gmx::InternalError(errorMessage + " " + getDeviceErrorString(deviceError)));
    }
}

/*! \brief Helper function to ensure no pending error silently
 * disrupts error handling.
 *
 * Asserts in a debug build if an unhandled error is present. Issues a
 * warning at run time otherwise.
 *
 * \param[in]  errorMessage  Undecorated error message.
 */
inline void ensureNoPendingDeviceError(const std::string& errorMessage)
{
    // Ensure there is no pending error that would otherwise affect
    // the behaviour of future error handling.
    hipError_t deviceError = hipGetLastError();
    if (deviceError == hipSuccess)
    {
        return;
    }

    // If we would find an error in a release build, we do not know
    // what is appropriate to do about it, so assert only for debug
    // builds.
    const std::string fullErrorMessage =
            errorMessage + " An unhandled error from a previous HIP operation was detected. "
            + gmx::getDeviceErrorString(deviceError);
    GMX_ASSERT(deviceError == hipSuccess, fullErrorMessage.c_str());
    // TODO When we evolve a better logging framework, use that
    // for release-build error reporting.
    gmx_warning("%s", fullErrorMessage.c_str());
}

} // namespace
} // namespace gmx

/*! \brief  Returns true if all tasks in \p deviceStream have completed.
 *
 *  \param[in] deviceStream HIP stream to check.
 *
 *  \returns True if all tasks enqueued in the stream \p deviceStream (at the time of this call) have completed.
 */
static inline bool haveStreamTasksCompleted(const DeviceStream& deviceStream)
{
    hipError_t stat = hipStreamQuery(deviceStream.stream());

    if (stat == hipErrorNotReady)
    {
        // work is still in progress in the stream
        return false;
    }

    GMX_ASSERT(stat != hipErrorInvalidHandle,
               ("Stream identifier not valid. " + gmx::getDeviceErrorString(stat)).c_str());

    // hipSuccess and hipErrorNotReady are the expected return values
    gmx::checkDeviceError(stat, "Unexpected hipStreamQuery failure. ");

    GMX_ASSERT(stat == hipSuccess,
               ("Values other than hipSuccess should have been explicitly handled. "
                + gmx::getDeviceErrorString(stat))
                       .c_str());

    return true;
}

/* Kernel launch helpers */

/*! \brief
 * A function for setting up a single HIP kernel argument.
 * This is the tail of the compile-time recursive function below.
 * It has to be seen by the compiler first.
 *
 * \tparam        totalArgsCount  Number of the kernel arguments
 * \tparam        KernelPtr       Kernel function handle type
 * \param[in]     argIndex        Index of the current argument
 */
template<size_t totalArgsCount, typename KernelPtr>
void prepareGpuKernelArgument(KernelPtr /*kernel*/,
                              std::array<void*, totalArgsCount>* /* kernelArgsPtr */,
                              size_t gmx_used_in_debug argIndex)
{
    GMX_ASSERT(argIndex == totalArgsCount, "Tail expansion");
}

/*! \brief
 * Compile-time recursive function for setting up a single HIP kernel argument.
 * This function copies a kernel argument pointer \p argPtr into \p kernelArgsPtr,
 * and calls itself on the next argument, eventually calling the tail function above.
 *
 * \tparam        CurrentArg      Type of the current argument
 * \tparam        RemainingArgs   Types of remaining arguments after the current one
 * \tparam        totalArgsCount  Number of the kernel arguments
 * \tparam        KernelPtr       Kernel function handle type
 * \param[in]     kernel          Kernel function handle
 * \param[in,out] kernelArgsPtr   Pointer to the argument array to be filled in
 * \param[in]     argIndex        Index of the current argument
 * \param[in]     argPtr          Pointer to the current argument
 * \param[in]     otherArgsPtrs   Pack of pointers to arguments remaining to process after the current one
 */
template<typename CurrentArg, typename... RemainingArgs, size_t totalArgsCount, typename KernelPtr>
void prepareGpuKernelArgument(KernelPtr                          kernel,
                              std::array<void*, totalArgsCount>* kernelArgsPtr,
                              size_t                             argIndex,
                              const CurrentArg*                  argPtr,
                              const RemainingArgs*... otherArgsPtrs)
{
    (*kernelArgsPtr)[argIndex] = const_cast<void*>(static_cast<const void*>(argPtr));
    prepareGpuKernelArgument(kernel, kernelArgsPtr, argIndex + 1, otherArgsPtrs...);
}

/*! \brief
 * A wrapper function for setting up all the HIP kernel arguments.
 * Calls the recursive functions above.
 *
 * \tparam    KernelPtr       Kernel function handle type
 * \tparam    Args            Types of all the kernel arguments
 * \param[in] kernel          Kernel function handle
 * \param[in] argsPtrs        Pointers to all the kernel arguments
 * \returns A prepared parameter pack to be used with launchGpuKernel() as the last argument.
 */
template<typename KernelPtr, typename... Args>
std::array<void*, sizeof...(Args)> prepareGpuKernelArguments(KernelPtr kernel,
                                                             const KernelLaunchConfig& /*config */,
                                                             const Args*... argsPtrs)
{
    std::array<void*, sizeof...(Args)> kernelArgs;
    prepareGpuKernelArgument(kernel, &kernelArgs, 0, argsPtrs...);
    return kernelArgs;
}

/*! \brief Launches the HIP kernel and handles the errors.
 *
 * \tparam    Args            Types of all the kernel arguments
 * \param[in] kernel          Kernel function handle
 * \param[in] config          Kernel configuration for launching
 * \param[in] deviceStream    GPU stream to launch kernel in
 * \param[in] kernelName      Human readable kernel description, for error handling only
 * \param[in] kernelArgs      Array of the pointers to the kernel arguments, prepared by
 *                            prepareGpuKernelArguments()
 * \throws gmx::InternalError on kernel launch failure
 */
template<typename... Args>
void launchGpuKernel(void (*kernel)(Args...),
                     const KernelLaunchConfig& config,
                     const DeviceStream&       deviceStream,
                     CommandEvent* /*timingEvent */,
                     const char*                               kernelName,
                     const std::array<void*, sizeof...(Args)>& kernelArgs)
{
    dim3 blockSize(config.blockSize[0], config.blockSize[1], config.blockSize[2]);
    dim3 gridSize(config.gridSize[0], config.gridSize[1], config.gridSize[2]);

    hipError_t stat = hipLaunchKernel(reinterpret_cast<void*>(kernel),
                                      gridSize,
                                      blockSize,
                                      const_cast<void**>(kernelArgs.data()),
                                      config.sharedMemorySize,
                                      deviceStream.stream());
    GMX_RELEASE_ASSERT(stat == hipSuccess,
                       ("GPU kernel (" + std::string(kernelName)
                        + ") failed to launch: " + gmx::getDeviceErrorString(stat))
                               .c_str());
}

#endif
