/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
#ifndef GMX_GPU_UTILS_CUDAUTILS_CUH
#define GMX_GPU_UTILS_CUDAUTILS_CUH

#include <cstdio>

#include <array>
#include <string>
#include <type_traits>

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gputraits.cuh"
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
 * \param[in]  deviceError  The error to assert cudaSuccess on.
 *
 * \returns A description of the API error. Returns '(CUDA error #0 (cudaSuccess): no error)' in case deviceError is cudaSuccess.
 */
inline std::string getDeviceErrorString(const cudaError_t deviceError)
{
    return formatString("CUDA error #%d (%s): %s.",
                        deviceError,
                        cudaGetErrorName(deviceError),
                        cudaGetErrorString(deviceError));
}

/*! \brief Check if API returned an error and throw an exception with information on it.
 *
 * \param[in]  deviceError  The error to assert cudaSuccess on.
 * \param[in]  errorMessage  Undecorated error message.
 *
 *  \throws InternalError if deviceError is not a success.
 */
inline void checkDeviceError(const cudaError_t deviceError, const std::string& errorMessage)
{
    if (deviceError != cudaSuccess)
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
    cudaError_t deviceError = cudaGetLastError();
    if (deviceError == cudaSuccess)
    {
        return;
    }

    // If we would find an error in a release build, we do not know
    // what is appropriate to do about it, so assert only for debug
    // builds.
    const std::string fullErrorMessage =
            errorMessage + " An unhandled error from a previous CUDA operation was detected. "
            + gmx::getDeviceErrorString(deviceError);
    GMX_ASSERT(deviceError == cudaSuccess, fullErrorMessage.c_str());
    // TODO When we evolve a better logging framework, use that
    // for release-build error reporting.
    gmx_warning("%s", fullErrorMessage.c_str());
}

} // namespace
} // namespace gmx

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
#    define CU_RET_ERR(deviceError, msg)                                                            \
        do                                                                                          \
        {                                                                                           \
            if ((deviceError) != cudaSuccess)                                                       \
            {                                                                                       \
                gmx_fatal(FARGS, "%s\n", ((msg) + gmx::getDeviceErrorString(deviceError)).c_str()); \
            }                                                                                       \
        } while (0)

#else /* CHECK_CUDA_ERRORS */

#    define CU_RET_ERR(status, msg) \
        do                          \
        {                           \
        } while (0)

#endif /* CHECK_CUDA_ERRORS */

// TODO: the 2 functions below are pretty much a constructor/destructor of a simple
// GPU table object. There is also almost self-contained fetchFromParamLookupTable()
// in cuda_kernel_utils.cuh. They could all live in a separate class/struct file.

/*! \brief  Returns true if all tasks in \p s have completed.
 *
 *  \param[in] deviceStream CUDA stream to check.
 *
 *  \returns True if all tasks enqueued in the stream \p deviceStream (at the time of this call) have completed.
 */
static inline bool haveStreamTasksCompleted(const DeviceStream& deviceStream)
{
    cudaError_t stat = cudaStreamQuery(deviceStream.stream());

    if (stat == cudaErrorNotReady)
    {
        // work is still in progress in the stream
        return false;
    }

    GMX_ASSERT(stat != cudaErrorInvalidResourceHandle,
               ("Stream identifier not valid. " + gmx::getDeviceErrorString(stat)).c_str());

    // cudaSuccess and cudaErrorNotReady are the expected return values
    CU_RET_ERR(stat, "Unexpected cudaStreamQuery failure. ");

    GMX_ASSERT(stat == cudaSuccess,
               ("Values other than cudaSuccess should have been explicitly handled. "
                + gmx::getDeviceErrorString(stat))
                       .c_str());

    return true;
}

/* Kernel launch helpers */

/*! \brief
 * A function for setting up a single CUDA kernel argument.
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
 * Compile-time recursive function for setting up a single CUDA kernel argument.
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
 * A wrapper function for setting up all the CUDA kernel arguments.
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

/*! \brief Launches the CUDA kernel and handles the errors.
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
    auto stat = cudaLaunchKernel(reinterpret_cast<void*>(kernel),
                                 gridSize,
                                 blockSize,
                                 const_cast<void**>(kernelArgs.data()),
                                 config.sharedMemorySize,
                                 deviceStream.stream());
    GMX_RELEASE_ASSERT(stat == cudaSuccess,
                       ("GPU kernel (" + std::string(kernelName)
                        + ") failed to launch: " + gmx::getDeviceErrorString(stat))
                               .c_str());
}

#endif
