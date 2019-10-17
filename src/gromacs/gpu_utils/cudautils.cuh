/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_GPU_UTILS_CUDAUTILS_CUH
#define GMX_GPU_UTILS_CUDAUTILS_CUH

#include <stdio.h>

#include <array>
#include <string>

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

/*! \brief Helper function to ensure no pending error silently
 * disrupts error handling.
 *
 * Asserts in a debug build if an unhandled error is present. Issues a
 * warning at run time otherwise.
 *
 * \todo This is similar to CU_CHECK_PREV_ERR, which should be
 * consolidated.
 */
static inline void ensureNoPendingCudaError(const char* errorMessage)
{
    // Ensure there is no pending error that would otherwise affect
    // the behaviour of future error handling.
    cudaError_t stat = cudaGetLastError();
    if (stat == cudaSuccess)
    {
        return;
    }

    // If we would find an error in a release build, we do not know
    // what is appropriate to do about it, so assert only for debug
    // builds.
    auto fullMessage = formatString(
            "%s An unhandled error from a previous CUDA operation was detected. %s: %s",
            errorMessage, cudaGetErrorName(stat), cudaGetErrorString(stat));
    GMX_ASSERT(stat == cudaSuccess, fullMessage.c_str());
    // TODO When we evolve a better logging framework, use that
    // for release-build error reporting.
    gmx_warning("%s", fullMessage.c_str());
}

} // namespace
} // namespace gmx

enum class GpuApiCallBehavior;

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
#    define CU_RET_ERR(status, msg)                                            \
        do                                                                     \
        {                                                                      \
            if (status != cudaSuccess)                                         \
            {                                                                  \
                gmx_fatal(FARGS, "%s: %s\n", msg, cudaGetErrorString(status)); \
            }                                                                  \
        } while (0)

/*! Check for any previously occurred uncaught CUDA error. */
#    define CU_CHECK_PREV_ERR()                                                           \
        do                                                                                \
        {                                                                                 \
            cudaError_t _CU_CHECK_PREV_ERR_status = cudaGetLastError();                   \
            if (_CU_CHECK_PREV_ERR_status != cudaSuccess)                                 \
            {                                                                             \
                gmx_warning(                                                              \
                        "Just caught a previously occurred CUDA error (%s), will try to " \
                        "continue.",                                                      \
                        cudaGetErrorString(_CU_CHECK_PREV_ERR_status));                   \
            }                                                                             \
        } while (0)

#else /* CHECK_CUDA_ERRORS */

#    define CU_RET_ERR(status, msg) \
        do                          \
        {                           \
        } while (0)
#    define CU_CHECK_PREV_ERR() \
        do                      \
        {                       \
        } while (0)

#endif /* CHECK_CUDA_ERRORS */

/*! \brief CUDA device information.
 *
 * The CUDA device information is queried and set at detection and contains
 * both information about the device/hardware returned by the runtime as well
 * as additional data like support status.
 */
struct gmx_device_info_t
{
    int            id;   /* id of the CUDA device */
    cudaDeviceProp prop; /* CUDA device properties */
    int            stat; /* result of the device check */
};

/*! Launches synchronous or asynchronous device to host memory copy.
 *
 *  The copy is launched in stream s or if not specified, in stream 0.
 */
int cu_copy_D2H(void* h_dest, void* d_src, size_t bytes, GpuApiCallBehavior transferKind, cudaStream_t /*s = nullptr*/);

/*! Launches synchronous host to device memory copy in stream 0. */
int cu_copy_D2H_sync(void* /*h_dest*/, void* /*d_src*/, size_t /*bytes*/);

/*! Launches asynchronous host to device memory copy in stream s. */
int cu_copy_D2H_async(void* /*h_dest*/, void* /*d_src*/, size_t /*bytes*/, cudaStream_t /*s = nullptr*/);

/*! Launches synchronous or asynchronous host to device memory copy.
 *
 *  The copy is launched in stream s or if not specified, in stream 0.
 */
int cu_copy_H2D(void* d_dest, const void* h_src, size_t bytes, GpuApiCallBehavior transferKind, cudaStream_t /*s = nullptr*/);

/*! Launches synchronous host to device memory copy. */
int cu_copy_H2D_sync(void* /*d_dest*/, const void* /*h_src*/, size_t /*bytes*/);

/*! Launches asynchronous host to device memory copy in stream s. */
int cu_copy_H2D_async(void* /*d_dest*/, const void* /*h_src*/, size_t /*bytes*/, cudaStream_t /*s = nullptr*/);

// TODO: the 2 functions below are pretty much a constructor/destructor of a simple
// GPU table object. There is also almost self-contained fetchFromParamLookupTable()
// in cuda_kernel_utils.cuh. They could all live in a separate class/struct file.

/*! \brief Initialize parameter lookup table.
 *
 * Initializes device memory, copies data from host and binds
 * a texture to allocated device memory to be used for parameter lookup.
 *
 * \tparam[in] T         Raw data type
 * \param[out] d_ptr     device pointer to the memory to be allocated
 * \param[out] texObj    texture object to be initialized
 * \param[in]  h_ptr     pointer to the host memory to be uploaded to the device
 * \param[in]  numElem   number of elements in the h_ptr
 */
template<typename T>
void initParamLookupTable(T*& d_ptr, cudaTextureObject_t& texObj, const T* h_ptr, int numElem);

// Add extern declarations so each translation unit understands that
// there will be a definition provided.
extern template void initParamLookupTable<int>(int*&, cudaTextureObject_t&, const int*, int);
extern template void initParamLookupTable<float>(float*&, cudaTextureObject_t&, const float*, int);

/*! \brief Destroy parameter lookup table.
 *
 * Unbinds texture object, deallocates device memory.
 *
 * \tparam[in] T         Raw data type
 * \param[in]  d_ptr     Device pointer to the memory to be deallocated
 * \param[in]  texObj    Texture object to be deinitialized
 */
template<typename T>
void destroyParamLookupTable(T* d_ptr, cudaTextureObject_t texObj);

// Add extern declarations so each translation unit understands that
// there will be a definition provided.
extern template void destroyParamLookupTable<int>(int*, cudaTextureObject_t);
extern template void destroyParamLookupTable<float>(float*, cudaTextureObject_t);

/*! \brief Add a triplets stored in a float3 to an rvec variable.
 *
 * \param[out]  a Rvec to increment
 * \param[in]   b Float triplet to increment with.
 */
static inline void rvec_inc(rvec a, const float3 b)
{
    rvec tmp = { b.x, b.y, b.z };
    rvec_inc(a, tmp);
}

/*! \brief Wait for all taks in stream \p s to complete.
 *
 * \param[in] s stream to synchronize with
 */
static inline void gpuStreamSynchronize(cudaStream_t s)
{
    cudaError_t stat = cudaStreamSynchronize(s);
    CU_RET_ERR(stat, "cudaStreamSynchronize failed");
}

/*! \brief  Returns true if all tasks in \p s have completed.
 *
 * \param[in] s stream to check
 *
 *  \returns     True if all tasks enqueued in the stream \p s (at the time of this call) have completed.
 */
static inline bool haveStreamTasksCompleted(cudaStream_t s)
{
    cudaError_t stat = cudaStreamQuery(s);

    if (stat == cudaErrorNotReady)
    {
        // work is still in progress in the stream
        return false;
    }

    GMX_ASSERT(stat != cudaErrorInvalidResourceHandle, "Stream idnetifier not valid");

    // cudaSuccess and cudaErrorNotReady are the expected return values
    CU_RET_ERR(stat, "Unexpected cudaStreamQuery failure");

    GMX_ASSERT(stat == cudaSuccess,
               "Values other than cudaSuccess should have been explicitly handled");

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
    (*kernelArgsPtr)[argIndex] = (void*)argPtr;
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
 * \param[in] kernelName      Human readable kernel description, for error handling only
 * \param[in] kernelArgs      Array of the pointers to the kernel arguments, prepared by
 * prepareGpuKernelArguments() \throws gmx::InternalError on kernel launch failure
 */
template<typename... Args>
void launchGpuKernel(void (*kernel)(Args...),
                     const KernelLaunchConfig& config,
                     CommandEvent* /*timingEvent */,
                     const char*                               kernelName,
                     const std::array<void*, sizeof...(Args)>& kernelArgs)
{
    dim3 blockSize(config.blockSize[0], config.blockSize[1], config.blockSize[2]);
    dim3 gridSize(config.gridSize[0], config.gridSize[1], config.gridSize[2]);
    cudaLaunchKernel((void*)kernel, gridSize, blockSize, const_cast<void**>(kernelArgs.data()),
                     config.sharedMemorySize, config.stream);

    cudaError_t status = cudaGetLastError();
    if (cudaSuccess != status)
    {
        const std::string errorMessage =
                "GPU kernel (" + std::string(kernelName)
                + ") failed to launch: " + std::string(cudaGetErrorString(status));
        GMX_THROW(gmx::InternalError(errorMessage));
    }
}

#endif
