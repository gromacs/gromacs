/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
#ifndef GMX_GPU_UTILS_DEVICEBUFFER_OCL_H
#define GMX_GPU_UTILS_DEVICEBUFFER_OCL_H

/*! \libinternal \file
 *  \brief Implements the DeviceBuffer type and routines for OpenCL.
 *  Should only be included directly by the main DeviceBuffer file devicebuffer.h.
 *  TODO: the intent is for DeviceBuffer to become a class.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *
 *  \inlibraryapi
 */

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gpu_utils.h" //only for GpuApiCallBehavior
#include "gromacs/gpu_utils/gputraits_ocl.h"
#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

/*! \libinternal \brief
 * Allocates a device-side buffer.
 * It is currently a caller's responsibility to call it only on not-yet allocated buffers.
 *
 * \tparam        ValueType            Raw value type of the \p buffer.
 * \param[in,out] buffer               Pointer to the device-side buffer.
 * \param[in]     numValues            Number of values to accommodate.
 * \param[in]     deviceContext        The buffer's device context-to-be.
 */
template<typename ValueType>
void allocateDeviceBuffer(DeviceBuffer<ValueType>* buffer, size_t numValues, const DeviceContext& deviceContext)
{
    GMX_ASSERT(buffer, "needs a buffer pointer");
    void*  hostPtr = nullptr;
    cl_int clError;
    *buffer = clCreateBuffer(
            deviceContext.context(), CL_MEM_READ_WRITE, numValues * sizeof(ValueType), hostPtr, &clError);
    GMX_RELEASE_ASSERT(clError == CL_SUCCESS,
                       gmx::formatString("clCreateBuffer failure (OpenCL error %d: %s)",
                                         clError,
                                         ocl_get_error_string(clError).c_str())
                               .c_str());
}

/*! \brief
 * Frees a device-side buffer.
 * This does not reset separately stored size/capacity integers,
 * as this is planned to be a destructor of DeviceBuffer as a proper class,
 * and no calls on \p buffer should be made afterwards.
 *
 * \param[in] buffer  Pointer to the buffer to free.
 */
template<typename DeviceBuffer>
void freeDeviceBuffer(DeviceBuffer* buffer)
{
    GMX_ASSERT(buffer, "needs a buffer pointer");
    if (*buffer)
    {
        cl_int clError = clReleaseMemObject(*buffer);
        GMX_RELEASE_ASSERT(clError == CL_SUCCESS,
                           gmx::formatString("clReleaseMemObject failed (OpenCL error %d: %s)",
                                             clError,
                                             ocl_get_error_string(clError).c_str())
                                   .c_str());
    }
}

/*! \brief
 * Performs the host-to-device data copy, synchronous or asynchronously on request.
 *
 * Note that synchronous copy will not synchronize the stream in case of zero \p numValues
 * because of the early return.
 *
 * \tparam        ValueType            Raw value type of the \p buffer.
 * \param[in,out] buffer               Pointer to the device-side buffer
 * \param[in]     hostBuffer           Pointer to the raw host-side memory, also typed \p ValueType
 * \param[in]     startingOffset       Offset (in values) at the device-side buffer to copy into.
 * \param[in]     numValues            Number of values to copy.
 * \param[in]     deviceStream         GPU stream to perform asynchronous copy in.
 * \param[in]     transferKind         Copy type: synchronous or asynchronous.
 * \param[out]    timingEvent          A pointer to the H2D copy timing event to be filled in.
 *                                     If the pointer is not null, the event can further be used
 *                                     to queue a wait for this operation or to query profiling information.
 */
template<typename ValueType>
void copyToDeviceBuffer(DeviceBuffer<ValueType>* buffer,
                        const ValueType*         hostBuffer,
                        size_t                   startingOffset,
                        size_t                   numValues,
                        const DeviceStream&      deviceStream,
                        GpuApiCallBehavior       transferKind,
                        CommandEvent*            timingEvent)
{
    if (numValues == 0)
    {
        return; // such calls are actually made with empty domains
    }
    GMX_ASSERT(buffer, "needs a buffer pointer");
    GMX_ASSERT(hostBuffer, "needs a host buffer pointer");
    cl_int       clError;
    const size_t offset = startingOffset * sizeof(ValueType);
    const size_t bytes  = numValues * sizeof(ValueType);
    switch (transferKind)
    {
        case GpuApiCallBehavior::Async:
            clError = clEnqueueWriteBuffer(
                    deviceStream.stream(), *buffer, CL_FALSE, offset, bytes, hostBuffer, 0, nullptr, timingEvent);
            GMX_RELEASE_ASSERT(
                    clError == CL_SUCCESS,
                    gmx::formatString("Asynchronous H2D copy failed (OpenCL error %d: %s)",
                                      clError,
                                      ocl_get_error_string(clError).c_str())
                            .c_str());
            break;

        case GpuApiCallBehavior::Sync:
            clError = clEnqueueWriteBuffer(
                    deviceStream.stream(), *buffer, CL_TRUE, offset, bytes, hostBuffer, 0, nullptr, timingEvent);
            GMX_RELEASE_ASSERT(
                    clError == CL_SUCCESS,
                    gmx::formatString("Synchronous H2D copy failed (OpenCL error %d: %s)",
                                      clError,
                                      ocl_get_error_string(clError).c_str())
                            .c_str());
            break;

        default: throw;
    }
}

/*! \brief
 * Performs the device-to-host data copy, synchronous or asynchronously on request.
 *
 * Note that synchronous copy will not synchronize the stream in case of zero \p numValues
 * because of the early return.
 *
 * \tparam        ValueType            Raw value type of the \p buffer.
 * \param[in,out] hostBuffer           Pointer to the raw host-side memory, also typed \p ValueType
 * \param[in]     buffer               Pointer to the device-side buffer
 * \param[in]     startingOffset       Offset (in values) at the device-side buffer to copy from.
 * \param[in]     numValues            Number of values to copy.
 * \param[in]     deviceStream         GPU stream to perform asynchronous copy in.
 * \param[in]     transferKind         Copy type: synchronous or asynchronous.
 * \param[out]    timingEvent          A pointer to the H2D copy timing event to be filled in.
 *                                     If the pointer is not null, the event can further be used
 *                                     to queue a wait for this operation or to query profiling information.
 */
template<typename ValueType>
void copyFromDeviceBuffer(ValueType*               hostBuffer,
                          DeviceBuffer<ValueType>* buffer,
                          size_t                   startingOffset,
                          size_t                   numValues,
                          const DeviceStream&      deviceStream,
                          GpuApiCallBehavior       transferKind,
                          CommandEvent*            timingEvent)
{
    if (numValues == 0)
    {
        return;
    }
    GMX_ASSERT(buffer, "needs a buffer pointer");
    GMX_ASSERT(hostBuffer, "needs a host buffer pointer");
    cl_int       clError;
    const size_t offset = startingOffset * sizeof(ValueType);
    const size_t bytes  = numValues * sizeof(ValueType);
    switch (transferKind)
    {
        case GpuApiCallBehavior::Async:
            clError = clEnqueueReadBuffer(
                    deviceStream.stream(), *buffer, CL_FALSE, offset, bytes, hostBuffer, 0, nullptr, timingEvent);
            GMX_RELEASE_ASSERT(
                    clError == CL_SUCCESS,
                    gmx::formatString("Asynchronous D2H copy failed (OpenCL error %d: %s)",
                                      clError,
                                      ocl_get_error_string(clError).c_str())
                            .c_str());
            break;

        case GpuApiCallBehavior::Sync:
            clError = clEnqueueReadBuffer(
                    deviceStream.stream(), *buffer, CL_TRUE, offset, bytes, hostBuffer, 0, nullptr, timingEvent);
            GMX_RELEASE_ASSERT(
                    clError == CL_SUCCESS,
                    gmx::formatString("Synchronous D2H copy failed (OpenCL error %d: %s)",
                                      clError,
                                      ocl_get_error_string(clError).c_str())
                            .c_str());
            break;

        default: throw;
    }
}

/*! \brief
 * Performs the device-to-device data copy, synchronous or asynchronously on request.
 *
 * \tparam        ValueType                Raw value type of the \p buffer.
 */
template<typename ValueType>
void copyBetweenDeviceBuffers(DeviceBuffer<ValueType>* /* destinationDeviceBuffer */,
                              DeviceBuffer<ValueType>* /* sourceDeviceBuffer */,
                              size_t /* numValues */,
                              const DeviceStream& /* deviceStream */,
                              GpuApiCallBehavior /* transferKind */,
                              CommandEvent* /*timingEvent*/)
{
    // OpenCL-TODO
    gmx_fatal(FARGS, "D2D copy stub was called. Not yet implemented in OpenCL.");
}

/*! \brief
 * Clears the device buffer asynchronously.
 *
 * \tparam        ValueType       Raw value type of the \p buffer.
 * \param[in,out] buffer          Pointer to the device-side buffer
 * \param[in]     startingOffset  Offset (in values) at the device-side buffer to start clearing at.
 * \param[in]     numValues       Number of values to clear.
 * \param[in]     deviceStream    GPU stream.
 */
template<typename ValueType>
void clearDeviceBufferAsync(DeviceBuffer<ValueType>* buffer,
                            size_t                   startingOffset,
                            size_t                   numValues,
                            const DeviceStream&      deviceStream)
{
    if (numValues == 0)
    {
        return;
    }

    /* Apple OpenCL breaks if clEnqueueFillBuffer does not return an event.
     * Creating and releasing an event should not be necessary and potentially
     * causes minor overhead, so we don't do it on other platforms.
     * See #4852. */
#ifdef __APPLE__
    static constexpr bool sc_needAppleWorkaround = true;
#else
    static constexpr bool sc_needAppleWorkaround = false;
#endif

    GMX_ASSERT(buffer, "needs a buffer pointer");
    const size_t    offset        = startingOffset * sizeof(ValueType);
    const size_t    bytes         = numValues * sizeof(ValueType);
    const int       pattern       = 0;
    const cl_uint   numWaitEvents = 0;
    const cl_event* waitEvents    = nullptr;
    cl_event        retEvent;

    cl_int clError = clEnqueueFillBuffer(deviceStream.stream(),
                                         *buffer,
                                         &pattern,
                                         sizeof(pattern),
                                         offset,
                                         bytes,
                                         numWaitEvents,
                                         waitEvents,
                                         sc_needAppleWorkaround ? &retEvent : nullptr);
    GMX_RELEASE_ASSERT(clError == CL_SUCCESS,
                       gmx::formatString("Couldn't clear the device buffer (OpenCL error %d: %s)",
                                         clError,
                                         ocl_get_error_string(clError).c_str())
                               .c_str());
    if (sc_needAppleWorkaround)
    {
        clReleaseEvent(retEvent);
    }
}

CLANG_DIAGNOSTIC_IGNORE("-Wunused-template")

/*! \brief Check the validity of the device buffer.
 *
 * Checks if the buffer is not nullptr and if its allocation is big enough.
 *
 * \param[in] buffer        Device buffer to be checked.
 * \param[in] requiredSize  Number of elements that the buffer will have to accommodate.
 *
 * \returns Whether the device buffer can be set.
 */
template<typename T>
static bool checkDeviceBuffer(DeviceBuffer<T> buffer, int requiredSize)
{
    const size_t requiredSizeBytes = requiredSize * sizeof(T);
    size_t       sizeBytes;
    cl_int retval = clGetMemObjectInfo(buffer, CL_MEM_SIZE, sizeof(sizeBytes), &sizeBytes, nullptr);
    GMX_ASSERT(retval == CL_SUCCESS,
               gmx::formatString("clGetMemObjectInfo failed with error code #%d", retval).c_str());
    GMX_ASSERT(sizeBytes >= requiredSizeBytes,
               "Number of atoms in device buffer is smaller then required size.");
    return retval == CL_SUCCESS && sizeBytes >= requiredSizeBytes;
}

//! Device texture wrapper.
using DeviceTexture = void*;

/*! \brief Create a texture object for an array of type ValueType.
 *
 * Creates the device buffer and copies read-only data for an array of type ValueType.
 *
 * \todo Decide if using image2d is most efficient.
 *
 * \tparam      ValueType      Raw data type.
 *
 * \param[out]  deviceBuffer   Device buffer to store data in.
 * \param[in]   hostBuffer     Host buffer to get date from.
 * \param[in]   numValues      Number of elements in the buffer.
 * \param[in]   deviceContext  GPU device context.
 */
template<typename ValueType>
void initParamLookupTable(DeviceBuffer<ValueType>* deviceBuffer,
                          DeviceTexture* /* deviceTexture */,
                          const ValueType*     hostBuffer,
                          int                  numValues,
                          const DeviceContext& deviceContext)
{
    GMX_ASSERT(hostBuffer, "Host buffer pointer can not be null");
    const size_t bytes = numValues * sizeof(ValueType);
    cl_int       clError;
    *deviceBuffer = clCreateBuffer(deviceContext.context(),
                                   CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY | CL_MEM_COPY_HOST_PTR,
                                   bytes,
                                   const_cast<ValueType*>(hostBuffer),
                                   &clError);

    GMX_RELEASE_ASSERT(clError == CL_SUCCESS,
                       gmx::formatString("Constant memory allocation failed (OpenCL error %d: %s)",
                                         clError,
                                         ocl_get_error_string(clError).c_str())
                               .c_str());
}

/*! \brief Release the OpenCL device buffer.
 *
 * \tparam        ValueType     Raw data type.
 *
 * \param[in,out] deviceBuffer  Device buffer to store data in.
 */
template<typename ValueType>
void destroyParamLookupTable(DeviceBuffer<ValueType>* deviceBuffer, const DeviceTexture& /* deviceTexture*/)
{
    freeDeviceBuffer(deviceBuffer);
}

CLANG_DIAGNOSTIC_RESET

template<typename ValueType>
ValueType* asMpiPointer(DeviceBuffer<ValueType>& /*buffer*/)
{
    assert(false);
    return nullptr;
}

#endif
