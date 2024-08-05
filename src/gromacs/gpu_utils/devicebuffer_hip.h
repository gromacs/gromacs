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
#ifndef GMX_GPU_UTILS_DEVICEBUFFER_HIP_H
#define GMX_GPU_UTILS_DEVICEBUFFER_HIP_H

/*! \libinternal \file
 *  \brief Implements the DeviceBuffer type and routines for HIP.
 *  Should only be included directly by the main DeviceBuffer file devicebuffer.h.
 *  TODO: the intent is for DeviceBuffer to become a class.
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *  \author Julio Maia <julio.maia@amd.com>
 *
 *  \inlibraryapi
 */

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gpu_utils.h" //only for GpuApiCallBehavior
#include "gromacs/gpu_utils/gputraits_hip.h"
#include "gromacs/gpu_utils/hiputils.h"
#include "gromacs/utility/gmxassert.h"

/*! \brief
 * Allocates a device-side buffer.
 * It is currently a caller's responsibility to call it only on not-yet allocated buffers.
 *
 * \tparam        ValueType            Raw value type of the \p buffer.
 * \param[in,out] buffer               Pointer to the device-side buffer.
 * \param[in]     numValues            Number of values to accommodate.
 */
template<typename ValueType>
void allocateDeviceBuffer(DeviceBuffer<ValueType>* buffer, size_t numValues, const DeviceContext& /* deviceContext */)
{
    GMX_ASSERT(buffer, "needs a buffer pointer");
    hipError_t stat = hipMalloc(buffer, numValues * sizeof(ValueType));
    GMX_RELEASE_ASSERT(
            stat == hipSuccess,
            ("Allocation of the device buffer failed. " + gmx::getDeviceErrorString(stat)).c_str());
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
        hipError_t stat = hipFree(*buffer);
        GMX_RELEASE_ASSERT(
                stat == hipSuccess,
                ("Freeing of the device buffer failed. " + gmx::getDeviceErrorString(stat)).c_str());
    }
}

/*! \brief
 * Performs the host-to-device data copy, synchronous or asynchronously on request.
 *
 * \tparam        ValueType            Raw value type of the \p buffer.
 * \param[in,out] buffer               Pointer to the device-side buffer
 * \param[in]     hostBuffer           Pointer to the raw host-side memory, also typed \p ValueType
 * \param[in]     startingOffset       Offset (in values) at the device-side buffer to copy into.
 * \param[in]     numValues            Number of values to copy.
 * \param[in]     deviceStream         GPU stream to perform asynchronous copy in.
 * \param[in]     transferKind         Copy type: synchronous or asynchronous.
 */
template<typename ValueType>
void copyToDeviceBuffer(DeviceBuffer<ValueType>* buffer,
                        const ValueType*         hostBuffer,
                        size_t                   startingOffset,
                        size_t                   numValues,
                        const DeviceStream&      deviceStream,
                        GpuApiCallBehavior       transferKind,
                        CommandEvent* /*timingEvent*/)
{
    if (numValues == 0)
    {
        return;
    }
    GMX_ASSERT(buffer, "needs a buffer pointer");
    GMX_ASSERT(hostBuffer, "needs a host buffer pointer");
    hipError_t   stat;
    const size_t bytes = numValues * sizeof(ValueType);
    switch (transferKind)
    {
        case GpuApiCallBehavior::Async:
            GMX_ASSERT(isHostMemoryPinned(hostBuffer), "Source host buffer was not pinned for HIP");
            stat = hipMemcpyAsync(*reinterpret_cast<ValueType**>(buffer) + startingOffset,
                                  hostBuffer,
                                  bytes,
                                  hipMemcpyHostToDevice,
                                  deviceStream.stream());
            GMX_RELEASE_ASSERT(
                    stat == hipSuccess,
                    ("Asynchronous H2D copy failed. " + gmx::getDeviceErrorString(stat)).c_str());
            break;

        case GpuApiCallBehavior::Sync:
            stat = hipMemcpy(*reinterpret_cast<ValueType**>(buffer) + startingOffset,
                             hostBuffer,
                             bytes,
                             hipMemcpyHostToDevice);
            GMX_RELEASE_ASSERT(
                    stat == hipSuccess,
                    ("Synchronous H2D copy failed. " + gmx::getDeviceErrorString(stat)).c_str());
            break;

        default: throw;
    }
}

/*! \brief
 * Performs the device-to-host data copy, synchronous or asynchronously on request.
 *
 * \tparam        ValueType            Raw value type of the \p buffer.
 * \param[in,out] hostBuffer           Pointer to the raw host-side memory, also typed \p ValueType
 * \param[in]     buffer               Pointer to the device-side buffer
 * \param[in]     startingOffset       Offset (in values) at the device-side buffer to copy from.
 * \param[in]     numValues            Number of values to copy.
 * \param[in]     deviceStream         GPU stream to perform asynchronous copy in.
 * \param[in]     transferKind         Copy type: synchronous or asynchronous.
 */
template<typename ValueType>
void copyFromDeviceBuffer(ValueType*               hostBuffer,
                          DeviceBuffer<ValueType>* buffer,
                          size_t                   startingOffset,
                          size_t                   numValues,
                          const DeviceStream&      deviceStream,
                          GpuApiCallBehavior       transferKind,
                          CommandEvent* /*timingEvent*/)
{
    if (numValues == 0)
    {
        return;
    }
    GMX_ASSERT(buffer, "needs a buffer pointer");
    GMX_ASSERT(hostBuffer, "needs a host buffer pointer");

    hipError_t   stat;
    const size_t bytes = numValues * sizeof(ValueType);
    switch (transferKind)
    {
        case GpuApiCallBehavior::Async:
            GMX_ASSERT(isHostMemoryPinned(hostBuffer),
                       "Destination host buffer was not pinned for HIP");
            stat = hipMemcpyAsync(hostBuffer,
                                  *reinterpret_cast<ValueType**>(buffer) + startingOffset,
                                  bytes,
                                  hipMemcpyDeviceToHost,
                                  deviceStream.stream());
            GMX_RELEASE_ASSERT(
                    stat == hipSuccess,
                    ("Asynchronous D2H copy failed. " + gmx::getDeviceErrorString(stat)).c_str());
            break;

        case GpuApiCallBehavior::Sync:
            stat = hipMemcpy(hostBuffer,
                             *reinterpret_cast<ValueType**>(buffer) + startingOffset,
                             bytes,
                             hipMemcpyDeviceToHost);
            GMX_RELEASE_ASSERT(
                    stat == hipSuccess,
                    ("Synchronous D2H copy failed. " + gmx::getDeviceErrorString(stat)).c_str());
            break;

        default: throw;
    }
}

/*! \brief
 * Performs the device-to-device data copy, synchronous or asynchronously on request.
 *
 * \tparam        ValueType                Raw value type of the \p buffer.
 * \param[in,out] destinationDeviceBuffer  Device-side buffer to copy to
 * \param[in]     sourceDeviceBuffer       Device-side buffer to copy from
 * \param[in]     numValues                Number of values to copy.
 * \param[in]     deviceStream             GPU stream to perform asynchronous copy in.
 * \param[in]     transferKind             Copy type: synchronous or asynchronous.
 */
template<typename ValueType>
void copyBetweenDeviceBuffers(DeviceBuffer<ValueType>* destinationDeviceBuffer,
                              DeviceBuffer<ValueType>* sourceDeviceBuffer,
                              size_t                   numValues,
                              const DeviceStream&      deviceStream,
                              GpuApiCallBehavior       transferKind,
                              CommandEvent* /*timingEvent*/)
{
    if (numValues == 0)
    {
        return;
    }
    GMX_ASSERT(destinationDeviceBuffer, "needs a destination buffer pointer");
    GMX_ASSERT(sourceDeviceBuffer, "needs a source buffer pointer");

    hipError_t   stat;
    const size_t bytes = numValues * sizeof(ValueType);
    switch (transferKind)
    {
        case GpuApiCallBehavior::Async:
            stat = hipMemcpyAsync(*destinationDeviceBuffer,
                                  *sourceDeviceBuffer,
                                  bytes,
                                  hipMemcpyDeviceToDevice,
                                  deviceStream.stream());
            GMX_RELEASE_ASSERT(
                    stat == hipSuccess,
                    ("Asynchronous D2D copy failed. " + gmx::getDeviceErrorString(stat)).c_str());
            break;

        case GpuApiCallBehavior::Sync:
            stat = hipMemcpy(*destinationDeviceBuffer, *sourceDeviceBuffer, bytes, hipMemcpyDeviceToDevice);
            GMX_RELEASE_ASSERT(
                    stat == hipSuccess,
                    ("Synchronous D2D copy failed. " + gmx::getDeviceErrorString(stat)).c_str());
            break;

        default: throw;
    }
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
    GMX_ASSERT(buffer, "needs a buffer pointer");
    const size_t bytes   = numValues * sizeof(ValueType);
    const char   pattern = 0;

    hipError_t stat = hipMemsetAsync(
            *reinterpret_cast<ValueType**>(buffer) + startingOffset, pattern, bytes, deviceStream.stream());

    GMX_RELEASE_ASSERT(stat == hipSuccess,
                       ("Couldn't clear the device buffer. " + gmx::getDeviceErrorString(stat)).c_str());
}

/*! \brief Check the validity of the device buffer.
 *
 * Checks if the buffer is not nullptr.
 *
 * \todo Add checks on the buffer size when it will be possible.
 *
 * \param[in] buffer        Device buffer to be checked.
 * \param[in] requiredSize  Number of elements that the buffer will have to accommodate.
 *
 * \returns Whether the device buffer can be set.
 */
template<typename T>
gmx_unused static bool checkDeviceBuffer(DeviceBuffer<T> buffer, gmx_unused int requiredSize)
{
    return buffer != nullptr;
}

//! Device texture wrapper.
using DeviceTexture = hipTextureObject_t;

/*! \brief Create a texture object for an array of type ValueType.
 *
 * Creates the device buffer, copies data and binds texture object for an array of type ValueType.
 *
 * \todo Test if using textures is still relevant on modern hardware.
 *
 * \tparam      ValueType      Raw data type.
 *
 * \param[out]  deviceBuffer   Device buffer to store data in.
 * \param[in]   hostBuffer     Host buffer to get date from
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
    if (numValues == 0)
    {
        return;
    }
    GMX_ASSERT(hostBuffer, "Host buffer should be specified.");

    allocateDeviceBuffer<ValueType>(deviceBuffer, numValues, deviceContext);

    const size_t sizeInBytes = numValues * sizeof(ValueType);

    hipError_t stat = hipMemcpy(
            *reinterpret_cast<ValueType**>(deviceBuffer), hostBuffer, sizeInBytes, hipMemcpyHostToDevice);

    GMX_RELEASE_ASSERT(stat == hipSuccess,
                       ("Synchronous H2D copy failed. " + gmx::getDeviceErrorString(stat)).c_str());
}

/*! \brief Unbind the texture and release the HIP texture object.
 *
 * \tparam         ValueType      Raw data type
 *
 * \param[in,out]  deviceBuffer   Device buffer to store data in.
 */
template<typename ValueType>
void destroyParamLookupTable(DeviceBuffer<ValueType>* deviceBuffer,
                             const DeviceTexture* /* deviceTexture */)
{
    freeDeviceBuffer(deviceBuffer);
}

template<typename ValueType>
ValueType* asMpiPointer(DeviceBuffer<ValueType>& buffer)
{
    return buffer;
}

#endif
