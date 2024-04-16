/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
#ifndef GMX_GPU_UTILS_DEVICEBUFFER_SYCL_H
#define GMX_GPU_UTILS_DEVICEBUFFER_SYCL_H

/*! \libinternal \file
 *  \brief Implements the DeviceBuffer type and routines for SYCL.
 *  Should only be included directly by the main DeviceBuffer file devicebuffer.h.
 *  TODO: the intent is for DeviceBuffer to become a class.
 *
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *  \author Erik Lindahl <erik.lindahl@gmail.com>
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *
 *  \inlibraryapi
 */

#include <utility>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/gpu_utils.h" //only for GpuApiCallBehavior
#include "gromacs/gpu_utils/gputraits_sycl.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#ifndef DOXYGEN

template<typename T>
class DeviceBuffer<T>::SyclBufferWrapper
{
public:
    T*            ptr_;
    sycl::context context_;
    SyclBufferWrapper(T* ptr, sycl::context context) : ptr_(ptr), context_(std::move(context)) {}
    operator T*() { return ptr_; }
};

template<typename T>
using SyclBufferWrapper = typename DeviceBuffer<T>::SyclBufferWrapper;

//! Constructor.
template<typename T>
DeviceBuffer<T>::DeviceBuffer() : buffer_(nullptr)
{
}

//! Dummy assignment operator to allow compilation of cross-platform code
template<typename T>
DeviceBuffer<T>::DeviceBuffer(std::nullptr_t nullPtr) : buffer_(nullPtr)
{
}

//! Destructor.
template<typename T>
DeviceBuffer<T>::~DeviceBuffer() = default;

//! Copy constructor (references the same underlying SYCL memory).
template<typename T>
DeviceBuffer<T>::DeviceBuffer(DeviceBuffer<T> const& src)
{
    if (src.buffer_)
    {
        buffer_ = std::make_unique<SyclBufferWrapper>(*src.buffer_);
    }
    else
    {
        buffer_ = nullptr;
    }
}

//! Move constructor.
template<typename T>
DeviceBuffer<T>::DeviceBuffer(DeviceBuffer<T>&& src) noexcept = default;

//! Copy assignment (references the same underlying SYCL buffer).
template<typename T>
DeviceBuffer<T>& DeviceBuffer<T>::operator=(DeviceBuffer<T> const& src)
{
    if (src.buffer_)
    {
        buffer_ = std::make_unique<SyclBufferWrapper>(*src.buffer_);
    }
    else
    {
        buffer_.reset(nullptr);
    }
    return *this;
}

//! Move assignment.
template<typename T>
DeviceBuffer<T>& DeviceBuffer<T>::operator=(DeviceBuffer<T>&& src) noexcept = default;

/*! \brief Dummy assignment operator to allow compilation of some cross-platform code.
 *
 * A hacky way to make SYCL implementation of DeviceBuffer compatible with details of CUDA and
 * OpenCL implementations.
 *
 * \todo Should be removed after DeviceBuffer refactoring.
 *
 * \tparam T Type of buffer content.
 * \param nullPtr \c std::nullptr. Not possible to assign any other pointers.
 */
template<typename T>
DeviceBuffer<T>& DeviceBuffer<T>::operator=(std::nullptr_t nullPtr)
{
    buffer_.reset(nullPtr);
    return *this;
}

//! Get the underlying device const pointer
template<typename T>
const T* DeviceBuffer<T>::get_pointer() const
{
    return buffer_ ? buffer_->ptr_ : nullptr;
}

//! Get the underlying device pointer
template<typename T>
T* DeviceBuffer<T>::get_pointer()
{
    return buffer_ ? buffer_->ptr_ : nullptr;
}

#endif // #ifndef DOXYGEN

/*! \brief Check the validity of the device buffer.
 *
 * Checks if the buffer is valid and if its allocation is big enough.
 *
 * \param[in] buffer        Device buffer to be checked.
 * \param[in] requiredSize  Number of elements that the buffer will have to accommodate.
 *
 * \returns Whether the device buffer exists and has enough capacity.
 */
template<typename T>
static gmx_unused bool checkDeviceBuffer(const DeviceBuffer<T>& buffer, int gmx_unused requiredSize)
{
    return buffer.buffer_ && buffer.buffer_->ptr_;
}

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
    ValueType* ptr = sycl::malloc_device<ValueType>(
            numValues, deviceContext.deviceInfo().syclDevice, deviceContext.context());
    buffer->buffer_.reset(new SyclBufferWrapper<ValueType>(ptr, deviceContext.context()));
}

/*! \brief
 * Frees a device-side buffer.
 * This does not reset separately stored size/capacity integers,
 * as this is planned to be a destructor of DeviceBuffer as a proper class,
 * and no calls on \p buffer should be made afterwards.
 *
 * \param[in] buffer  Pointer to the buffer to free.
 */
template<typename ValueType>
void freeDeviceBuffer(DeviceBuffer<ValueType>* buffer)
{
    if (buffer->buffer_ && buffer->buffer_->ptr_)
    {
        sycl::free(buffer->buffer_->ptr_, buffer->buffer_->context_);
    }
    buffer->buffer_.reset(nullptr);
}

/*! \brief
 * Performs the host-to-device data copy, synchronous or asynchronously on request.
 *
 * Unlike in CUDA and OpenCL, synchronous call does not guarantee that all previously
 * submitted operations are complete, only the ones that are required for \p buffer consistency.
 *
 * \tparam        ValueType            Raw value type of the \p buffer.
 * \param[in,out] buffer               Pointer to the device-side buffer.
 * \param[in]     hostBuffer           Pointer to the raw host-side memory, also typed \p ValueType.
 * \param[in]     startingOffset       Offset (in values) at the device-side buffer to copy into.
 * \param[in]     numValues            Number of values to copy.
 * \param[in]     deviceStream         GPU stream to perform asynchronous copy in.
 * \param[in]     transferKind         Copy type: synchronous or asynchronous.
 * \param[out]    timingEvent          A pointer to the H2D copy timing event to be filled in.
 *                                     Ignored in SYCL.
 */
template<typename ValueType>
void copyToDeviceBuffer(DeviceBuffer<ValueType>* buffer,
                        const ValueType*         hostBuffer,
                        size_t                   startingOffset,
                        size_t                   numValues,
                        const DeviceStream&      deviceStream,
                        GpuApiCallBehavior       transferKind,
                        CommandEvent* gmx_unused timingEvent)
{
    if (numValues == 0)
    {
        return; // such calls are actually made with empty domains
    }
    GMX_ASSERT(buffer, "needs a buffer pointer");
    GMX_ASSERT(hostBuffer, "needs a host buffer pointer");

    GMX_ASSERT(checkDeviceBuffer(*buffer, startingOffset + numValues),
               "buffer too small or not initialized");

    if (transferKind == GpuApiCallBehavior::Async)
    {
        using sycl::usm::alloc;
        GMX_ASSERT(sycl::get_pointer_type(hostBuffer, buffer->buffer_->context_) == alloc::host,
                   "Trying to launch async copy from unpinned host buffer");
    }

    ValueType*   dstPtr = buffer->buffer_->ptr_ + startingOffset;
    const size_t size   = numValues * sizeof(ValueType);
    if (transferKind == GpuApiCallBehavior::Sync)
    {
        deviceStream.stream()
                .submit([&](sycl::handler& cgh) { cgh.memcpy(dstPtr, hostBuffer, size); })
                .wait_and_throw();
    }
    else
    {
        deviceStream.stream().submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
            cgh.memcpy(dstPtr, hostBuffer, size);
        });
    }
}

/*! \brief
 * Performs the device-to-host data copy, synchronous or asynchronously on request.
 *
 * Unlike in CUDA and OpenCL, synchronous call does not guarantee that all previously
 * submitted operations are complete, only the ones that are required for \p buffer consistency.
 *
 * \tparam        ValueType            Raw value type of the \p buffer.
 * \param[in,out] hostBuffer           Pointer to the raw host-side memory, also typed \p ValueType
 * \param[in]     buffer               Pointer to the device-side buffer.
 * \param[in]     startingOffset       Offset (in values) at the device-side buffer to copy from.
 * \param[in]     numValues            Number of values to copy.
 * \param[in]     deviceStream         GPU stream to perform asynchronous copy in.
 * \param[in]     transferKind         Copy type: synchronous or asynchronous.
 * \param[out]    timingEvent          A pointer to the H2D copy timing event to be filled in.
 *                                     Ignored in SYCL.
 */
template<typename ValueType>
void copyFromDeviceBuffer(ValueType*               hostBuffer,
                          DeviceBuffer<ValueType>* buffer,
                          size_t                   startingOffset,
                          size_t                   numValues,
                          const DeviceStream&      deviceStream,
                          GpuApiCallBehavior       transferKind,
                          CommandEvent* gmx_unused timingEvent)
{
    if (numValues == 0)
    {
        return; // such calls are actually made with empty domains
    }
    GMX_ASSERT(buffer, "needs a buffer pointer");
    GMX_ASSERT(hostBuffer, "needs a host buffer pointer");

    GMX_ASSERT(checkDeviceBuffer(*buffer, startingOffset + numValues),
               "buffer too small or not initialized");

    if (transferKind == GpuApiCallBehavior::Async)
    {
        using sycl::usm::alloc;
        GMX_ASSERT(sycl::get_pointer_type(hostBuffer, buffer->buffer_->context_) == alloc::host,
                   "Trying to launch async copy to unpinned host buffer");
    }

    const ValueType* srcPtr = buffer->buffer_->ptr_ + startingOffset;
    const size_t     size   = numValues * sizeof(ValueType);
    if (transferKind == GpuApiCallBehavior::Sync)
    {
        deviceStream.stream()
                .submit([&](sycl::handler& cgh) { cgh.memcpy(hostBuffer, srcPtr, size); })
                .wait_and_throw();
    }
    else
    {
        deviceStream.stream().submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
            cgh.memcpy(hostBuffer, srcPtr, size);
        });
    }
}

/*! \brief
 * Performs the device-to-device data copy, synchronous or asynchronously on request.
 *
 * \tparam        ValueType                Raw value type of the \p buffer.
 */
template<typename ValueType>
void copyBetweenDeviceBuffers(DeviceBuffer<ValueType>* destinationDeviceBuffer,
                              DeviceBuffer<ValueType>* sourceDeviceBuffer,
                              size_t                   numValues,
                              const DeviceStream&      deviceStream,
                              GpuApiCallBehavior       transferKind,
                              CommandEvent* gmx_unused timingEvent)
{
    if (numValues == 0)
    {
        return; // such calls are actually made with empty domains
    }
    GMX_ASSERT(destinationDeviceBuffer, "needs a destination buffer pointer");
    GMX_ASSERT(sourceDeviceBuffer, "needs a source buffer pointer");

    const ValueType* srcPtr = sourceDeviceBuffer->buffer_->ptr_;
    ValueType*       dstPtr = destinationDeviceBuffer->buffer_->ptr_;
    const size_t     size   = numValues * sizeof(ValueType);
    if (transferKind == GpuApiCallBehavior::Sync)
    {
        deviceStream.stream().submit([&](sycl::handler& cgh) { cgh.memcpy(dstPtr, srcPtr, size); }).wait_and_throw();
    }
    else
    {
        deviceStream.stream().submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
            cgh.memcpy(dstPtr, srcPtr, size);
        });
    }
}

/*! \brief
 * Clears the device buffer asynchronously.
 *
 * \tparam        ValueType       Raw value type of the \p buffer.
 * \param[in,out] buffer          Pointer to the device-side buffer.
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

    GMX_ASSERT(checkDeviceBuffer(*buffer, startingOffset + numValues),
               "buffer too small or not initialized");

    deviceStream.stream().submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        cgh.memset(buffer->buffer_->ptr_ + startingOffset, 0, numValues * sizeof(ValueType));
    });
}

/*! \brief Create a texture object for an array of type ValueType.
 *
 * Creates the device buffer and copies read-only data for an array of type ValueType.
 * Like OpenCL, does not really do anything with textures, simply creates a buffer
 * and initializes it.
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
    GMX_ASSERT(hostBuffer, "Host buffer should be specified.");
    GMX_ASSERT(deviceBuffer, "Device buffer should be specified.");

    allocateDeviceBuffer<ValueType>(deviceBuffer, numValues, deviceContext);
    /* Not perfect, but we call this function only on simulation initialization, so the
     * overhead of a queue creation should be manageable. */
    DeviceStream temporaryStream(deviceContext, DeviceStreamPriority::Normal, false);
    copyToDeviceBuffer(
            deviceBuffer, hostBuffer, 0, numValues, temporaryStream, GpuApiCallBehavior::Sync, nullptr);
}

/*! \brief Release the underlying device allocations.
 *
 * \tparam        ValueType     Raw data type.
 *
 * \param[in,out] deviceBuffer  Device buffer to store data in.
 */
template<typename ValueType>
void destroyParamLookupTable(DeviceBuffer<ValueType>* deviceBuffer, DeviceTexture* /* deviceTexture */)
{
    freeDeviceBuffer(deviceBuffer);
}

template<typename ValueType>
ValueType* asMpiPointer(DeviceBuffer<ValueType>& buffer)
{
    return buffer.get_pointer();
}

#endif // GMX_GPU_UTILS_DEVICEBUFFER_SYCL_H
