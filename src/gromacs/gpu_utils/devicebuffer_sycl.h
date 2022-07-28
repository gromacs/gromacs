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

#    if !GMX_SYCL_USE_USM
template<typename T>
class DeviceBuffer<T>::SyclBufferWrapper : public sycl::buffer<T, 1>
{
    using sycl::buffer<T, 1>::buffer; // Get all the constructors
};
#    else
template<typename T>
class DeviceBuffer<T>::SyclBufferWrapper
{
public:
    T*            ptr_;
    sycl::context context_;
    SyclBufferWrapper(T* ptr, sycl::context context) : ptr_(ptr), context_(std::move(context)) {}
    operator T*() { return ptr_; }
};
#    endif

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

#    if GMX_SYCL_USE_USM

/** \brief
 * Thin wrapper around placeholder accessor that allows implicit construction from \c DeviceBuffer.
 *
 * "Placeholder accessor" is an indicator of the intent to create an accessor for certain buffer
 * of a certain type, that is not yet bound to a specific command group handler (device). Such
 * accessors can be created outside SYCL kernels, which is helpful if we want to pass them as
 * function arguments.
 *
 * \tparam T Type of buffer content.
 * \tparam mode Access mode.
 */
template<class T, sycl::access_mode mode>
class DeviceAccessor
{
public:
    //! Construct DeviceAccessor from DeviceBuffer (must be initialized)
    DeviceAccessor(DeviceBuffer<T>& buffer) : ptr_(getPointer(buffer)) {}
    //! Construct read-only DeviceAccessor from a const DeviceBuffer (must be initialized)
    DeviceAccessor(const DeviceBuffer<T>& buffer) : ptr_(getPointer(buffer))
    {
        static_assert(mode == sycl::access_mode::read,
                      "Can not create non-read-only accessor from a const DeviceBuffer");
    }
    //! Construct DeviceAccessor from a raw pointer
    DeviceAccessor(T* ptr) : ptr_(ptr) {}
    __attribute__((always_inline)) void bind(sycl::handler& /*cgh*/)
    {
        // Do nothing
    }
    using ValueType = std::conditional_t<mode == sycl::access_mode::read, const T, T>;
    __attribute__((always_inline)) ValueType* get_pointer() const noexcept { return ptr_; }
    __attribute__((always_inline)) ValueType& operator[](sycl::id<1> index) const
    {
        return ptr_[index.get(0)];
    }
    __attribute__((always_inline)) ValueType& operator[](size_t index) const { return ptr_[index]; }

private:
    //! Helper function to get sycl:global_ptr object from DeviceBuffer wrapper, with a sanity check.
    static inline sycl::global_ptr<T> getPointer(const DeviceBuffer<T>& buffer)
    {
        GMX_ASSERT(bool(buffer), "Trying to use an uninitialized buffer");
        return buffer.buffer_->ptr_;
    }
    T* ptr_;
};
#    else
namespace gmx::internal
{
//! Shorthand alias to create a placeholder SYCL accessor with chosen data type and access mode.
template<class T, sycl::access_mode mode>
using PlaceholderAccessor =
        sycl::accessor<T, 1, mode, sycl::target::global_buffer, sycl::access::placeholder::true_t>;
} // namespace gmx::internal

/** \brief
 * Thin wrapper around placeholder accessor that allows implicit construction from \c DeviceBuffer.
 *
 * "Placeholder accessor" is an indicator of the intent to create an accessor for certain buffer
 * of a certain type, that is not yet bound to a specific command group handler (device). Such
 * accessors can be created outside SYCL kernels, which is helpful if we want to pass them as
 * function arguments.
 *
 * \tparam T Type of buffer content.
 * \tparam mode Access mode.
 */
template<class T, sycl::access_mode mode>
class DeviceAccessor : public gmx::internal::PlaceholderAccessor<T, mode>
{
public:
    // Inherit all the constructors
    using gmx::internal::PlaceholderAccessor<T, mode>::PlaceholderAccessor;
    //! Construct Accessor from DeviceBuffer (must be initialized)
    DeviceAccessor(DeviceBuffer<T>& buffer) :
        gmx::internal::PlaceholderAccessor<T, mode>(getSyclBuffer(buffer))
    {
    }
    //! Construct read-only Accessor from a const DeviceBuffer (must be initialized)
    DeviceAccessor(const DeviceBuffer<T>& buffer) :
        gmx::internal::PlaceholderAccessor<T, mode>(getSyclBuffer(const_cast<DeviceBuffer<T>&>(buffer)))
    {
        /* There were some discussions about making it possible to create read-only sycl::accessor
         * from a const sycl::buffer (https://github.com/KhronosGroup/SYCL-Docs/issues/10), but
         * it did not make it into the SYCL2020 standard. So, we have to use const_cast above */
        /* Using static_assert to ensure that only mode::read accessors can be created from a
         * const DeviceBuffer. static_assert provides better error messages than std::enable_if. */
        static_assert(mode == sycl::access_mode::read,
                      "Can not create non-read-only accessor from a const DeviceBuffer");
    }
    void bind(sycl::handler& cgh) { cgh.require(*this); }

private:
    //! Helper function to get sycl:buffer object from DeviceBuffer wrapper, with a sanity check.
    static inline sycl::buffer<T, 1>& getSyclBuffer(DeviceBuffer<T>& buffer)
    {
        GMX_ASSERT(bool(buffer), "Trying to construct accessor from an uninitialized buffer");
        return *buffer.buffer_;
    }
};
#    endif // GMX_SYCL_USE_USM

namespace gmx::internal
{
//! A non-functional class that can be used instead of real accessors
template<class T>
struct NullAccessor
{
    NullAccessor(const DeviceBuffer<T>& /*buffer*/) {}
    //! Allow casting to nullptr
    constexpr operator std::nullptr_t() const { return nullptr; }
    //! Placeholder implementation of \c sycl::accessor::get_pointer.
    T* get_pointer() const noexcept { return nullptr; }
    //! Placeholder for \c sycl::handler::require.
    void bind(sycl::handler& /*cgh*/) { assert(false); }
};

} // namespace gmx::internal

/** \brief
 * Helper class to be used as function argument. Will either correspond to a device accessor, or an empty class.
 *
 * Example usage:
 * \code
    template <bool doFoo>
    void getBarKernel(handler& cgh, OptionalAccessor<float, mode::read, doFoo> a_fooPrms)
    {
        if constexpr (doFoo)
            a_fooPrms.bind(cgh);
        // Can only use a_fooPrms if doFoo == true
    }

    template <bool doFoo>
    void callBar(DeviceBuffer<float> b_fooPrms)
    {
        // If doFoo is false, b_fooPrms will be ignored (can be not initialized).
        // Otherwise, an accessor will be built (b_fooPrms must be a valid buffer).
        auto kernel = getBarKernel<doFoo>(b_fooPrms);
        // If the accessor in not enabled, anything can be passed as its ctor argument.
        auto kernel2 = getBarKernel<false>(nullptr_t);
    }
 * \endcode
 *
 * \tparam T Data type of the underlying buffer
 * \tparam mode Access mode of the accessor
 * \tparam enabled Compile-time flag indicating whether we want to actually create an accessor.
 */
template<class T, sycl::access_mode mode, bool enabled>
using OptionalAccessor =
        std::conditional_t<enabled, DeviceAccessor<T, mode>, gmx::internal::NullAccessor<T>>;

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
static gmx_unused bool checkDeviceBuffer(const DeviceBuffer<T>& buffer, int requiredSize)
{
#if GMX_SYCL_USE_USM
    GMX_UNUSED_VALUE(requiredSize); // Can't be checked with USM
    return buffer.buffer_ && buffer.buffer_->ptr_;
#else
    return buffer.buffer_ && (static_cast<int>(buffer.buffer_->size()) >= requiredSize);
#endif
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
#if GMX_SYCL_USE_USM
    ValueType* ptr = sycl::malloc_device<ValueType>(
            numValues, deviceContext.deviceInfo().syclDevice, deviceContext.context());
    buffer->buffer_.reset(new SyclBufferWrapper<ValueType>(ptr, deviceContext.context()));
#else
    /* SYCL does not require binding buffer to a specific context or device. The ::context_bound
     * property only enforces the use of only given context, and possibly offers some optimizations */
    const sycl::property_list bufferProperties{ sycl::property::buffer::context_bound(
            deviceContext.context()) };
    buffer->buffer_.reset(new SyclBufferWrapper<ValueType>(sycl::range<1>(numValues), bufferProperties));
#endif
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
#if GMX_SYCL_USE_USM
    if (buffer->buffer_ && buffer->buffer_->ptr_)
    {
        sycl::free(buffer->buffer_->ptr_, buffer->buffer_->context_);
    }
#endif
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

    sycl::event ev;
#if GMX_SYCL_USE_USM
    ev = deviceStream.stream().submit([&](sycl::handler& cgh) {
        cgh.memcpy(buffer->buffer_->ptr_ + startingOffset, hostBuffer, numValues * sizeof(ValueType));
    });
#else
    sycl::buffer<ValueType>& syclBuffer = *buffer->buffer_;

    ev                                  = deviceStream.stream().submit([&](sycl::handler& cgh) {
        /* Here and elsewhere in this file, accessor constructor is user instead of a more common
         * buffer::get_access, since the compiler (icpx 2021.1-beta09) occasionally gets confused
         * by all the overloads */
        auto d_bufferAccessor = sycl::accessor<ValueType, 1, sycl::access_mode::write>{
            syclBuffer, cgh, sycl::range(numValues), sycl::id(startingOffset), { sycl::no_init }
        };
        cgh.copy(hostBuffer, d_bufferAccessor);
    });
#endif
    if (transferKind == GpuApiCallBehavior::Sync)
    {
        ev.wait_and_throw();
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

    sycl::event ev;
#if GMX_SYCL_USE_USM
    ev = deviceStream.stream().submit([&](sycl::handler& cgh) {
        cgh.memcpy(hostBuffer, buffer->buffer_->ptr_ + startingOffset, numValues * sizeof(ValueType));
    });
#else
    sycl::buffer<ValueType>& syclBuffer = *buffer->buffer_;

    ev = deviceStream.stream().submit([&](sycl::handler& cgh) {
        const auto d_bufferAccessor = sycl::accessor<ValueType, 1, sycl::access_mode::read>{
            syclBuffer, cgh, sycl::range(numValues), sycl::id(startingOffset)
        };
        cgh.copy(d_bufferAccessor, hostBuffer);
    });
#endif
    if (transferKind == GpuApiCallBehavior::Sync)
    {
        ev.wait_and_throw();
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
    // SYCL-TODO
    gmx_fatal(FARGS, "D2D copy stub was called. Not yet implemented in SYCL.");
}


namespace gmx::internal
{
/*! \brief Helper function to clear device buffer.
 *
 * Not applicable to GROMACS's Float3 (a.k.a. gmx::RVec) and other custom types.
 * From SYCL specs: "T must be a scalar value or a SYCL vector type."
 */
template<typename ValueType>
sycl::event fillSyclBufferWithNull(sycl::buffer<ValueType, 1>& buffer,
                                   size_t                      startingOffset,
                                   size_t                      numValues,
                                   sycl::queue                 queue)
{
    using sycl::access_mode;
    const sycl::range<1> range(numValues);
    const sycl::id<1>    offset(startingOffset);
    const ValueType      pattern = ValueType(0); // SYCL vectors support initialization by scalar

    return queue.submit([&](sycl::handler& cgh) {
        auto d_bufferAccessor =
                sycl::accessor<ValueType, 1, access_mode::write>{ buffer, cgh, range, offset, sycl::no_init };
        cgh.fill(d_bufferAccessor, pattern);
    });
}

//! \brief Helper function to clear device buffer of type Float3.
template<>
inline sycl::event fillSyclBufferWithNull(sycl::buffer<Float3, 1>& buffer,
                                          size_t                   startingOffset,
                                          size_t                   numValues,
                                          sycl::queue              queue)
{
    constexpr bool usingHipSycl =
#ifdef __HIPSYCL__
            true;
#else
            false;
#endif

    if constexpr (usingHipSycl)
    {
        // hipSYCL does not support reinterpret but allows using Float3 directly.
        using sycl::access_mode;
        const sycl::range<1> range(numValues);
        const sycl::id<1>    offset(startingOffset);
        const Float3         pattern{ 0, 0, 0 };

        return queue.submit([&](sycl::handler& cgh) {
            auto d_bufferAccessor =
                    sycl::accessor<Float3, 1, access_mode::write>{ buffer, cgh, range, offset, sycl::no_init };
            cgh.fill(d_bufferAccessor, pattern);
        });
    }
    else // When not using hipSYCL, reinterpret as a flat float array
    {
#ifndef __HIPSYCL__
        sycl::buffer<float, 1> bufferAsFloat = buffer.reinterpret<float, 1>(buffer.size() * DIM);
        return fillSyclBufferWithNull<float>(
                bufferAsFloat, startingOffset * DIM, numValues * DIM, std::move(queue));
#endif
    }
}

//! \brief Helper function to clear device memory
template<typename ValueType>
sycl::event fillSyclUsmWithNull(ValueType* buffer, size_t startingOffset, size_t numValues, sycl::queue queue)
{
    return queue.memset(buffer + startingOffset, 0, numValues * sizeof(ValueType));
}

} // namespace gmx::internal

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

#if GMX_SYCL_USE_USM
    gmx::internal::fillSyclUsmWithNull<ValueType>(
            buffer->buffer_->ptr_, startingOffset, numValues, deviceStream.stream());
#else
    sycl::buffer<ValueType>& syclBuffer = *(buffer->buffer_);

    gmx::internal::fillSyclBufferWithNull<ValueType>(
            syclBuffer, startingOffset, numValues, deviceStream.stream());
#endif
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

#if GMX_SYCL_USE_USM
    allocateDeviceBuffer<ValueType>(deviceBuffer, numValues, deviceContext);
    /* Not perfect, but we call this function only on simulation initialization, so the
     * overhead of a queue creation should be manageable. */
    DeviceStream temporaryStream(deviceContext, DeviceStreamPriority::Normal, false);
    copyToDeviceBuffer(
            deviceBuffer, hostBuffer, 0, numValues, temporaryStream, GpuApiCallBehavior::Sync, nullptr);
#else
    /* Constructing buffer with sycl::buffer(T* data, size_t size) will take ownership
     * of this memory region making it unusable, which might lead to side-effects.
     * On the other hand, sycl::buffer(InputIterator<T> begin, InputIterator<T> end) will
     * initialize the buffer without affecting ownership of the memory, although
     * it will consume extra memory on host. */
    const sycl::property_list bufferProperties{ sycl::property::buffer::context_bound(
            deviceContext.context()) };
    deviceBuffer->buffer_.reset(
            new SyclBufferWrapper<ValueType>(hostBuffer, hostBuffer + numValues, bufferProperties));
#endif
}

/*! \brief Release the OpenCL device buffer.
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

#endif // GMX_GPU_UTILS_DEVICEBUFFER_SYCL_H
