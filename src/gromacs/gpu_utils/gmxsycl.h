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
/*! \libinternal \file
 * \brief
 * Wraps the complexity of including SYCL in GROMACS.
 *
 * \inlibraryapi
 */

#ifndef GMX_GPU_UTILS_GMXSYCL_H
#define GMX_GPU_UTILS_GMXSYCL_H

#include "config.h"

#include <type_traits>

#include <sycl/sycl.hpp>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{
#if GMX_SYCL_ACPP
namespace internal
{
static const sycl::property_list sc_syclDiscardEventProperty_list{
#    if defined(ACPP_EXT_COARSE_GRAINED_EVENTS) // Since ACpp 24.06
    sycl::property::command_group::AdaptiveCpp_coarse_grained_events()
#    else
    sycl::property::command_group::hipSYCL_coarse_grained_events()
#    endif
};
} // namespace internal
#endif

#if GMX_SYCL_ENABLE_HANDLER_FREE_SUBMISSION
using CommandGroupHandler = std::nullptr_t;
#else
using CommandGroupHandler = sycl::handler&;
#endif

namespace internal
{

/*! \brief Internal helper function to submit to a queue via a free
 * function, avoiding returning an event.
 */
template<typename KernelNameType, typename Queue, typename Kernel, typename Range>
static inline void syclSubmitWithFreeFunction(Queue&& queue, Kernel&& kernel, Range&& range)
{
#if GMX_SYCL_ENABLE_HANDLER_FREE_SUBMISSION
    // GROMACS sometimes expresses kernel work ranges with sycl::range
    // and other times with sycl::nd_range. These map to handler-free
    // submission functions with different names, so we use the type
    // and dimensions of the range to decide which function to call.
    using UnderlyingRange       = std::remove_const_t<std::remove_reference_t<Range>>;
    constexpr int sc_dimensions = UnderlyingRange::dimensions;
    if constexpr (std::is_same_v<UnderlyingRange, sycl::range<sc_dimensions>>)
    {
        sycl::ext::oneapi::experimental::parallel_for<KernelNameType>(queue, range, kernel);
    }
    else if constexpr (std::is_same_v<UnderlyingRange, sycl::nd_range<sc_dimensions>>)
    {
        sycl::ext::oneapi::experimental::nd_launch<KernelNameType>(queue, range, kernel);
    }
    else
    {
        static_assert(std::is_same_v<Range, std::nullptr_t>, "Unknown range type!");
    }
#else
    GMX_RELEASE_ASSERT(false, "Function should not be called");
    GMX_UNUSED_VALUE(queue);
    GMX_UNUSED_VALUE(kernel);
    GMX_UNUSED_VALUE(range);
#endif
}

} // namespace internal

/*! \brief Helper function to submit to a queue via a handler,
 * avoiding returning an event (if possible).
 */
template<typename Queue, typename CommandGroupFunction>
static inline void syclSubmitWithCghWithoutEvent(Queue&& queue, CommandGroupFunction&& commandGroupFunction)
{
#if GMX_SYCL_ENABLE_HANDLER_FREE_SUBMISSION
    GMX_RELEASE_ASSERT(false, "Function should not be called");
    GMX_UNUSED_VALUE(queue);
    GMX_UNUSED_VALUE(commandGroupFunction);
#elif GMX_SYCL_ENABLE_EXPERIMENTAL_SUBMIT_API
    sycl::ext::oneapi::experimental::submit(queue, std::move(commandGroupFunction));
#elif GMX_SYCL_ACPP
    queue.submit(internal::sc_syclDiscardEventProperty_list, std::move(commandGroupFunction));
#else
    // Generic fallback that actually does make and ignore a sycl::event
    queue.submit(std::move(commandGroupFunction));
#endif
}

/*! \brief Helper function to submit a SYCL operation without
 * returning an event (if possible).
 *
 * Gives some nice performance optimizations, especially on AMD and NVIDIA devices.
 *
 * \param[in]  queue                  The sycl::queue to use
 * \param[in]  kernelFunctionBuilder  A function that returns a device kernel as a lambda
 * \param[in]  range                  The range over which the device kernel runs
 * \param[in]  args                   Parameter pack of device-kernel arguments
 *
 * In ACpp, it relies on the ACPP_EXT_CG_PROPERTY_* and ACPP_EXT_COARSE_GRAINED_EVENTS extensions.
 *
 * Two implementations are provided for DPC++. Both rely (when
 * explicitly enabled) on SYCL_EXT_ONEAPI_ENQUEUE_FUNCTIONS.  The
 * GMX_SYCL_ENABLE_EXPERIMENTAL_SUBMIT_API implementation continues to
 * use a command-group-handler, but no longer returns a
 * sycl::event. The GMX_SYCL_ENABLE_HANDLER_FREE_SUBMISSION implementation
 * uses free functions to submit work on a queue, requiring no handler
 * and permitting the caller to choose whether to receive a
 * sycl::event as a return value.
 *
 * Falls back to the default submit otherwise.
 *
 * The function aims to avoid the overhead associated with
 * creating/recording/destroying events and (depending on the chosen
 * implementation) command-group handlers, while always permitting
 * fallback standard execution if other capabilities are not
 * available. */
template<typename KernelNameType, typename KernelFunctionBuilder, typename Queue, typename Range, class... Args>
static inline void syclSubmitWithoutEvent(Queue&&               queue,
                                          KernelFunctionBuilder kernelFunctionBuilder,
                                          Range&&               range,
                                          Args&&... args)
{
#if GMX_SYCL_ENABLE_HANDLER_FREE_SUBMISSION
    // Make a dummy command-group handler so that SYCL kernel-function
    // builders have the same call signature, whether or not the
    // command-group handler will actually be used.
    CommandGroupHandler dummyCgh = nullptr;
    auto                kernel   = kernelFunctionBuilder(dummyCgh, std::forward<Args>(args)...);
    internal::syclSubmitWithFreeFunction<KernelNameType>(queue, kernel, range);
#else
    // Use a proper command-group handler and pass it to the
    // kernel-function builder.
    auto commandGroupFunction = [&](sycl::handler& cgh)
    {
        auto kernel = kernelFunctionBuilder(cgh, std::forward<Args>(args)...);
        cgh.parallel_for<KernelNameType>(range, kernel);
    };
    syclSubmitWithCghWithoutEvent(queue, std::move(commandGroupFunction));
#endif
}

/*! \brief Helper function to submit a SYCL operation without returning an event.
 *
 * Differs from \c syclSubmitWithoutEvent() in that if the submission
 * protocol must use a command-group handler, it is not passed to the
 * \c kernelFunctionBuilder.
 */
template<typename KernelNameType, typename KernelFunctionBuilder, typename Queue, typename Range, class... Args>
static inline void syclSubmitWithoutCghOrEvent(Queue&&               queue,
                                               KernelFunctionBuilder kernelFunctionBuilder,
                                               Range&&               range,
                                               Args&&... args)
{
#if GMX_SYCL_ENABLE_HANDLER_FREE_SUBMISSION
    auto kernel = kernelFunctionBuilder(std::forward<Args>(args)...);
    internal::syclSubmitWithFreeFunction<KernelNameType>(queue, kernel, range);
#else
    // Generic fallback that does use a proper command-group handler
    // as part of submission (despite the function name), but does not
    // use it in the call to the kernel-function builder.
    auto commandGroupFunction = [&](sycl::handler& cgh)
    {
        auto kernel = kernelFunctionBuilder(std::forward<Args>(args)...);
        cgh.parallel_for<KernelNameType>(range, kernel);
    };
    syclSubmitWithCghWithoutEvent(queue, std::move(commandGroupFunction));
#endif
}

/*! \brief Helper function to submit a SYCL async memcpy without
 * returning an event (if possible).
 *
 * See \c syclSubmitWithoutEvent() for general details. */
template<typename Queue>
static inline void syclMemcpyWithoutEvent(Queue&& queue, void* destination, const void* source, size_t numBytes)
{
#if GMX_SYCL_ENABLE_HANDLER_FREE_SUBMISSION
    sycl::ext::oneapi::experimental::memcpy(queue, destination, source, numBytes);
#else
    // Use a proper command-group handler
    auto commandGroupFunction = [&](sycl::handler& cgh)
    { cgh.memcpy(destination, source, numBytes); };
    syclSubmitWithCghWithoutEvent(queue, std::move(commandGroupFunction));
#endif
}

/*! \brief Helper function to submit a SYCL async memset without
 * returning an event (if possible).
 *
 * See \c syclSubmitWithoutEvent() for general details. */
static inline void syclMemsetWithoutEvent(sycl::queue&& queue, void* pointer, int value, size_t numBytes)
{
#if GMX_SYCL_ENABLE_HANDLER_FREE_SUBMISSION
    sycl::ext::oneapi::experimental::memset(queue, pointer, value, numBytes);
#else
    // Use a proper command-group handler
    auto commandGroupFunction = [&](sycl::handler& cgh) { cgh.memset(pointer, value, numBytes); };
    syclSubmitWithCghWithoutEvent(queue, std::move(commandGroupFunction));
#endif
}

/*! \brief Helper function to add a custom operation to the SYCL handler.
 *
 * In ACpp, it relies on the ACPP_EXT_ENQUEUE_CUSTOM_OPERATION extension.
 * Should not be called when the extension is not available.
 */
template<typename CommandGroupFunc>
static inline void syclEnqueueCustomOp(sycl::handler& cgh, CommandGroupFunc&& cgf)
{
#if defined(ACPP_EXT_ENQUEUE_CUSTOM_OPERATION)
    cgh.AdaptiveCpp_enqueue_custom_operation(std::move(cgf));
#elif defined(HIPSYCL_EXT_ENQUEUE_CUSTOM_OPERATION)
    cgh.hipSYCL_enqueue_custom_operation(std::move(cgf));
#elif defined(SYCL_EXT_ONEAPI_ENQUEUE_NATIVE_COMMAND) && GMX_SYCL_ENABLE_EXPERIMENTAL_SUBMIT_API
    cgh.ext_codeplay_enqueue_native_command([=](sycl::interop_handle h) { cgf(h); });
#else
    GMX_UNUSED_VALUE(cgh);
    GMX_UNUSED_VALUE(cgf);
    GMX_RELEASE_ASSERT(false, "Function called with unsupported backend");
#endif
}

} // namespace gmx

#endif
