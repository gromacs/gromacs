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

/*! \brief Helper function to submit a SYCL operation without returning an event.
 *
 * Gives some nice performance optimizations, especially on AMD and NVIDIA devices.
 *
 * In ACpp, it relies on the ACPP_EXT_CG_PROPERTY_* and ACPP_EXT_COARSE_GRAINED_EVENTS extensions.
 *
 * In DPC++, it relies (when explicitly enabled) on SYCL_EXT_ONEAPI_ENQUEUE_FUNCTIONS
 * Falls back to the default submit otherwise.
 *
 * The function aims to avoid the overhead associated with creating/recording/destroying events.
 */
template<typename Queue, typename CommandGroupFunc>
static inline void syclSubmitWithoutEvent(Queue&& queue, CommandGroupFunc&& cgf)
{
#if GMX_SYCL_ACPP
    queue.submit(gmx::internal::sc_syclDiscardEventProperty_list, std::move(cgf));
#elif defined(SYCL_EXT_ONEAPI_ENQUEUE_FUNCTIONS) && GMX_SYCL_ENABLE_EXPERIMENTAL_SUBMIT_API
    sycl::ext::oneapi::experimental::submit(queue, std::move(cgf));
#else
    queue.submit(std::move(cgf));
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
