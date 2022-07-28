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
/*! \internal \file
 *  \brief Pinned memory allocation routines for SYCL
 *
 * This module is a direct wrapper around the sycl::malloc_host / sycl::free, except for the
 * management of the default context.
 *
 * Unlike in CUDA, pinning memory in SYCL requires a context. It can be passed explicitly to
 * \c pmalloc and \c pfree, but that is not straightforward in some calling code. Therefore,
 * we use \c pmallocSetDefaultDeviceContext and \c pmallocClearDefaultDeviceContext to manage
 * the default context (stored in \c g_threadDefaultContext).
 *
 * That puts a constraint on allocation order: we shall free all the pinned memory before resetting
 * the default context. This is easy to achieve in normal use, but hard to guarantee during
 * stack unwinding when handling an exception. Therefore, we introduce \c g_threadAllocCount to
 * count the number of allocations that are using the default context. If \c pmallocClearDefaultDeviceContext
 * is called while handling an exception, we check \c g_threadAllocCount, and, if there are any
 * remaining allocations, set \c g_threadDelayContextClearing to defer the context resetting until
 * all the allocations are freed. We also use \c g_threadAllocCount to control the correctness of
 * memory allocations in normal runs.
 *
 * GROMACS (at least in 2022 and earlier) uses separate contexts for each rank. Since we support
 * threadMPI, the context management is per-thread, and all the static variables are \c thread_local.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 */
#include "gmxpre.h"

#include <exception>
#include <optional>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "pmalloc.h"

//! Default context to use for pinning memory.
static thread_local std::optional<const sycl::context> g_threadDefaultContext = std::nullopt;
//! Count the number of memory allocations in the default context.
static thread_local int g_threadAllocCount = 0;
//! Whether we should delay resetting the default context because there is still memory allocated there.
static thread_local bool g_threadDelayContextClearing = false;

/*! \brief Allocates \p nbytes of host memory. Use \c pfree to free memory allocated with this function.
 *
 * \param[in,out] h_ptr   Pointer where to store the address of the newly allocated buffer.
 * \param[in]     nbytes  Size in bytes of the buffer to be allocated.
 * \param[in]     deviceContext SYCL context to use. Will use the default one (see \c pmallocSetDefaultDeviceContext) if not set.
 */
void pmalloc(void** h_ptr, size_t nbytes, const DeviceContext* deviceContext)
{
    GMX_RELEASE_ASSERT(deviceContext || g_threadDefaultContext.has_value(),
                       "Calling pmalloc without a context");
    if (!deviceContext)
    {
        g_threadAllocCount++;
    }
    auto** h_ptrAsBytes = reinterpret_cast<std::byte**>(h_ptr);

    *h_ptrAsBytes = sycl::malloc_host<std::byte>(
            nbytes, deviceContext ? deviceContext->context() : *g_threadDefaultContext);
}

/*! \brief Frees memory allocated with pmalloc.
 *
 * \param[in] h_ptr         Buffer allocated with \c pmalloc that needs to be freed.
 * \param[in] deviceContext SYCL context to use. Will use the default one (see \c pmallocSetDefaultDeviceContext) if not set.
 */
void pfree(void* h_ptr, const DeviceContext* deviceContext)
{
    GMX_RELEASE_ASSERT(deviceContext || g_threadDefaultContext.has_value(),
                       "Calling pfree without a context");
    if (h_ptr)
    {
        if (!deviceContext)
        {
            GMX_RELEASE_ASSERT(g_threadAllocCount > 0,
                               "Calling pfree when we don't expect to have any allocations in the "
                               "default context");
            g_threadAllocCount--;
        }
        const sycl::context& ctx = deviceContext ? deviceContext->context() : *g_threadDefaultContext;
        GMX_ASSERT(sycl::get_pointer_type(h_ptr, ctx) == sycl::usm::alloc::host,
                   "Calling pfree on non-pinned pointer (or with wrong context)");
        sycl::free(h_ptr, ctx);
        if (g_threadDelayContextClearing && g_threadAllocCount == 0)
        {
            g_threadDefaultContext       = std::nullopt;
            g_threadDelayContextClearing = false;
        }
    }
}

void pmallocSetDefaultDeviceContext(const DeviceContext* deviceContext)
{
    GMX_RELEASE_ASSERT(g_threadAllocCount == 0,
                       "Changing default context when there are allocations in it");
    GMX_RELEASE_ASSERT(!g_threadDelayContextClearing,
                       "Can't set the new context while still trying to clear the old one");
    g_threadDefaultContext.emplace(deviceContext->context());
}

void pmallocClearDefaultDeviceContext()
{
    const bool inStackUnwinding             = std::uncaught_exceptions() > 0;
    const bool havePinnedMemoryInTheContext = g_threadAllocCount > 0;
    if (inStackUnwinding && havePinnedMemoryInTheContext)
    {
        /* We are now handling an exception, meaning this function was likely called from a destructor.
         * In Mdrunner::mdrunner(), and probably in other places, we care about the order in which the destructors
         * are called, yet we have no way of enforcing it if an exception is thrown. To avoid leaking memory,
         * we keep the default context, but set this flag to clear the context as soon as all allocations in it are freed.
         */
        g_threadDelayContextClearing = true;
    }
    else
    {
        GMX_RELEASE_ASSERT(g_threadAllocCount == 0,
                           "Clearing default context when there are allocations in it");
        g_threadDefaultContext = std::nullopt;
    }
}
