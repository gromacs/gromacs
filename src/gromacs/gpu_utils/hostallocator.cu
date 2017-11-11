/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief Implements gmx::HostAllocationPolicy for allocating memory
 * suitable for GPU transfers on CUDA.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "hostallocator.h"

#include <cstddef>

#include <memory>

#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

//! Private implementation class.
class HostAllocationPolicy::Impl
{
    public:
        //! Pointer to actual storage, nullptr if no allocation exists.
        void       *pointer_ = nullptr;
        //! Number of bytes in the allocation.
        std::size_t numBytes_ = 0;
        //! Whether this object is in mode where pages will be pinned.
        bool        usingPinningMode_ = false;
        //! The pointer to any storage that has been pinned, or nullptr if none has been pinned.
        void       *pinnedPointer_ = nullptr;
};

HostAllocationPolicy::HostAllocationPolicy() : impl_(std::make_shared<Impl>())
{
}

void *HostAllocationPolicy::malloc(std::size_t bytes) const
{
    // A container could have a pinned allocation that is being
    // extended, in which case we must un-pin while we still know the
    // old pinned vector, and which also ensures we don't pin two
    // buffers at the same time. If there's no allocation, or it isn't
    // pinned then attempting to unpin it is OK, too.
    unpin();
    impl_->pointer_ = PageAlignedAllocationPolicy::malloc(bytes);
    if (impl_->pointer_ != nullptr)
    {
        impl_->numBytes_ = bytes;
    }
    pin();
    return impl_->pointer_;
}

void HostAllocationPolicy::free(void *buffer) const
{
    unpin();
    PageAlignedAllocationPolicy::free(buffer);
    impl_->pointer_  = nullptr;
    impl_->numBytes_ = 0;
}

void HostAllocationPolicy::activatePinningMode()
{
    impl_->usingPinningMode_ = true;
}

void HostAllocationPolicy::deactivatePinningMode()
{
    impl_->usingPinningMode_ = false;
}

//! Is \c ptr aligned on a boundary that is a multiple of \c bytes.
static inline bool isAligned(const void *ptr, size_t bytes)
{
    return (reinterpret_cast<intptr_t>(ptr) % bytes) == 0;
}

namespace
{

//! Ensure our error handling works as well as possible.
static void ensureNoPendingCudaError(const char *descriptiveErrorMessage)
{
    // Ensure there is no pending error that would otherwise affect
    // the behaviour of future error handling.
    cudaError_t stat = cudaGetLastError();
    // If we would find an error in a release build, we do not know
    // what is appropriate to do about it, so assert only for debug
    // builds.
    GMX_ASSERT(stat == cudaSuccess,
               formatString("%s An unhandled error from a previous CUDA operation was detected.", descriptiveErrorMessage).c_str());
}

}   // namespace

void HostAllocationPolicy::pin() const
{
    if (!impl_->usingPinningMode_ || impl_->pointer_ == nullptr)
    {
        // Do nothing if we're not in pinning mode, or the allocation
        // is empty.
        return;
    }
    const char *descriptiveErrorMessage = "Could not register the host memory for page locking for GPU transfers.";

    GMX_ASSERT(impl_->pinnedPointer_ == nullptr,
               formatString("%s Expected no previous pinned allocation to be still pinned", descriptiveErrorMessage).c_str());
    GMX_ASSERT(isAligned(impl_->pointer_, PageAlignedAllocationPolicy::alignment()),
               formatString("%s Host memory needs to be page aligned.", descriptiveErrorMessage).c_str());

    ensureNoPendingCudaError(descriptiveErrorMessage);
    cudaError_t stat = cudaHostRegister(impl_->pointer_, impl_->numBytes_, cudaHostRegisterDefault);

    // These errors can only arise from a coding error somewhere.
    GMX_RELEASE_ASSERT(stat != cudaErrorInvalidValue && stat != cudaErrorNotSupported && stat != cudaErrorHostMemoryAlreadyRegistered,
                       formatString("%s %s: %s", cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());
    // Remaining errors come from correct use of the CUDA API, but it
    // could not do what we asked of it. That should only be
    // cudaErrorMemoryAllocation.
    if (stat != cudaSuccess)
    {
        // We always handle the error, but if it's a type we didn't
        // expect (e.g. because CUDA changes the set of errors it
        // returns) then we should get an assertion in Debug mode so
        // we know to fix our expectations.
        GMX_ASSERT(stat != cudaErrorMemoryAllocation,
                   formatString("%s %s: %s which was an unexpected error", cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());
        GMX_THROW(InternalError(formatString("%s %s: %s", descriptiveErrorMessage, cudaGetErrorName(stat), cudaGetErrorString(stat))));
    }
    impl_->pinnedPointer_ = impl_->pointer_;
}

void HostAllocationPolicy::unpin() const
{
    if (impl_->pinnedPointer_ == nullptr)
    {
        return;
    }

    // If the caller deactivated pinning mode, we still want to be
    // able to unpin if the allocation is still pinned.
    const char *descriptiveErrorMessage = "Could not unregister pinned host memory used for GPU transfers.";

    GMX_ASSERT(impl_->pointer_ != nullptr, formatString("%s pointer should not be nullptr when pinned.", descriptiveErrorMessage).c_str());

    ensureNoPendingCudaError(descriptiveErrorMessage);
    cudaError_t stat = cudaHostUnregister(impl_->pinnedPointer_);
    // These errors can only arise from a coding error somewhere.
    GMX_RELEASE_ASSERT(stat != cudaErrorInvalidValue && stat != cudaErrorHostMemoryNotRegistered,
                       formatString("%s %s: %s", cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());
    // If there's an error whose type we didn't expect (e.g. because a
    // future CUDA changes the set of errors it returns) then we
    // should assert, because our code is wrong.
    //
    // The approach differs from that in pin() because we might
    // unpin() from a destructor, in which case any attempt to throw
    // an uncaught exception would anyway terminate the program. A
    // release assertion is a better behaviour than that.
    GMX_RELEASE_ASSERT(stat == cudaSuccess,
                       formatString("%s %s: %s which was an unexpected error", cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());

    impl_->pinnedPointer_ = nullptr;
}

} // namespace gmx
