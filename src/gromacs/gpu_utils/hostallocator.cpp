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
 * suitable for e.g. GPU transfers on CUDA.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "hostallocator.h"

#include "config.h"

#include <cstddef>

#include <memory>

#include "gromacs/gpu_utils/pinning.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

//! Private implementation class.
class HostAllocationPolicy::Impl
{
    public:
        /*! \brief Pointer to the last unfreed allocation, or nullptr
         * if no allocation exists.
         *
         * Note that during e.g. std::vector.resize() a call to its
         * allocator's allocate() function precedes the call to its
         * allocator's deallocate() function for freeing the old
         * buffer after the data has been copied from it. So in
         * general, pointer_ will not match the argument received by
         * free(). */
        void         *pointer_ = nullptr;
        //! Number of bytes in the last unfreed allocation.
        std::size_t   numBytes_ = 0;
        //! The pointer to any storage that has been pinned, or nullptr if none has been pinned.
        void         *pinnedPointer_ = nullptr;
        //! Whether this object is in mode where new allocations will be pinned by default.
        PinningPolicy pinningPolicy_ = PinningPolicy::CannotBePinned;
};

HostAllocationPolicy::HostAllocationPolicy() : impl_(std::make_shared<Impl>())
{
}

std::size_t HostAllocationPolicy::alignment()
{
    return (impl_->pinningPolicy_ == PinningPolicy::CanBePinned ?
            PageAlignedAllocationPolicy::alignment() :
            AlignedAllocationPolicy::alignment());
}
void *HostAllocationPolicy::malloc(std::size_t bytes) const noexcept
{
    // A container could have a pinned allocation that is being
    // extended, in which case we must un-pin while we still know the
    // old pinned vector, and which also ensures we don't pin two
    // buffers at the same time. If there's no allocation, or it isn't
    // pinned, then attempting to unpin it is OK, too.
    unpin();
    impl_->pointer_ = (impl_->pinningPolicy_ == PinningPolicy::CanBePinned ?
                       PageAlignedAllocationPolicy::malloc(bytes) :
                       AlignedAllocationPolicy::malloc(bytes));

    if (impl_->pointer_ != nullptr)
    {
        impl_->numBytes_ = bytes;
    }
    pin();
    return impl_->pointer_;
}

void HostAllocationPolicy::free(void *buffer) const noexcept
{
    unpin();
    if (buffer == nullptr)
    {
        // Nothing to do
        return;
    }
    if (impl_->pinningPolicy_ == PinningPolicy::CanBePinned)
    {
        PageAlignedAllocationPolicy::free(buffer);
    }
    else
    {
        AlignedAllocationPolicy::free(buffer);
    }
    impl_->pointer_  = nullptr;
    impl_->numBytes_ = 0;
}

PinningPolicy HostAllocationPolicy::pinningPolicy() const
{
    return impl_->pinningPolicy_;
}

void HostAllocationPolicy::setPinningPolicy(PinningPolicy pinningPolicy)
{
    if (GMX_GPU != GMX_GPU_CUDA)
    {
        GMX_RELEASE_ASSERT(pinningPolicy == PinningPolicy::CannotBePinned,
                           "A suitable build of GROMACS (e.g. with CUDA) is required for a "
                           "HostAllocationPolicy to be set to a mode that produces pinning.");
    }
    impl_->pinningPolicy_ = pinningPolicy;
}

void HostAllocationPolicy::pin() const noexcept
{
    if (impl_->pinningPolicy_ == PinningPolicy::CannotBePinned ||
        impl_->pointer_ == nullptr ||
        impl_->pinnedPointer_ != nullptr)
    {
        // Do nothing if we're not in pinning mode, or the allocation
        // is empty, or it is already pinned.
        return;
    }
#if GMX_GPU == GMX_GPU_CUDA
    pinBuffer(impl_->pointer_, impl_->numBytes_);
    impl_->pinnedPointer_ = impl_->pointer_;
#else
    const char *errorMessage = "Could not register the host memory for pinning.";

    GMX_RELEASE_ASSERT(impl_->pinningPolicy_ == PinningPolicy::CannotBePinned,
                       formatString("%s This build configuration must only have pinning policy "
                                    "that leads to no pinning.", errorMessage).c_str());
#endif
}

void HostAllocationPolicy::unpin() const noexcept
{
    if (impl_->pinnedPointer_ == nullptr)
    {
        return;
    }

#if GMX_GPU == GMX_GPU_CUDA
    // Note that if the caller deactivated pinning mode, we still want
    // to be able to unpin if the allocation is still pinned.

    unpinBuffer(impl_->pointer_);
    impl_->pinnedPointer_ = nullptr;
#else
    GMX_RELEASE_ASSERT(impl_->pinnedPointer_ == nullptr,
                       "Since the build configuration does not support pinning, then "
                       "the pinned pointer must be nullptr.");
#endif
}

} // namespace gmx
