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
/*! \libinternal \file
 * \brief Declares gmx::HostAllocationPolicy, gmx::HostAllocator, and
 * gmx::HostVector, which are used to make/be standard library
 * containers that can allocate memory suitable for GPU transfers.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_HOSTALLOCATOR_H
#define GMX_GPU_UTILS_HOSTALLOCATOR_H

#include <cstddef>

#include <memory>
#include <vector>

#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

/*! \libinternal
 * \brief Policy class for configuring gmx::Allocator, to manage
 * allocations of memory that is suitable for GPU transfers.
 *
 * This allocator has state, so is most useful in cases where it is
 * not known at compile time whether the allocated memory will be
 * transferred to a GPU. It will increase the size of containers that
 * use it. Memory allocated will always be aligned by the GPU
 * framework, or by AlignedAllocationPolicy.
 *
 * The intended use is to configure gmx::Allocator as its policy
 * class, and then to use std::vector::get_allocator().getPolicy() to
 * control whether the allocation policy should activate its pinning
 * mode. The policy object can also be used to explicitly pin() and
 * unpin() the buffer. The policy object is returned by value (as
 * required by the C++ standard for get_allocator(), which copies a
 * std::shared_ptr, so the policy object should be retrieved
 * sparingly, e.g. upon resize of the allocation. (Normal operation of
 * the vector incurs only the cost of the pointer indirection needed
 * to consult the current state of the allocation policy.)
 *
 * \todo Consider also having a stateless version of this policy,
 * which might be slightly faster or more convenient to use in the
 * cases where it is known at compile time that the allocation will be
 * used to transfer to a GPU.
 */
class HostAllocationPolicy
{
    public:
        //! Default constructor.
        HostAllocationPolicy();
        /*! \brief Allocate and perhaps pin page-aligned memory suitable for GPU transfers.
         *
         * Before attempting to allocate, unpin() is called. After a successful allocation,
         * pin() is called.
         *
         *  \param bytes Amount of memory (bytes) to allocate. It is valid to ask for
         *               0 bytes, which will return a non-null pointer that is properly
         *               aligned and padded (but that you should not use).
         *
         *  \return Valid pointer if the allocation worked, otherwise nullptr.
         *
         *  \note Memory allocated with this routine must be released
         *        with gmx::HostAllocationPolicy::free(), and
         *        absolutely not the system free().
         *
         * Does not throw.
         */
        void *malloc(std::size_t bytes) const;
        /*! \brief Free the memory, after unpinning (if appropriate).
         *
         *  \param buffer  Memory pointer previously returned from gmx::HostAllocationPolicy::malloc()
         *
         *  \note This routine should only be called with pointers
         *        obtained from gmx:HostAllocationPolicy::malloc(),
         *        and absolutely not any pointers obtained the system
         *        malloc().
         *
         * Does not throw.
         */
        void free(void *buffer) const;
        /*! \brief Pin the allocation to physical memory, if appropriate.
         *
         * If the allocation policy is not in pinning mode, or the
         * allocation is empy, ot the allocation is already pinned,
         * then do nothing.
         *
         * \throws  InternalError  Upon any unexpected error from an underlying API.
         */
        void pin() const;
        /*! \brief Unpin the allocation, if appropriate.
         *
         * Regardless of the allocation policy, unpin the memory if
         * previously pinned, otherwise do nothing.
         *
         * Does not throw.
         */
        void unpin() const;
        /*! \brief Activate pinning mode.
         *
         * Does not pin the current buffer.
         *
         * Does not throw.
         */
        void activatePinningMode();
        /*! \brief Deactivate pinning mode.
         *
         * Does not unpin an already pinned buffer.
         *
         * Does not throw.
         */
        void deactivatePinningMode();
    private:
        //! Private implementation class.
        class Impl;
        /*! \brief State of the allocator.
         *
         * This could change through assignment of one policy to
         * another, such as returned so isn't const. */
        std::shared_ptr<Impl> impl_;
};

/*! \brief Memory allocator for host-side memory for GPU transfers.
 *
 *  \tparam T          Type of objects to allocate
 *
 * This convenience partial specialization can be used for the
 * optional allocator template parameter in standard library
 * containers whose memory will be used for GPU transfers. The memory
 * will always be allocated according to the behavior of
 * HostAllocationPolicy.
 */
template <class T>
using HostAllocator = Allocator<T, HostAllocationPolicy>;

//! Convenience alias for host-side vector suited for efficient GPU transfers.
template <class T>
using HostVector = std::vector<T, HostAllocator<T> >;

}      // namespace gmx

#endif
