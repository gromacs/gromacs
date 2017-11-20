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
/*! \file
 * \brief Declares gmx::HostAllocationPolicy, gmx::HostAllocator, and
 * gmx::HostVector, which are used to make/be standard library
 * containers that can allocate memory suitable for transfers.
 * Currently the only supported transfers using pinned memory are
 * to CUDA GPUs, but other possibilities exist in future.
 *
 * \todo This should not be in the public API, but it needs to be
 * for the moment because state.h is in that API.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inpublicapi
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

/*! \brief Helper enum for pinning policy of the allocation of
 * HostAllocationPolicy.
 *
 * For an efficient non-blocking transfer (e.g. to a GPU), the memory
 * pages for a buffer need to be pinned to a physical page. Aligning
 * such buffers to a physical page should miminize the number of pages
 * that need to be pinned. However, some buffers that may be used for
 * such transfers may also be used in either GROMACS builds or run
 * paths that cannot use such a device, so the policy can be
 * configured so that the resource consumption is no higher than
 * required for correct, efficient operation in all cases. */
enum class PinningPolicy : int
{
    CannotBePinned,     // Memory is not known to be suitable for pinning.
    CanBePinned,        // Memory is suitable for efficient pinning, e.g. because it is
                        // allocated to be page aligned, and will be pinned when non-empty.
};

//! Forward declaration of host allocation policy class.
class HostAllocationPolicy;

/*! \brief Memory allocator that uses HostAllocationPolicy.
 *
 *  \tparam T          Type of objects to allocate
 *
 * This convenience partial specialization can be used for the
 * optional allocator template parameter in standard library
 * containers whose memory may be used for e.g. GPU transfers. The
 * memory will always be allocated according to the behavior of
 * HostAllocationPolicy.
 */
template <class T>
using HostAllocator = Allocator<T, HostAllocationPolicy>;

//! Convenience alias for std::vector that uses HostAllocator.
template <class T>
using HostVector = std::vector<T, HostAllocator<T> >;

/*! \libinternal
 * \brief Policy class for configuring gmx::Allocator, to manage
 * allocations of memory that may be needed for e.g. GPU transfers.
 *
 * This allocator has state, so is most useful in cases where it is
 * not known at compile time whether the allocated memory will be
 * transferred to some device. It will increase the size of containers
 * that use it. If the GROMACS build is configured with CUDA support,
 * then memory will be allocated with PageAlignedAllocator, and that
 * page pinned to physical memory if the pinning mode has been
 * activated. If pinning mode is deactivated, or the GROMACS build
 * does not support CUDA, then the memory will be allocated with
 * AlignedAllocator. The pin() and unpin() methods work with the CUDA
 * build, and silently do nothing otherwise. In future, we may modify
 * or generalize this to work differently in other cases.
 *
 * The intended use is to configure gmx::Allocator with this class as
 * its policy class, and then to use e.g.
 * std::vector::get_allocator().getPolicy() to control whether the
 * allocation policy should activate its pinning mode. The policy
 * object can also be used to explicitly pin() and unpin() the buffer
 * when it is using PinningPolicy::CanBePinned. The policy object is
 * returned by value (as required by the C++ standard for
 * get_allocator(), which copies a std::shared_ptr, so the policy
 * object should be retrieved sparingly, e.g. only upon resize of the
 * allocation. (Normal operation of the vector, e.g. during resize,
 * incurs only the cost of the pointer indirection needed to consult
 * the current state of the allocation policy.)
 *
 * \todo As a minor optimization, consider also having a stateless
 * version of this policy, which might be slightly faster or more
 * convenient to use in the cases where it is known at compile time
 * that the allocation will be used to transfer to a GPU.
 */
class HostAllocationPolicy
{
    public:
        //! Default constructor.
        HostAllocationPolicy();
        /*! \brief Return the alignment size currently used by the active pinning policy. */
        std::size_t alignment();
        /*! \brief Allocate and perhaps pin page-aligned memory suitable for
         * e.g. GPU transfers.
         *
         * Before attempting to allocate, unpin() is called. After a
         * successful allocation, pin() is called. (Whether these do
         * things depends on the PinningPolicy that is in effect.)
         *
         *  \param bytes Amount of memory (bytes) to allocate. It is valid to ask for
         *               0 bytes, which will return a non-null pointer that is properly
         *               aligned and padded (but that you should not use).
         *
         *  \return Valid pointer if the allocation+optional pinning worked, otherwise nullptr.
         *
         *  \note Memory allocated with this routine must be released
         *        with gmx::HostAllocationPolicy::free(), and
         *        absolutely not the system free().
         *
         * Does not throw.
         */
        void *malloc(std::size_t bytes) const noexcept;
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
        void free(void *buffer) const noexcept;
        /*! \brief Pin the allocation to physical memory, if appropriate.
         *
         * If the allocation policy is not in pinning mode, or the
         * allocation is empty, ot the allocation is already pinned,
         * then do nothing.
         *
         * Does not throw.
         */
        void pin() const noexcept;
        /*! \brief Unpin the allocation, if appropriate.
         *
         * Regardless of the allocation policy, unpin the memory if
         * previously pinned, otherwise do nothing.
         *
         * Does not throw.
         */
        void unpin() const noexcept;
        /*! \brief Return the current pinning policy (which is semi-independent
         * of whether the buffer is actually pinned).
         *
         * Does not throw.
         */
        PinningPolicy pinningPolicy() const;
        //! Specify an allocator trait so that the stateful allocator should propagate.
        using propagate_on_container_copy_assignment = std::true_type;
        //! Specify an allocator trait so that the stateful allocator should propagate.
        using propagate_on_container_move_assignment = std::true_type;
        //! Specify an allocator trait so that the stateful allocator should propagate.
        using propagate_on_container_swap = std::true_type;
    private:
        /*! \brief Set the current pinning policy.
         *
         * Does not pin any current buffer. Use changePinningPolicy to
         * orchestrate the necessary unpin, allocate, copy, pin for
         * effectively changing the pinning policy of a HostVector.
         *
         * Does not throw.
         */
        // cppcheck-suppress unusedPrivateFunction
        void setPinningPolicy(PinningPolicy pinningPolicy);
        /*! \brief Declare as a friend function the only supported way
         * to change the pinning policy.
         *
         * When the pinning policy changes, we want the state of the
         * allocation to match the new policy. However, that requires
         * a copy and swap of the buffers, which can only take place
         * at the level of the container. So we wrap the required
         * operations in a helper friend function.
         *
         * Of course, if there is no allocation because the vector is
         * empty, then nothing will change. */
        template <class T> friend
        void changePinningPolicy(HostVector<T> *v, PinningPolicy pinningPolicy);
        //! Private implementation class.
        class Impl;
        /*! \brief State of the allocator.
         *
         * This could change through move- or copy-assignment of one
         * policy to another, so isn't const. */
        std::shared_ptr<Impl> impl_;
};

/*! \brief Helper function for changing the pinning policy of a HostVector.
 *
 * If the vector has contents, then a full reallocation and buffer
 * copy are needed if the policy change requires tighter restrictions,
 * and desirable even if the policy change requires looser
 * restrictions. That cost is OK, because GROMACS will do this
 * operation very rarely (e.g. when auto-tuning and deciding to switch
 * whether a task will run on a GPU, or not). */
template <class T>
void changePinningPolicy(HostVector<T> *v, PinningPolicy pinningPolicy)
{
    // Do we have anything to do?
    HostAllocationPolicy vAllocationPolicy = v->get_allocator().getPolicy();
    if (pinningPolicy == vAllocationPolicy.pinningPolicy())
    {
        return;
    }
    // Make sure we never have two allocated buffers that are both pinned.
    vAllocationPolicy.unpin();

    // Construct a new vector that has the requested
    // allocation+pinning policy, to swap into *v. If *v is empty,
    // then no real work is done.
    HostAllocator<T> newAllocator;
    newAllocator.getPolicy().setPinningPolicy(pinningPolicy);
    HostVector<T>    newV(v->begin(), v->end(), newAllocator);
    // Replace the contents of *v, including the stateful allocator.
    v->swap(newV);
    // The destructor of newV cleans up the memory formerly managed by *v.
}

}      // namespace gmx

#endif
