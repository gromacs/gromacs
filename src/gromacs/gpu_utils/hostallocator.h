/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
/*! \file
 * \brief Declares gmx::HostAllocationPolicy, gmx::HostAllocator,
 * gmx::HostVector and gmx::PaddedHostVector, which are used to make/be
 * standard library containers that can allocate memory suitable for transfers.
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
#include <type_traits>
#include <vector>

#include "gromacs/math/paddedvector.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/allocator.h"
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
    CannotBePinned,    // Memory is not known to be suitable for pinning.
    PinnedIfSupported, // Memory is suitable for efficient pinning, e.g. because it is
                       // allocated to be page aligned, and will be pinned when supported.
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
template<class T>
using HostAllocator = Allocator<T, HostAllocationPolicy>;

//! Convenience alias for std::vector that uses HostAllocator.
template<class T>
using HostVector = std::vector<T, HostAllocator<T>>;

//! Convenience alias for PaddedVector that uses HostAllocator.
template<class T>
using PaddedHostVector = PaddedVector<T, HostAllocator<T>>;

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
 * when it is using PinningPolicy::PinnedIfSupported. The policy object is
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
    /*! \brief Constructor
     *
     * \param[in] policy
     *                Whether to pin the allocation
     * \param[in] propagateDuringContainerCopyConstruction
     *                Default is chosen to be consistent with copy assignment
     */
    HostAllocationPolicy(PinningPolicy policy = PinningPolicy::CannotBePinned,
                         bool          propagateDuringContainerCopyConstruction = false);
    /*! \brief Return the alignment size currently used by the active pinning policy. */
    std::size_t alignment() const noexcept;
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
    void* malloc(std::size_t bytes) const noexcept;
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
    void free(void* buffer) const noexcept;
    /*! \brief Return the active pinning policy.
     *
     * Does not throw.
     */
    PinningPolicy pinningPolicy() const { return pinningPolicy_; }
    /*! \brief Do not propagate the allocator for copy assignment
     *
     * We choose that the allocator is a property of the container,
     * and should be changed explicitly as required with e.g. \c
     * changePinningPolicy if the usage dictates that the copy adopts
     * the policy of the original container or a specific policy.
     */
    using propagate_on_container_copy_assignment = std::false_type;
    //! Propagate the allocator for move assignment
    using propagate_on_container_move_assignment = std::true_type;
    //! Propagate the allocator during swap
    using propagate_on_container_swap = std::true_type;
    //! \brief Return the policy the container should use for copy construction
    // NOLINTNEXTLINE readability-convert-member-functions-to-static
    HostAllocationPolicy select_on_container_copy_construction() const
    {
        if (propagateDuringContainerCopyConstruction_)
        {
            return *this;
        }
        else
        {
            return {};
        }
    }

private:
    //! Pinning policy
    PinningPolicy pinningPolicy_;
    //! Whether to propagate the allocator during copy construction by a container.
    bool propagateDuringContainerCopyConstruction_;
};

/*! \brief Return true if two allocators are identical
 *
 * True if pinning policy is the same.
 */
template<class T1, class T2>
bool operator==(const Allocator<T1, HostAllocationPolicy>& a, const Allocator<T2, HostAllocationPolicy>& b)
{
    return a.pinningPolicy() == b.pinningPolicy();
}

/*! \brief Helper function for changing the pinning policy of a pinnable vector.
 *
 * If the vector has contents, then a full reallocation and buffer
 * copy are needed if the policy change requires tighter restrictions,
 * and desirable even if the policy change requires looser
 * restrictions. That cost is OK, because GROMACS will do this
 * operation very rarely (e.g. when auto-tuning and deciding to switch
 * whether a task will run on a GPU, or not). */
template<typename PinnableVector>
void changePinningPolicy(PinnableVector* v, PinningPolicy pinningPolicy)
{
    // Force reallocation by element-wise move (because policy is
    // different container is forced to realloc). Does nothing if
    // policy is the same.
    *v = PinnableVector(std::move(*v), { pinningPolicy });
}

//! Convenience type for vector with aligned memory
template<typename T>
using AlignedVector = std::vector<T, AlignedAllocator<T>>;

} // namespace gmx

#endif
