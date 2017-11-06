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
 * \brief Declares gmx::HostAllocationPolicy and gmx::HostAllocator,
 * which are used to make standard library containers that can
 * allocate memory suitable for GPU transfers.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_HOSTALLOCATOR_H
#define GMX_GPU_UTILS_HOSTALLOCATOR_H

#include <cstddef>

#include "gromacs/utility/allocator.h"

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
 * \todo Consider also having a stateless version of this policy,
 * which might be slightly faster or more convenient to use in the
 * cases where it is known at compile time that the allocation will be
 * used to transfer to a GPU.
 */
class HostAllocationPolicy
{
    public:
        //! Helper construction enum
        enum class Impl : int
        {
            AllocateAligned  = 0,
            AllocateForGpu   = 1
        };
        //! Constructor.
        explicit HostAllocationPolicy(Impl s = Impl::AllocateAligned);
        /*! \brief Allocate GPU memory
         *
         *  \param bytes Amount of memory (bytes) to allocate. It is
         *               valid to ask for 0 bytes, which will return a
         *               non-null pointer that is properly aligned in
         *               page-locked memory (but that you should not
         *               use). TODO check this.
         *
         * \return Valid pointer if the allocation worked, otherwise nullptr.
         *
         * The memory will always be allocated according to the requirements
         * of the acceleration platform in use (e.g. CUDA).
         *
         *  \note Memory allocated with this routine must be released
         *        with gmx::HostAllocationPolicy::free(), and
         *        absolutely not the system free().
         */
        void *
        malloc(std::size_t bytes) const;
        /*! \brief Free GPU memory
         *
         *  \param buffer  Memory pointer previously returned from gmx::HostAllocationPolicy::malloc()
         *
         *  \note This routine should only be called with pointers
         *        obtained from gmx:HostAllocationPolicy::malloc(),
         *        and absolutely not any pointers obtained the system
         *        malloc().
         */
        void
        free(void *buffer) const;
    private:
        /*! \brief State of the allocator.
         *
         * This could change through assignment of one policy to
         * another, so isn't const. */
        Impl allocateForGpu_;
};

/*! \brief Convenience function
 *
 * The default construction is for non-GPU allocation, and this
 * function makes it less verbose to get allocation intended for use
 * with a GPU. */
HostAllocationPolicy makeHostAllocationPolicyForGpu();

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

}      // namespace gmx

#endif
