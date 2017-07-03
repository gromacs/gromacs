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
 * \brief Declares PageLockedMemory for appropriately handling locking
 * previously allocated page-aligned memory to a physical memory page
 * for use in GPU transfers.
 *
 * For CUDA, handles calling cudaHostRegister and
 * cudaHostUnregister in RAII style.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 */

#ifndef GMX_GPU_UTILS_PAGELOCKEDMEMORY_H
#define GMX_GPU_UTILS_PAGELOCKEDMEMORY_H

#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*! \libinternal
 * \brief Manage "ownership" of just the page locking of memory whose
 * storage allocation and lifetime is managed by another component.
 *
 * The caller must ensure that the host-side memory is aligned to a
 * page boundary (e.g. uses PageAlignedAllocationPolicy), and that the
 * lifetime of the host-side memory exceeds that of this object.
 *
 * Only one constructor is implemented, and there are no semantics for
 * copying, or for taking or releasing the lock to another component
 * because this is not useful. However the lock may be moved or
 * swapped.
 *
 * Note that within GROMACS, only the CUDA path currently implements
 * such page-locking; the OpenCL path has a null implementation (and
 * thus all OpenCL transfers between device and host are subject to
 * buffer-copy overhead).
 *
 * The constructor is templated so that no cast is needed at the point
 * of call. Extern template declarations should be made in source
 * files that need particular instantiations, which are matched with
 * explicit template instantiations in the source file, so that we
 * only compile each constructor once.
 *
 * \todo When implementing support for OpenCL, manual use of mmap()
 * will probably be required, and it is unclear if such an
 * implementation will ever be properly portable.
 */
class PageLockedMemory
{
    private:
        //! Pointer to host-side memory. Once locked, the type is no longer needed.
        const void *memory_;
    public:
        //! No default constructor
        PageLockedMemory() = delete;
        /*! \brief Constructor that locks non-empty host memory in \c arrayRef to its physical page.
         *
         * The caller must ensure this memory has not been locked previously.
         *
         * \tparam    T         Type of the objects in the memory region to lock.
         * \param[in] arrayRef  View of the memory region to lock.
         *
         * \throws  InternalError        When page-locking fails.
         *          NotImplementedError  When called for a non-CUDA build.
         */
        template <typename T>
        explicit PageLockedMemory(ConstArrayRef<T> arrayRef);
        //! Destructor. Does not throw.
        ~PageLockedMemory();
        //! No copy constructor.
        PageLockedMemory(const PageLockedMemory &) = delete;
        //! No copy assignment.
        PageLockedMemory &operator=(const PageLockedMemory &) = delete;
        //! Default move constructor.
        PageLockedMemory(PageLockedMemory &&origin) = default;
        //! Default move assignment.
        PageLockedMemory &operator=(PageLockedMemory &&other) = default;
        /*! \brief Swaps referenced memory with the other object.
         *
         * The actual memory areas are not modified, only the references are
         * swapped.
         *
         * Does not throw.
         */
        void swap(PageLockedMemory &other);
        //! Getter, only intended for use in testing.
        const void *memory() const;
};

} // namespace gmx

#endif
