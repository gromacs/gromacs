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
 * \brief Declares gmx::HostVector, which are used to make
 * standard-conforming containers that can allocate memory suitable
 * for GPU transfers.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_HOSTVECTOR_H
#define GMX_GPU_UTILS_HOSTVECTOR_H

#include <cstddef>

#include <vector>

#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \libinternal
 * \brief Vector-like class that can be configured to
 * be suitable for efficient GPU transfers.
 */
template <typename T>
class HostVector
{
    private:
        //! The type of the underlying vector.
        using VectorType = std::vector < T, PageAlignedAllocator < T>>;
        //! The type representing an amount of memory.
        using size_type = typename VectorType::size_type;
        //! The type of a value contained.
        using value_type = typename VectorType::size_type;
        //! The type of a reference to a value contained.
        using reference = typename VectorType::reference;
        //! The type of a const reference to a value contained.
        using const_reference = typename VectorType::const_reference;
        //! The type of a pointer to a value contained.
        using pointer = typename VectorType::pointer;
        //! The type of a const pointer to a value contained.
        using const_pointer = typename VectorType::const_pointer;
        //! The type of a iterator to a value contained.
        using iterator = typename VectorType::iterator;
        //! The type of a const iterator to a value contained.
        using const_iterator = typename VectorType::const_iterator;
        //! The type of a reverse iterator to a value contained.
        using reverse_iterator = typename VectorType::reverse_iterator;
        //! The type of a const reverse iterator to a value contained.
        using const_reverse_iterator = typename VectorType::const_reverse_iterator;
        //! Vector of actual storage.
        VectorType v_;
        //! Whether this object is in mode where pages will be locked.
        bool       usingLockingMode_;
        //! Whether the allocation is locked to a page that will suit efficient GPU transfer.
        bool       isLocked_;
    public:
        //! Default constructor.
        HostVector() : v_(), usingLockingMode_(false), isLocked_(false) {};
        // TODO any more constructors we need?
        //! Default destructor.
        ~HostVector();
        //! Getter
        pointer data()
        {
            return v_.data();
        }
        //! Getter
        const_pointer data() const
        {
            return v_.data();
        }
        //! Subscript accessor.
        reference operator[](size_type n) { return v_[n]; };
        //! Subscript const accessor.
        const_reference operator[](size_type n) const { return v_[n]; };
        //!  Returns a read/write iterator to the first element.
        iterator begin() { return v_.begin(); };
        //!  Returns a read-only (const) iterator to the first element.
        const_iterator begin() const { return v_.begin(); };
        //! Returns a read/write iterator that points one past the last element.
        iterator end() { return v_.end(); };
        //!  Returns a read-only (const) iterator that points one past the last element.
        const_iterator end() const { return v_.end(); };
        //!  Returns a read/write reverse iterator to the first element.
        reverse_iterator rbegin() { return v_.rbegin(); };
        //!  Returns a read-only (const) reverse iterator to the first element.
        const_reverse_iterator rbegin() const { return v_.rbegin(); };
        //! Returns a read/write reverse iterator that points one past the last element.
        reverse_iterator rend() { return v_.rend(); };
        //!  Returns a read-only (const) reverse iterator that points one past the last element.
        const_reverse_iterator rend() const { return v_.rend(); };
        //!  Returns the capacity of elements.
        size_type capacity() const { return v_.capacity(); };
        //!  Returns the number of elements.
        size_type size() const { return v_.size(); };
        //!  Returns the largest possible size() without reallocation.
        size_type max_size() const { return v_.max_size(); };
        //! Is the allocation empty?
        bool empty() const { return v_.empty(); };
        //! Set the memory-management mode to use.
        void useLockingMode(bool newMode);
        //! Lock the allocation to physical memory, if it is not already locked.
        void lock();
        //! Unlock the allocation, if it is locked.
        void unlock();
        //! Reserve capacity of \c newCapacity elements.
        void reserve(size_type newCapacity);
        //! Resize the allocation.
        void resize(size_type newSize);
        //! Free the allocation.
        void free();
};

// Doxygen does not understand extern templates.
#if !defined(DOXYGEN)
extern template class HostVector<real>;
#endif
}      // namespace gmx

#endif
