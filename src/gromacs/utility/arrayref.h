/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * \brief
 * Declares gmx::ConstArrayRef.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_ARRAYREF_H
#define GMX_UTILITY_ARRAYREF_H

#include <cstddef>

#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>

#include "gmxassert.h"

namespace gmx
{

/*! \brief
 * STL container non-mutable interface for a C array (or part of a std::vector).
 *
 * \tparam T  Value type of elements.
 *
 * This class provides an interface similar to \c std::vector<T>, with the
 * following main differences:
 *  - This class does not have its own storage.  Instead, it references an
 *    existing array of values (either a C-style array or part of an existing
 *    std::vector<T>).
 *  - Only const methods are provided to access the stored values.
 *    It is not possible to alter the referenced array.
 *  - Copying objects of this type is cheap, and the copies behave identically
 *    to the original object: the copy references the same set of values.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
template <typename T>
class ConstArrayRef
{
    public:
        //! Type of values stored in the container.
        typedef T         value_type;
        //! Type for representing size of the container.
        typedef size_t    size_type;
        //! Type for representing difference between two container indices.
        typedef ptrdiff_t difference_type;
        //! Const reference to a container element.
        typedef const T  &const_reference;
        //! Const pointer to a container element.
        typedef const T  *const_pointer;
        //! Const iterator type for the container.
        typedef const T  *const_iterator;
        //! Equal to \a const_reference since changes are not allowed.
        typedef const_reference reference;
        //! Equal to \a const_pointer since changes are not allowed.
        typedef const_pointer   pointer;
        //! Equal to \a const_iterator since changes are not allowed.
        typedef const_iterator  iterator;
        //! Standard reverse iterator.
        typedef std::reverse_iterator<iterator>       reverse_iterator;
        //! Standard reverse iterator.
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

        /*! \brief
         * Constructs an empty reference.
         */
        ConstArrayRef() : begin_(NULL), end_(NULL) {}
        /*! \brief
         * Constructs a reference to a particular range.
         *
         * \param[in] begin  Pointer to the beginning of a range.
         * \param[in] end    Pointer to the end of a range.
         *
         * Passed pointers must remain valid for the lifetime of this object.
         */
        ConstArrayRef(const_pointer begin, const_pointer end)
            : begin_(begin), end_(end)
        {
            GMX_ASSERT(end >= begin, "Invalid range");
        }
        /*! \brief
         * Constructs a reference to a particular rangein a std::vector.
         *
         * \param[in] begin  Pointer to the beginning of a range.
         * \param[in] end    Pointer to the end of a range.
         *
         * The referenced vector must remain valid and not be reallocated for
         * the lifetime of this object.
         */
        ConstArrayRef(typename std::vector<T>::const_iterator begin,
                      typename std::vector<T>::const_iterator end)
            : begin_((begin != end) ? &*begin : NULL),
              end_(begin_+(end-begin))
        {
            GMX_ASSERT(end >= begin, "Invalid range");
        }
        /*! \brief
         * Constructs a reference to an array.
         *
         * \param[in] begin  Pointer to the beginning of the array.
         *      May be NULL if \p size is zero.
         * \param[in] size   Number of elements in the array.
         *
         * Passed pointer must remain valid for the lifetime of this object.
         */
        ConstArrayRef(const_pointer begin, size_type size)
            : begin_(begin), end_(begin + size)
        {
        }

        //! Returns an interator to the beginning of the container.
        const_iterator begin() const { return begin_; }
        //! Returns an interator to the end of the container.
        const_iterator end() const { return end_; }
        //! Returns an interator to the reverse beginning of the container.
        const_iterator rbegin() const { return reverse_iterator(end()); }
        //! Returns an interator to the reverse end of the container.
        const_iterator rend() const { return reverse_iterator(begin()); }

        //! Returns the size of the container.
        size_type size() const { return end_ - begin_; }
        //! Identical to size().
        size_type capacity() const { return end_ - begin_; }
        //! Whether the container is empty.
        bool empty() const { return begin_ == end_; }

        //! Access container element.
        const_reference operator[](size_type n) const { return begin_[n]; }
        //! Access container element (throws on out-of-range error).
        const_reference at(size_type n) const
        {
            if (n >= size())
            {
                throw std::out_of_range("Vector index out of range");
            }
            return begin_[n];
        }
        //! Returns the first element in the container.
        const_reference front() const { return *begin_; }
        //! Returns the last element in the container.
        const_reference back() const { return *(end_ - 1); }

        //! Returns a raw pointer to the contents of the array.
        const_pointer data() const { return begin_; }

        /*! \brief
         * Swaps referenced memory with the other object.
         *
         * The actual memory areas are not modified, only the references are
         * swapped.
         */
        void swap(ConstArrayRef<T> &other)
        {
            std::swap(begin_, other.begin_);
            std::swap(end_, other.end_);
        }

    private:
        const_pointer           begin_;
        const_pointer           end_;
};

/*! \brief
 * Simple swap method for ConstArrayRef objects.
 *
 * \see ConstArrayRef::swap()
 */
template <typename T>
void swap(ConstArrayRef<T> &a, ConstArrayRef<T> &b)
{
    a.swap(b);
}

} // namespace gmx

#endif
