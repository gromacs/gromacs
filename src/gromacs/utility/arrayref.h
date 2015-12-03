/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
 * \brief
 * Declares gmx::ArrayRef and gmx::ConstArrayRef.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
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

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief
 * Tag type to initialize empty array references.
 *
 * This type (together with appropriate constructors in ArrayRef and
 * ConstArrayRef) allows initializing any array reference to an empty value
 * without explicitly specifying its type.  This is convenient when calling
 * a function that takes an array reference, where constructing an empty
 * reference explicitly would otherwise require specifying the full array
 * reference type, including the template parameter.
 */
struct EmptyArrayRef {};

/*! \brief
 * STL-like container for an interface to a C array (or part of a std::vector).
 *
 * \tparam T  Value type of elements.
 *
 * This class provides an interface similar to \c std::vector<T>, with the
 * following main differences:
 *  - This class does not have its own storage.  Instead, it references an
 *    existing array of values (either a C-style array or part of an existing
 *    std::vector<T>).
 *  - It is only possible to modify the values themselves through ArrayRef;
 *    it is not possible to add or remove values.
 *  - Copying objects of this type is cheap, and the copies behave identically
 *    to the original object: the copy references the same set of values.
 *
 * This class is useful for writing wrappers that expose a different view of
 * the internal data stored as a single vector/array.
 *
 * Methods in this class do not throw, except where indicated.
 *
 * Note that due to a Doxygen limitation, the constructor that takes a C array
 * whose size is known at compile time does not appear in the documentation.
 *
 * \todo
 * This class is not complete.  At least, it should be possible to convert an
 * ArrayRef to a ConstArrayRef.  There are likely also methods missing (not
 * required for current usage).
 *
 * \inpublicapi
 * \ingroup module_utility
 */
template <typename T>
class ArrayRef
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
        //! Reference to a container element.
        typedef T        &reference;
        //! Pointer to a container element.
        typedef T        *pointer;
        //! Iterator type for the container.
        typedef T        *iterator;
        //! Standard reverse iterator.
        typedef std::reverse_iterator<iterator>       reverse_iterator;
        //! Standard reverse iterator.
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

        /*! \brief
         * Constructs a reference to a particular range from two pointers.
         *
         * \param[in] begin  Pointer to the beginning of a range.
         * \param[in] end    Pointer to the end of a range.
         *
         * Passed pointers must remain valid for the lifetime of this object.
         */
        static ArrayRef<value_type>
        fromPointers(value_type *begin, value_type *end)
        {
            return ArrayRef<value_type>(begin, end);
        }
        /*! \brief
         * Constructs a reference to an array.
         *
         * \param[in] begin  Pointer to the beginning of the array.
         *                   May be NULL if \p size is zero.
         * \param[in] size   Number of elements in the array.
         *
         * Passed pointer must remain valid for the lifetime of this object.
         */
        static ArrayRef<value_type>
        fromArray(value_type *begin, size_t size)
        {
            return ArrayRef<value_type>(begin, begin+size);
        }
        /*! \brief
         * Constructs a reference to a particular range in a std::vector.
         *
         * \param[in] begin  Iterator to the beginning of a range.
         * \param[in] end    Iterator to the end of a range.
         *
         * The referenced vector must remain valid and not be reallocated for
         * the lifetime of this object.
         */
        static ArrayRef<value_type>
        fromVector(typename std::vector<value_type>::iterator begin,
                   typename std::vector<value_type>::iterator end)
        {
            value_type *p_begin = (begin != end) ? &*begin : NULL;
            value_type *p_end   = p_begin + (end-begin);
            return ArrayRef<value_type>(p_begin, p_end);
        }

        /*! \brief
         * Constructs an empty reference.
         */
        ArrayRef() : begin_(NULL), end_(NULL) {}
        /*! \brief
         * Constructs an empty reference.
         *
         * This is provided for convenience, such that EmptyArrayRef can be
         * used to initialize any ArrayRef, without specifying the template
         * type.  It is not explicit to enable that usage.
         */
        ArrayRef(const EmptyArrayRef &) : begin_(NULL), end_(NULL) {}
        /*! \brief
         * Constructs a reference to a particular range.
         *
         * \param[in] begin  Pointer to the beginning of a range.
         * \param[in] end    Pointer to the end of a range.
         *
         * Passed pointers must remain valid for the lifetime of this object.
         *
         * \note For clarity, use the non-member function arrayRefFromPointers
         * instead.
         */
        ArrayRef(pointer begin, pointer end)
            : begin_(begin), end_(end)
        {
            GMX_ASSERT(end >= begin, "Invalid range");
        }
        /*! \brief
         * Constructs a reference to a whole vector.
         *
         * \param[in] v  Vector to reference.
         *
         * Passed vector must remain valid and not be reallocated for the
         * lifetime of this object.
         *
         * This constructor is not explicit to allow directly passing
         * std::vector to a method that takes ArrayRef.
         */
        ArrayRef(std::vector<T> &v)
            : begin_((!v.empty()) ? &v[0] : NULL),
              end_((!v.empty()) ? &v[0] + v.size() : NULL)
        {
        }
        //! \cond
        // Doxygen 1.8.5 doesn't parse the declaration correctly...
        /*! \brief
         * Constructs a reference to a C array.
         *
         * \param[in] array  C array to reference.
         * \tparam    count  Deduced number of elements in \p array.
         *
         * This constructor can only be used with a real array (not with a
         * pointer).  It constructs a reference to the whole array, without
         * a need to pass the number of elements explicitly.  The compiler
         * must be able to deduce the array size.
         *
         * Passed array must remain valid for the lifetime of this object.
         *
         * This constructor is not explicit to allow directly passing
         * a C array to a function that takes an ArrayRef parameter.
         *
         * xlc on BG/Q compiles wrong code if the C array is a struct
         * field, unless value_type is char or unsigned char. There's
         * no good way to assert on this before C++11 (which that
         * compiler will never support).
         */
        template <size_t count>
        ArrayRef(value_type (&array)[count])
            : begin_(array), end_(array + count)
        {
        }
        //! \endcond

        //! Returns an iterator to the beginning of the container.
        iterator begin() { return begin_; }
        //! Returns an iterator to the beginning of the container.
        const_iterator begin() const { return begin_; }
        //! Returns an iterator to the end of the container.
        iterator end() { return end_; }
        //! Returns an iterator to the end of the container.
        const_iterator end() const { return end_; }
        //! Returns an iterator to the reverse beginning of the container.
        iterator rbegin() { return reverse_iterator(end()); }
        //! Returns an iterator to the reverse beginning of the container.
        const_iterator rbegin() const { return reverse_iterator(end()); }
        //! Returns an iterator to the reverse end of the container.
        iterator rend() { return reverse_iterator(begin()); }
        //! Returns an iterator to the reverse end of the container.
        const_iterator rend() const { return reverse_iterator(begin()); }

        //! Returns the size of the container.
        size_type size() const { return end_ - begin_; }
        //! Identical to size().
        size_type capacity() const { return end_ - begin_; }
        //! Whether the container is empty.
        bool empty() const { return begin_ == end_; }

        //! Access container element.
        reference operator[](size_type n) { return begin_[n]; }
        //! Access container element.
        const_reference operator[](size_type n) const { return begin_[n]; }
        //! Access container element (throws on out-of-range error).
        reference at(size_type n)
        {
            if (n >= size())
            {
                throw std::out_of_range("Vector index out of range");
            }
            return begin_[n];
        }
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
        reference front() { return *begin_; }
        //! Returns the first element in the container.
        const_reference front() const { return *begin_; }
        //! Returns the last element in the container.
        reference back() { return *(end_ - 1); }
        //! Returns the last element in the container.
        const_reference back() const { return *(end_ - 1); }

        //! Returns a raw pointer to the contents of the array.
        pointer data() { return begin_; }
        //! Returns a raw pointer to the contents of the array.
        const_pointer data() const { return begin_; }

        /*! \brief
         * Swaps referenced memory with the other object.
         *
         * The actual memory areas are not modified, only the references are
         * swapped.
         */
        void swap(ArrayRef<T> &other)
        {
            std::swap(begin_, other.begin_);
            std::swap(end_, other.end_);
        }

    private:
        pointer           begin_;
        pointer           end_;
};



/*! \brief
 * STL-like container for non-mutable interface to a C array (or part of a
 * std::vector).
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
 * This class is useful for writing wrappers that expose a different view of
 * the internal data stored as a single vector/array.
 *
 * Methods in this class do not throw, except where indicated.
 *
 * Note that due to a Doxygen limitation, the constructor that takes a C array
 * whose size is known at compile time does not appear in the documentation.
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

        //! \copydoc ArrayRef::fromPointers()
        static ConstArrayRef<value_type>
        fromPointers(const value_type *begin, const value_type *end)
        {
            return ConstArrayRef<value_type>(begin, end);
        }
        //! \copydoc ArrayRef::fromArray()
        static ConstArrayRef<value_type>
        fromArray(const value_type *begin, size_t size)
        {
            return ConstArrayRef<value_type>(begin, begin+size);
        }
        //! \copydoc ArrayRef::fromVector()
        static ConstArrayRef<value_type>
        fromVector(typename std::vector<value_type>::const_iterator begin,
                   typename std::vector<value_type>::const_iterator end)
        {
            const value_type *p_begin = (begin != end) ? &*begin : NULL;
            const value_type *p_end   = p_begin + (end-begin);
            return ConstArrayRef<value_type>(p_begin, p_end);
        }

        /*! \brief
         * Constructs an empty reference.
         */
        ConstArrayRef() : begin_(NULL), end_(NULL) {}
        /*! \brief
         * Constructs an empty reference.
         *
         * This is provided for convenience, such that EmptyArrayRef can be
         * used to initialize any Const ArrayRef, without specifying the
         * template type.  It is not explicit to enable that usage.
         */
        ConstArrayRef(const EmptyArrayRef &) : begin_(NULL), end_(NULL) {}
        /*! \brief
         * Constructs a reference to a particular range.
         *
         * \param[in] begin  Pointer to the beginning of a range.
         * \param[in] end    Pointer to the end of a range.
         *
         * Passed pointers must remain valid for the lifetime of this object.
         *
         * \note For clarity, use the non-member function constArrayRefFromPointers
         * instead.
         */
        ConstArrayRef(const_pointer begin, const_pointer end)
            : begin_(begin), end_(end)
        {
            GMX_ASSERT(end >= begin, "Invalid range");
        }
        /*! \brief
         * Constructs a reference to a whole vector.
         *
         * \param[in] v  Vector to reference.
         *
         * Passed vector must remain valid and not be reallocated for the
         * lifetime of this object.
         *
         * This constructor is not explicit to allow directly passing
         * std::vector to a method that takes ConstArrayRef.
         */
        ConstArrayRef(const std::vector<T> &v)
            : begin_((!v.empty()) ? &v[0] : NULL),
              end_((!v.empty()) ? &v[0] + v.size() : NULL)
        {
        }
        //! \cond
        // Doxygen 1.8.5 doesn't parse the declaration correctly...
        /*! \brief
         * Constructs a reference to a C array.
         *
         * \param[in] array  C array to reference.
         * \tparam    count  Deduced number of elements in \p array.
         *
         * This constructor can only be used with a real array (not with a
         * pointer).  It constructs a reference to the whole array, without
         * a need to pass the number of elements explicitly.  The compiler
         * must be able to deduce the array size.
         *
         * Passed array must remain valid for the lifetime of this object.
         *
         * This constructor is not explicit to allow directly passing
         * a C array to a function that takes a ConstArrayRef parameter.
         *
         * xlc on BG/Q compiles wrong code if the C array is a struct
         * field, unless value_type is char or unsigned char. There's
         * no good way to assert on this before C++11 (which that
         * compiler will never support).
         */
        template <size_t count>
        ConstArrayRef(const value_type (&array)[count])
            : begin_(array), end_(array + count)
        {
        }
        //! \endcond

        //! Returns an iterator to the beginning of the container.
        const_iterator begin() const { return begin_; }
        //! Returns an iterator to the end of the container.
        const_iterator end() const { return end_; }
        //! Returns an iterator to the reverse beginning of the container.
        const_iterator rbegin() const { return reverse_iterator(end()); }
        //! Returns an iterator to the reverse end of the container.
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


//! \copydoc ArrayRef::fromPointers()
//! \related ArrayRef
template <typename T>
ArrayRef<T> arrayRefFromPointers(T *begin, T *end)
{
    return ArrayRef<T>::fromPointers(begin, end);
}
//! \copydoc ArrayRef::fromArray()
//! \related ArrayRef
template <typename T>
ArrayRef<T> arrayRefFromArray(T *begin, size_t size)
{
    return ArrayRef<T>::fromArray(begin, size);
}
//! \copydoc ArrayRef::fromVector()
//! \related ArrayRef
template <typename T>
ArrayRef<T> arrayRefFromVector(typename std::vector<T>::iterator begin,
                               typename std::vector<T>::iterator end)
{
    return ArrayRef<T>::fromVector(begin, end);
}


//! \copydoc ConstArrayRef::fromPointers()
//! \related ConstArrayRef
template <typename T>
ConstArrayRef<T> constArrayRefFromPointers(const T *begin, const T *end)
{
    return ConstArrayRef<T>::fromPointers(begin, end);
}
//! \copydoc ConstArrayRef::fromArray()
//! \related ConstArrayRef
template <typename T>
ConstArrayRef<T> constArrayRefFromArray(const T *begin, size_t size)
{
    return ConstArrayRef<T>::fromArray(begin, size);
}
//! \copydoc ConstArrayRef::fromVector()
//! \related ConstArrayRef
template <typename T>
ConstArrayRef<T> constArrayRefFromVector(typename std::vector<T>::const_iterator begin,
                                         typename std::vector<T>::const_iterator end)
{
    return ConstArrayRef<T>::fromVector(begin, end);
}

/*! \brief
 * Simple swap method for ArrayRef objects.
 *
 * \see ArrayRef::swap()
 *
 * \ingroup module_utility
 */
template <typename T>
void swap(ArrayRef<T> &a, ArrayRef<T> &b)
{
    a.swap(b);
}

/*! \brief
 * Simple swap method for ConstArrayRef objects.
 *
 * \see ConstArrayRef::swap()
 *
 * \ingroup module_utility
 */
template <typename T>
void swap(ConstArrayRef<T> &a, ConstArrayRef<T> &b)
{
    a.swap(b);
}

} // namespace gmx

#endif
