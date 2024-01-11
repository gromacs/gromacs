/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 * \brief
 * Declares gmx::ArrayRef
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Roland Schulz <roland.schulz@intel.com>
 * \author Berk Hess <hess@kth.se>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_ARRAYREF_H
#define GMX_UTILITY_ARRAYREF_H

#include <cassert>
#include <cstddef>

#include <array>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>

#if __has_include("external/boost/stl_interfaces/iterator_interface.hpp")
#    include "external/boost/stl_interfaces/iterator_interface.hpp"
#else // fallback for installed headers
#    include <gromacs/external/boost/stl_interfaces/iterator_interface.hpp>
#endif

namespace gmx
{

template<class T>
struct ArrayRefIter :
    gmx::boost::stl_interfaces::iterator_interface<ArrayRefIter<T>, std::random_access_iterator_tag, T>
{
    // This default constructor does not initialize it_
    constexpr ArrayRefIter() noexcept {}
    constexpr explicit ArrayRefIter(T* it) noexcept : it_(it) {}
    // TODO: Use std::is_const_v when CUDA 11 is a requirement.
    template<class T2 = T, class = std::enable_if_t<std::is_const<T2>::value>>
    constexpr ArrayRefIter(ArrayRefIter<std::remove_const_t<T2>> it) noexcept : it_(&*it)
    {
    }
    constexpr T*            data() const noexcept { return it_; }
    constexpr T&            operator*() const noexcept { return *it_; }
    constexpr ArrayRefIter& operator+=(std::ptrdiff_t i) noexcept
    {
        it_ += i;
        return *this;
    }
    constexpr auto operator-(ArrayRefIter other) const noexcept { return it_ - other.it_; }

private:
    T* it_ = nullptr;
};

/*! \brief STL-like interface to a C array of T (or part
 * of a std container of T).
 *
 * \tparam T  Value type of elements.
 *
 * This class provides an interface similar to \c std::vector<T, A>, with the
 * following main differences:
 *  - This class does not have its own storage.  Instead, it references an
 *    existing array of values (either a C-style array or part of an existing
 *    std::vector<T, A> or std::array<T>).
 *  - It is only possible to modify the values themselves through ArrayRef;
 *    it is not possible to add or remove values.
 *  - Copying objects of this type is cheap, and the copies behave identically
 *    to the original object: the copy references the same set of values.
 *
 * This class is useful for writing wrappers that expose a view of the
 * internal data stored as a single vector/array, which can be a whole
 * or part of the underlying storage.
 *
 * Methods in this class do not throw, except where indicated.
 *
 * Note that due to a Doxygen limitation, the constructor that takes a C array
 * whose size is known at compile time does not appear in the documentation.
 *
 * To refer to const data of type T, ArrayRef<const T> is used. For both const
 * and non-const std::vector and std::array an ArrayRef view can be created.
 * Attempting to create a non-const ArrayRef of a const vector/array will result
 * in a compiler error in the respective constructor.
 *
 * For SIMD types there is template specialization available
 * (e.g. ArrayRef<SimdReal>) in gromacs/simd/simd_memory.h which should have
 * the same functionality as much as possible.
 *
 * \todo
 * This class is not complete. There are likely also methods missing (not
 * required for current usage).
 *
 * \inpublicapi
 * \ingroup module_utility
 */
template<typename T>
class ArrayRef
{
public:
    //! Type of values stored in the reference.
    typedef T value_type;
    //! Type for representing size of the reference.
    typedef size_t size_type;
    //! Type for representing difference between two indices.
    typedef ptrdiff_t difference_type;
    //! Const reference to an element.
    typedef const T& const_reference;
    //! Const pointer to an element.
    typedef const T* const_pointer;
    //! Const iterator type to an element.
    typedef ArrayRefIter<const T> const_iterator;
    //! Reference to an element.
    typedef T& reference;
    //! Pointer to an element.
    typedef T* pointer;
    //! Iterator type to an element.
    typedef ArrayRefIter<T> iterator;
    //! Standard reverse iterator.
    typedef std::reverse_iterator<iterator> reverse_iterator;
    //! Standard reverse iterator.
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    /*! \brief
     * Constructs an empty reference.
     */
    ArrayRef() : begin_(nullptr), end_(nullptr) {}
    /*! \brief
     * Constructs a reference to a container or reference
     *
     * \param[in] o container to reference.
     *
     * Can be used to create a reference to a whole vector, std::array or
     * an ArrayRef. The destination has to have a convertible pointer type
     * (identical besides const or base class).
     *
     * Passed container must remain valid and not be reallocated for the
     * lifetime of this object.
     *
     * This constructor is not explicit to allow directly passing
     * a container to a method that takes ArrayRef.
     *
     * \todo Use std::is_convertible_v when CUDA 11 is a requirement.
     */
    template<typename U, typename = std::enable_if_t<std::is_convertible<typename std::remove_reference_t<U>::pointer, pointer>::value>>
    ArrayRef(U&& o) : begin_(o.data()), end_(o.data() + o.size())
    {
    }
    /*! \brief
     * Constructs a reference to a particular range.
     *
     * \param[in] begin  Pointer to the beginning of a range.
     * \param[in] end    Pointer to the end of a range.
     *
     * Passed pointers must remain valid for the lifetime of this object.
     */
    ArrayRef(pointer begin, pointer end) : begin_(begin), end_(end)
    {
        assert((end >= begin && "Invalid range"));
        assert((begin != nullptr || (begin == nullptr && end == nullptr))
               && "If begin is nullptr, end needs to be nullptr as well");
    }
    /*! \brief
     * Constructs a reference to a particular range.
     *
     * \param[in] begin  Iterator to the beginning of a range.
     * \param[in] end    iterator to the end of a range.
     *
     * Passed iterators must remain valid for the lifetime of this object.
     */
    ArrayRef(iterator begin, iterator end) : begin_(begin), end_(end)
    {
        assert((end >= begin && "Invalid range"));
        assert((begin.data() != nullptr || (begin.data() == nullptr && end.data() == nullptr))
               && "If begin is nullptr, end needs to be nullptr as well");
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
     */
    template<size_t count>
    ArrayRef(value_type (&array)[count]) : begin_(array), end_(array + count)
    {
    }
    //! \endcond

    //! Returns a reference to part of the memory.
    ArrayRef subArray(size_type start, size_type count) const
    {
        return { begin_ + start, begin_ + start + count };
    }
    //! Returns an iterator to the beginning of the reference.
    iterator begin() const { return iterator(begin_); }
    //! Returns an iterator to the end of the reference.
    iterator end() const { return iterator(end_); }
    //! Returns an iterator to the reverse beginning of the reference.
    reverse_iterator rbegin() const { return reverse_iterator(end()); }
    //! Returns an iterator to the reverse end of the reference.
    reverse_iterator rend() const { return reverse_iterator(begin()); }

    /*! \brief Returns the size of the reference.
     *
     * \note Use ssize for any expression involving arithmetic operations
         (including loop indices).
     */
    size_type size() const { return end_ - begin_; }
    //! Returns the signed size of the reference.
    difference_type ssize() const { return size(); }
    //! Identical to size().
    size_type capacity() const { return end_ - begin_; }
    //! Whether the reference refers to no memory.
    bool empty() const { return begin_ == end_; }

    //! Access an element.
    reference operator[](size_type n) const { return begin_[n]; }
    //! Access an element (throws on out-of-range error).
    reference at(size_type n) const
    {
        if (n >= size())
        {
            throw std::out_of_range("Vector index out of range");
        }
        return begin_[n];
    }
    //! Returns the first element.
    reference front() const { return *(begin_); }
    //! Returns the first element.
    reference back() const { return *(end_ - 1); }

    //! Returns a raw pointer to the contents of the array.
    pointer data() const { return begin_.data(); }

    /*! \brief
     * Swaps referenced memory with the other object.
     *
     * The actual memory areas are not modified, only the references are
     * swapped.
     */
    void swap(ArrayRef<T>& other)
    {
        std::swap(begin_, other.begin_);
        std::swap(end_, other.end_);
    }

private:
    iterator begin_;
    iterator end_;
};

/*! \brief
 * Constructs a reference to a C array.
 *
 * \param[in] begin  Pointer to the beginning of array.
 * \param[in] size   Number of elements in array.
 *
 * Passed array must remain valid for the lifetime of this object.
 * If \c begin is nullptr, return an empty ArrayRef.
 */
//! \related ArrayRef
template<typename T>
ArrayRef<T> arrayRefFromArray(T* begin, size_t size)
{
    return (begin != nullptr) ? ArrayRef<T>(begin, begin + size) : ArrayRef<T>{};
}

//! \copydoc arrayRefFromArray
//! \related ArrayRef
template<typename T>
ArrayRef<const T> constArrayRefFromArray(const T* begin, size_t size)
{
    return (begin != nullptr) ? ArrayRef<const T>(begin, begin + size) : ArrayRef<const T>{};
}

/*! \brief
 * Create ArrayRef from container with type deduction
 *
 * \see ArrayRef
 *
 * \todo Use std::is_const_v when CUDA 11 is a requirement.
 */
template<typename T>
ArrayRef<std::conditional_t<std::is_const<T>::value, const typename T::value_type, typename T::value_type>>
makeArrayRef(T& c)
{
    return c;
}

/*! \brief
 * Create ArrayRef to const T from container with type deduction
 *
 * \see ArrayRef
 */
template<typename T>
ArrayRef<const typename T::value_type> makeConstArrayRef(const T& c)
{
    return c;
}

/*! \brief
 * Simple swap method for ArrayRef objects.
 *
 * \see ArrayRef::swap()
 *
 * \ingroup module_utility
 */
template<typename T>
void swap(ArrayRef<T>& a, ArrayRef<T>& b)
{
    a.swap(b);
}

/*! \brief Return a vector that is a copy of an ArrayRef.
 *
 * This makes it convenient, clear, and performant (the compiler will
 * either do RVO to elide the temporary, or invoke the move constructor
 * taking the unnamed temporary) to write a declaration like
 *
 *   auto v = copyOf(arrayRef);
 *
 * \ingroup module_utility
 */
template<typename T>
std::vector<T> copyOf(const ArrayRef<const T>& arrayRef)
{
    return std::vector<T>(arrayRef.begin(), arrayRef.end());
}

} // namespace gmx

#endif
