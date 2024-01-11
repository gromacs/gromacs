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
/*! \libinternal \file
 * \brief Declares SimdArrayRef
 *
 * \author Roland Schulz <roland.schulz@intel.com>
 * \inlibraryapi
 * \ingroup module_simd
 */
#ifndef GMX_SIMD_SIMD_MEMORY_H
#define GMX_SIMD_SIMD_MEMORY_H

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{
namespace internal
{

template<typename T>
class SimdReference
{
private:
    using non_const_T = std::remove_const_t<T>;
    using pointer     = SimdTraitsT<T>*;

public:
    //! \brief Constructor
    explicit SimdReference(pointer m) : m_(m) {}
    //! \brief Conversion method that will execute load
    operator non_const_T() const { return load<non_const_T>(m_); }
    //! \brief Assignment operator that will execute store
    SimdReference operator=(T o) // NOLINT(misc-unconventional-assign-operator,cppcoreguidelines-c-copy-assignment-signature)
    {
        store(m_, o);
        return *this;
    }
    //! \brief Addition assignment operator that will execute load+store
    SimdReference operator+=(T o)
    {
        store(m_, load<non_const_T>(m_) + o);
        return *this;
    }
    //! \brief Subtraction assignment operator that will execute load+store
    SimdReference operator-=(T o)
    {
        store(m_, load<non_const_T>(m_) - o);
        return *this;
    }
    //! \brief Multiplication assignment operator that will execute load+store
    SimdReference operator*=(T o)
    {
        store(m_, load<non_const_T>(m_) * o);
        return *this;
    }

private:
    pointer m_; //!< The pointer used to load memory
};

template<typename T>
class SimdIterator :
    public gmx::boost::stl_interfaces::iterator_interface<SimdIterator<T>, std::random_access_iterator_tag, T, SimdReference<T>>
{
    using Base =
            gmx::boost::stl_interfaces::iterator_interface<SimdIterator<T>, std::random_access_iterator_tag, T, SimdReference<T>>;
    // pointer is T*
    using DataPointer = SimdTraitsT<T>*;

public:
    explicit SimdIterator(DataPointer p = 0) : p_(p)
    {
        GMX_ASSERT((reinterpret_cast<size_t>(p) / sizeof(*p)) % simdWidth == 0,
                   "Trying to create aligned iterator for non aligned address.");
    }
    SimdIterator& operator+=(typename Base::difference_type d)
    {
        p_ += simdWidth * d;
        return *this;
    }
    typename Base::difference_type operator-(SimdIterator o) { return (p_ - o.p_) / simdWidth; }
    typename Base::reference       operator*() const { return typename Base::reference(p_); }

private:
    DataPointer          p_;
    static constexpr int simdWidth = SimdTraits<T>::width;
};

/*! \internal
 * \brief STL-like container for aligned SIMD type. Used as ArrayRef<SimdReal>.
 *
 * Should provide the same interface as ArrayRef. Any missing functions should be
 * added as needed. The pointer type (used e.g. for initialization) is a real
 * pointer. The reference type (used e.g. for operator[] and iterator dereference)
 * is SimdReference which executes the aligned load/store as is appropriate. For
 * both iterator and element access, the access happens in blocks of SIMD width.
 * Meaning that a[1] refers to the 2nd SIMD vector and thus reals 8-15 for 8-wide
 * SIMD. The starting address has to be aligned and the length has to be multiple
 * of the SIMD width.
 *
 * \tparam T SIMD type (e.g. SimdReal)
 */
template<typename T>
class SimdArrayRef
{
public:
    //! Type for representing size of the container.
    using size_type = size_t;
    //! Type for representing difference between two container indices.
    using difference_type = std::ptrdiff_t;
    //! Type of values stored in the container.
    using value_type = T;
    //! Pointer to a container element.
    using pointer = SimdTraitsT<T>*;
    //! Reference to a container element.
    using reference = internal::SimdReference<T>;
    //! Iterator type for the container.
    using iterator = SimdIterator<T>;
    //! Standard reverse iterator.
    using reverse_iterator = std::reverse_iterator<iterator>;

    //! \copydoc ArrayRef::ArrayRef(pointer, pointer)
    SimdArrayRef(pointer begin, pointer end) : begin_(begin), end_(end)
    {
        GMX_ASSERT(end >= begin, "Invalid range");
        GMX_ASSERT((reinterpret_cast<size_type>(begin) / sizeof(*begin)) % simdWidth == 0,
                   "Aligned ArrayRef requires aligned starting address");
        GMX_ASSERT((reinterpret_cast<size_type>(end) / sizeof(*end)) % simdWidth == 0,
                   "Size of ArrayRef needs to be divisible by type size");
    }
    //! \copydoc ArrayRef::ArrayRef(U)
    template<typename U, typename = std::enable_if_t<std::is_convertible_v<typename std::remove_reference_t<U>::pointer, pointer>>>
    SimdArrayRef(U&& o) :
        begin_(reinterpret_cast<pointer>(o.data())),
        end_(reinterpret_cast<pointer>(o.data() + o.size()))
    {
    }
    // reinterpret_cast is only needed for const conversion of SimdArrayRef itself.
    // All other containers have type(o.data())==U::pointer (the cast does nothing).

    //! Returns the size of the container.
    size_type size() const { return (end_ - begin_) / simdWidth; }
    //! Whether the container is empty.
    bool empty() const { return begin_ == end_; }
    //! Returns an iterator to the beginning of the container.
    iterator begin() const { return iterator(begin_); }
    //! Returns an iterator to the end of the container.
    iterator end() const { return iterator(end_); }

    //! Access container element.
    reference operator[](size_type n) { return reference(begin_ + n * simdWidth); }

    //! Returns the first element in the container.
    reference front() const { return reference(begin_); }
    //! Returns the first element in the container.
    reference back() const { return reference(end_ - simdWidth); }

private:
    static constexpr int simdWidth = SimdTraits<T>::width;
    using pack_type                = SimdTraitsT<T>[simdWidth];
    // Private because dereferencing return value is undefined behavior (strict aliasing rule)
    // Only use is conversion constructor above which immediately casts it back.
    // Return type is not "pointer" because then data()+size() would be ill defined.
    // Has to be pack_type and not value_type in case
    // sizeof(value_type)/sizeof(pointer)!=simdWidth (e.g. int32 for double SSE2).
    pack_type* data() const { return reinterpret_cast<pack_type*>(begin_); }

    template<typename U>
    friend class SimdArrayRef;

    pointer begin_;
    pointer end_;
};

} // namespace internal

/* Specialize ArraryRef<SimdReal>
 * So far only an aligned version is implemented. The constructor verifies that
 * a ArrayRef<SimdReal> is constructed only for properly aligned data.
 */
#if GMX_SIMD_HAVE_FLOAT
template<>
class ArrayRef<SimdFloat> : public internal::SimdArrayRef<SimdFloat>
{
    using Base = internal::SimdArrayRef<SimdFloat>;
    using Base::Base;
};
template<>
class ArrayRef<const SimdFloat> : public internal::SimdArrayRef<const SimdFloat>
{
    using Base = internal::SimdArrayRef<const SimdFloat>;
    using Base::Base;
};
#endif
#if GMX_SIMD_HAVE_DOUBLE
template<>
class ArrayRef<SimdDouble> : public internal::SimdArrayRef<SimdDouble>
{
    using Base = internal::SimdArrayRef<SimdDouble>;
    using Base::Base;
};
template<>
class ArrayRef<const SimdDouble> : public internal::SimdArrayRef<const SimdDouble>
{
    using Base = internal::SimdArrayRef<const SimdDouble>;
    using Base::Base;
};
#endif

} // namespace gmx

#endif
