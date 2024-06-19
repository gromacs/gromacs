/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Declares gmx::FixedCapacityVector
 *
 * \author Berk Hess <hess@kth.se>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_FIXEDCAPACITYVECTOR_H
#define GMX_UTILITY_FIXEDCAPACITYVECTOR_H

#include <array>
#include <stdexcept>
#include <type_traits>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief Vector that behaves likes std::vector but has fixed capacity.
 *
 * \tparam T         Value type of elements, should be default constructible
 * \tparam capacity_ The maximum number of elements that can be stored.
 *
 * This class provides a variable size container, but with constant
 * memory usage and can be allocated on the stack and avoid the overhead
 * of dynamic allocation. This is especially useful for small vectors
 * which are set up frequently.
 *
 * The class supports all methods from \p std::array, but behaves more
 * like \p std::vector since it has variable size. In addition to the methods
 * from std::array, from \p std::vector the methods \p push_back(), \p pop_back(),
 * emplace_back() and \p clear() are supported. In particular, methods that
 * requires reordering, such as \p insert() and \p emplace() are not
 * supported to keep the code simple.
 *
 * The size is 0 at construction and elements can only be added with
 * \p push_back() and \p emplace_back().
 *
 * \note This class is very similar to the fixed_capacity_vector class
 * proposed for the C++ standard in document P0843r see:
 * http://open-std.org/JTC1/SC22/WG21/docs/papers/2018/p0843r1.html
 *
 * \inpublicapi
 * \ingroup module_utility
 */

template<typename T, size_t capacity_>
class FixedCapacityVector
{
    static_assert(std::is_default_constructible_v<T>);

public:
    //! Type of values stored in the vector
    using value_type = T;
    //! Type for representing size of the vector
    using size_type = size_t;
    //! Type for representing difference between two indices
    using difference_type = ptrdiff_t;
    //! Const reference to an element
    using const_reference = const T&;
    //! Const pointer to an element
    using const_pointer = const T*;
    //! Const iterator type to an element
    using const_iterator = const T*;
    //! Reference to an element
    using reference = T&;
    //! Pointer to an element
    using pointer = T*;
    //! Iterator type to an element
    using iterator = T*;
    //! Standard reverse iterator
    using reverse_iterator = std::reverse_iterator<iterator>;
    //! Standard reverse iterator
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    //! Returns a const iterator to the beginning
    const_iterator begin() const noexcept { return data(); }
    //! Returns an iterator to the beginning
    iterator begin() noexcept { return data(); }
    //! Returns a const iterator to the end
    const_iterator end() const noexcept { return data_.data() + size_; }
    //! Returns an iterator to the end
    iterator end() noexcept { return data_.data() + size_; }
    //! Returns a const iterator to the reverse beginning
    const_reverse_iterator rbegin() const noexcept { return reverse_iterator(end()); }
    //! Returns an iterator to the reverse beginning
    reverse_iterator rbegin() noexcept { return reverse_iterator(end()); }
    //! Returns a const iterator to the reverse end
    const_reverse_iterator rend() const noexcept { return reverse_iterator(begin()); }
    //! Returns an iterator to the reverse end
    reverse_iterator rend() noexcept { return reverse_iterator(begin()); }

    /*! \brief Returns the size
     *
     * \note Use ssize for any expression involving arithmetic operations
         (including loop indices).
     */
    constexpr size_type size() const noexcept { return size_; }
    //! Returns the signed size
    constexpr Index ssize() const noexcept { return size_; }
    //! Returns whether the vector is empty
    constexpr bool empty() const noexcept { return size_ == 0; }

    //! Returns the vector capacity (max. number of elements that can be stored)
    static constexpr size_type max_size() noexcept { return capacity_; }
    //! Returns the vector capacity (max. number of elements that can be stored)
    static constexpr size_type capacity() noexcept { return capacity_; }

    //! Const access an element
    constexpr const_reference operator[](size_type n) const noexcept
    {
        GMX_ASSERT(n < size(), "Index should be in range");
        return data_[n];
    }
    //! Access an element
    constexpr reference operator[](size_type n) noexcept
    {
        GMX_ASSERT(n < size(), "Index should be in range");
        return data_[n];
    }
    //! Const access an element, throws an out_of_range exception when out of range
    constexpr const_reference at(size_type n) const
    {
        if (n >= size())
        {
            throw std::out_of_range("Vector index out of range");
        }
        return data_[n];
    }
    //! Access an element, throws an out_of_range exception when out of range
    constexpr reference at(size_type n)
    {
        if (n >= size())
        {
            throw std::out_of_range("Vector index out of range");
        }
        return data_[n];
    }
    //! Returns the first element
    //! Returns the first element
    constexpr reference front() noexcept { return data_.front(); }
    //! Returns the first element (const version)
    constexpr const_reference front() const noexcept { return data_.front(); }
    //! Returns the last element
    constexpr reference back() noexcept { return data_[size_ - 1]; }
    //! Returns the last element (const version)
    constexpr const_reference back() const noexcept { return data_[size_ - 1]; }

    //! Returns a raw pointer to the contents of the array
    constexpr const T* data() const noexcept { return data_.data(); }

    //! Returns a raw pointer to the contents of the array
    constexpr T* data() noexcept { return data_.data(); }

    //! Adds element at the end
    constexpr void push_back(const T& value) noexcept
    {
        GMX_ASSERT(size() < capacity_, "Cannot add more elements than the capacity");
        data_[size_] = value;
        size_++;
    }

    //! Deletes last element
    constexpr void pop_back() noexcept
    {
        GMX_ASSERT(!empty(), "Can only delete last element when present");
        if constexpr (!std::is_trivially_destructible_v<T>)
        {
            ~back();
        }
        size_--;
    }

    //! Constructs an element at the end
    template<class... Args>
    constexpr reference emplace_back(Args&&... args)
    {
        GMX_ASSERT(size() < capacity_, "Cannot add more elements than the capacity");
        if constexpr (std::is_move_assignable<T>::value)
        {
            data_[size_] = std::move(T(args...));
        }
        else
        {
            data_[size_] = T(args...);
        }
        size_++;

        return back();
    }

    //! Clears content
    constexpr void clear() noexcept
    {
        if constexpr (!std::is_trivially_destructible_v<T>)
        {
            for (auto& entry : *this)
            {
                ~entry;
            }
        }

        size_ = 0;
    }

private:
    //! The elements, stored in a fixed size array
    std::array<T, capacity_> data_;
    //! The size of the vector
    size_type size_ = 0;
};

} // namespace gmx

#endif
