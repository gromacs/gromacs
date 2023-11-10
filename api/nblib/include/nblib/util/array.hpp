/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief  implementation of a compile-time size array that can be used on the host and device
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_UTIL_ARRAY_HPP
#define NBLIB_UTIL_ARRAY_HPP

#include <cmath>

#include <iterator>
#include <utility>

#include "nblib/util/annotation.hpp"

namespace util
{

template<class T>
constexpr int determineAlignment(int n)
{
    if (sizeof(T) * n % 16 == 0)
    {
        return 16;
    }
    else if (sizeof(T) * n % 8 == 0)
    {
        return 8;
    }
    else
    {
        return alignof(T);
    }
}

/*! \brief std::array-like compile-time size array
 * \tparam T element type
 * \tparam N number of elements
 *
 * The implementation corresponds to a device-qualified std::array minus support for length 0
 * plus arithmetic operations.
 */
template<class T, std::size_t N>
struct alignas(determineAlignment<T>(N)) array
{
    typedef T                                     value_type;
    typedef value_type*                           pointer;
    typedef const value_type*                     const_pointer;
    typedef value_type&                           reference;
    typedef const value_type&                     const_reference;
    typedef value_type*                           iterator;
    typedef const value_type*                     const_iterator;
    typedef std::size_t                           size_type;
    typedef std::ptrdiff_t                        difference_type;
    typedef std::reverse_iterator<iterator>       reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    // public for aggregate type, as in std::array
    T data_[N];

    // No explicit construct/copy/destroy for aggregate type.

    HOST_DEVICE_FUN constexpr iterator begin() noexcept { return iterator(data()); }

    HOST_DEVICE_FUN constexpr const_iterator begin() const noexcept
    {
        return const_iterator(data());
    }

    HOST_DEVICE_FUN constexpr iterator end() noexcept { return iterator(data() + N); }

    HOST_DEVICE_FUN constexpr const_iterator end() const noexcept
    {
        return const_iterator(data() + N);
    }

    HOST_DEVICE_FUN constexpr reverse_iterator rbegin() noexcept { return reverse_iterator(end()); }

    HOST_DEVICE_FUN constexpr const_reverse_iterator rbegin() const noexcept
    {
        return const_reverse_iterator(end());
    }

    HOST_DEVICE_FUN constexpr reverse_iterator rend() noexcept { return reverse_iterator(begin()); }

    HOST_DEVICE_FUN constexpr const_reverse_iterator rend() const noexcept
    {
        return const_reverse_iterator(begin());
    }

    HOST_DEVICE_FUN constexpr const_iterator cbegin() const noexcept
    {
        return const_iterator(data());
    }

    HOST_DEVICE_FUN constexpr const_iterator cend() const noexcept
    {
        return const_iterator(data() + N);
    }

    HOST_DEVICE_FUN constexpr const_reverse_iterator crbegin() const noexcept
    {
        return const_reverse_iterator(end());
    }

    HOST_DEVICE_FUN constexpr const_reverse_iterator crend() const noexcept
    {
        return const_reverse_iterator(begin());
    }

    HOST_DEVICE_FUN constexpr size_type size() const noexcept { return N; }

    HOST_DEVICE_FUN constexpr size_type max_size() const noexcept { return N; }

    [[nodiscard]] HOST_DEVICE_FUN constexpr bool empty() const noexcept { return size() == 0; }

    // Element access.
    HOST_DEVICE_FUN constexpr reference operator[](size_type n) noexcept { return data_[n]; }

    HOST_DEVICE_FUN constexpr const_reference operator[](size_type n) const noexcept
    {
        return data_[n];
    }

    HOST_DEVICE_FUN constexpr reference front() noexcept { return *begin(); }

    HOST_DEVICE_FUN constexpr const_reference front() const noexcept { return data_[0]; }

    HOST_DEVICE_FUN constexpr reference back() noexcept { return *(end() - 1); }

    HOST_DEVICE_FUN constexpr const_reference back() const noexcept { return data_[N - 1]; }

    HOST_DEVICE_FUN constexpr pointer data() noexcept { return data_; }

    HOST_DEVICE_FUN constexpr const_pointer data() const noexcept { return data_; }

    HOST_DEVICE_FUN constexpr array<T, N>& operator+=(const array<T, N>& rhs) noexcept
    {
        auto add = [](T a, T b) { return a + b; };
        assignImpl(data(), rhs.data(), add, std::make_index_sequence<N>{});
        return *this;
    }

    HOST_DEVICE_FUN constexpr array<T, N>& operator-=(const array<T, N>& rhs) noexcept
    {
        auto minus = [](T a, T b) { return a - b; };
        assignImpl(data(), rhs.data(), minus, std::make_index_sequence<N>{});
        return *this;
    }

    HOST_DEVICE_FUN constexpr array<T, N>& operator*=(const value_type& rhs) noexcept
    {
        auto mult = [](T a, T b) { return a * b; };
        assignImpl(data(), rhs, mult, std::make_index_sequence<N>{});
        return *this;
    }

    HOST_DEVICE_FUN constexpr array<T, N>& operator/=(const value_type& rhs) noexcept
    {
        auto divide = [](T a, T b) { return a / b; };
        assignImpl(data(), rhs, divide, std::make_index_sequence<N>{});
        return *this;
    }

    //! conversion to pointer for interoperability with gmx::BasicVector
    HOST_DEVICE_FUN constexpr T* as_vec() { return data_; }
    //! conversion to const pointer for interoperability with gmx::BasicVector
    HOST_DEVICE_FUN constexpr const T* as_vec() const { return data_; }


private:
    template<class F, std::size_t... Is>
    HOST_DEVICE_FUN constexpr static void assignImpl(T* a, const T* b, F&& f, std::index_sequence<Is...>) noexcept
    {
        [[maybe_unused]] std::initializer_list<int> list{ (a[Is] = f(a[Is], b[Is]), 0)... };
    }

    template<class F, std::size_t... Is>
    HOST_DEVICE_FUN constexpr static void assignImpl(T* a, const T& b, F&& f, std::index_sequence<Is...>) noexcept
    {
        [[maybe_unused]] std::initializer_list<int> list{ (a[Is] = f(a[Is], b), 0)... };
    }
};

template<std::size_t I, class T, std::size_t N>
HOST_DEVICE_FUN constexpr T& get(array<T, N>& a_)
{
    return a_[I];
}

template<std::size_t I, class T, std::size_t N>
HOST_DEVICE_FUN constexpr const T& get(const array<T, N>& a_)
{
    return a_[I];
}

template<class T, std::size_t N>
HOST_DEVICE_FUN constexpr array<T, N> operator+(const array<T, N>& a, const array<T, N>& b)
{
    auto ret = a;
    return ret += b;
}

namespace detail
{

template<class T, std::size_t... Is>
HOST_DEVICE_FUN constexpr array<T, sizeof...(Is)> negateImpl(const array<T, sizeof...(Is)>& a,
                                                             std::index_sequence<Is...>)
{
    return array<T, sizeof...(Is)>{ -a[Is]... };
}

} // namespace detail

template<class T, std::size_t N>
HOST_DEVICE_FUN constexpr array<T, N> operator-(const array<T, N>& a)
{
    return detail::negateImpl(a, std::make_index_sequence<N>{});
}

template<class T, std::size_t N>
HOST_DEVICE_FUN constexpr array<T, N> operator-(const array<T, N>& a, const array<T, N>& b)
{
    auto ret = a;
    return ret -= b;
}

template<class T, std::size_t N>
HOST_DEVICE_FUN constexpr array<T, N> operator*(const array<T, N>& a, const T& b)
{
    auto ret = a;
    return ret *= b;
}

template<class T, std::size_t N>
HOST_DEVICE_FUN constexpr array<T, N> operator*(const T& a, const array<T, N>& b)
{
    auto ret = b;
    return ret *= a;
}

namespace detail
{

template<class T, std::size_t... Is>
HOST_DEVICE_FUN constexpr bool eqImpl(const T* a, const T* b, std::index_sequence<Is...>)
{
    return ((a[Is] == b[Is]) && ...);
}

template<class T, std::size_t... Is>
HOST_DEVICE_FUN constexpr T dotImpl(const T* a, const T* b, std::index_sequence<Is...>)
{
    return ((a[Is] * b[Is]) + ...);
}

template<int N, int I = 0>
struct LexicographicalCompare
{
    template<class T, class F>
    HOST_DEVICE_FUN constexpr static auto loop(const T* lhs, const T* rhs, F&& compare)
    {
        if (compare(lhs[I], rhs[I]))
        {
            return true;
        }
        if (compare(rhs[I], lhs[I]))
        {
            return false;
        }
        return LexicographicalCompare<N, I + 1>::loop(lhs, rhs, compare);
    }
};

template<int N>
struct LexicographicalCompare<N, N>
{
    template<class T, class F>
    HOST_DEVICE_FUN constexpr static auto loop(const T*, const T*, F&&)
    {
        return false;
    }
};

} // namespace detail

template<class T, std::size_t N>
HOST_DEVICE_FUN constexpr bool operator==(const array<T, N>& a, const array<T, N>& b)
{
    return detail::eqImpl(a.data(), b.data(), std::make_index_sequence<N>{});
}

template<class T, std::size_t N>
HOST_DEVICE_FUN constexpr bool operator!=(const array<T, N>& a, const array<T, N>& b)
{
    return !(a == b);
}

template<class T, std::size_t N>
HOST_DEVICE_FUN constexpr bool operator<(const array<T, N>& a, const array<T, N>& b)
{
    auto less = [](T a, T b) { return a < b; };
    return detail::LexicographicalCompare<N>::loop(a.data(), b.data(), less);
}

template<class T, std::size_t N>
HOST_DEVICE_FUN constexpr bool operator>(const array<T, N>& a, const array<T, N>& b)
{
    auto greater = [](T a, T b) { return a > b; };
    return detail::LexicographicalCompare<N>::loop(a.data(), b.data(), greater);
}

template<class T, std::size_t N>
HOST_DEVICE_FUN constexpr T dot(const array<T, N>& a, const array<T, N>& b)
{
    return detail::dotImpl(a.data(), b.data(), std::make_index_sequence<N>{});
}

template<class T, std::size_t N>
HOST_DEVICE_FUN constexpr T norm2(const array<T, N>& a)
{
    return dot(a, a);
}

template<class T, std::size_t N>
HOST_DEVICE_FUN constexpr T norm(const array<T, N>& a)
{
    return std::sqrt(norm2(a));
}

template<class T>
HOST_DEVICE_FUN constexpr array<T, 3> cross(const array<T, 3>& a, const array<T, 3>& b)
{
    return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
}

} // namespace util

template<class Vector, class = void>
struct VectorValueType
{
    using type = typename Vector::value_type;
};

template<class Vector>
struct VectorValueType<Vector, std::void_t<typename Vector::RawArray>>
{
    using type = std::remove_all_extents_t<typename Vector::RawArray>;
};

template<class Vector>
using VectorValueType_t = typename VectorValueType<Vector>::type;

#endif // NBLIB_UTIL_ARRAY_HPP
