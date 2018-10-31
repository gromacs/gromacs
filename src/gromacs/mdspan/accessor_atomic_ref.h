/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef MDSPAN_ACCESSOR_ATOMIC_REF
#define MDSPAN_ACCESSOR_ATOMIC_REF

#include <type_traits>

template<class T>
class atomic_ref
{
};

template<class T>
struct accessor_atomic {

    static_assert(std::is_trivially_copyable<T>::value, "Failed assert: trivially copyable");

    using element_type = T;
    using reference    = atomic_ref<T>;
    using pointer      = T*;
    using offset       = accessor_atomic;

    constexpr accessor_atomic() noexcept = default;
    constexpr accessor_atomic( accessor_atomic && ) noexcept      = default;
    constexpr accessor_atomic( const accessor_atomic & ) noexcept = default;
    accessor_atomic &operator=( accessor_atomic && ) noexcept      = default;
    accessor_atomic &operator=( const accessor_atomic & ) noexcept = default;

    explicit constexpr accessor_atomic( pointer other ) noexcept
        : ptr(other)
    { assert( 0 == reinterpret_cast<uintptr_t>(ptr) % reference::required_alignment ); };

    constexpr reference operator[]( size_t i ) const noexcept
    { return reference(ptr[i]); }

    constexpr offset operator+( size_t i ) const noexcept
    { return offset(ptr+i); }

    constexpr operator pointer() const
    { assert(false /* cannot access raw data outside of atomic */); }
};
#endif /* end of include guard: MDSPAN_ACCESSOR_ATOMIC_REF */
