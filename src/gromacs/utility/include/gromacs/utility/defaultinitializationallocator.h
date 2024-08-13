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
/*! \libinternal \file
 * \brief Declares an allocator that can use default initialization instead
 * of values initialization. This is useful for improving performance of
 * resize() in standard vectors for buffers in performance critical code.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_DEFAULTINITIALIZATIONALLOCATOR_H
#define GMX_UTILITY_DEFAULTINITIALIZATIONALLOCATOR_H

#include <memory>
#include <vector>

namespace gmx
{

/*! \libinternal \brief Allocator adaptor that interposes construct() calls to
 * convert value initialization into default initialization.
 *
 * This can be used to avoid initialization e.g. on resize() in std::vector.
 */
template<typename T, typename A = std::allocator<T>>
class DefaultInitializationAllocator : public A
{
    typedef std::allocator_traits<A> a_t;

public:
    template<typename U>
    struct rebind
    {
        using other = DefaultInitializationAllocator<U, typename a_t::template rebind_alloc<U>>;
    };

    using A::A;

    /*! \brief Constructs an object and default initializes
     *
     * \todo Use std::is_nothrow_default_constructible_v when CUDA 11 is a requirement.
     */
    template<typename U>
    void construct(U* ptr) noexcept(std::is_nothrow_default_constructible<U>::value)
    {
        ::new (static_cast<void*>(ptr)) U;
    }

    /*! \brief Constructs an object and value initializes */
    template<typename U, typename... Args>
    void construct(U* ptr, Args&&... args)
    {
        a_t::construct(static_cast<A&>(*this), ptr, std::forward<Args>(args)...);
    }
};

//! Convenience type for vector that avoids initialization at resize()
template<typename T>
using FastVector = std::vector<T, DefaultInitializationAllocator<T>>;

} // namespace gmx

#endif // GMX_UTILITY_DEFAULTINITIALIZATIONALLOCATOR_H
