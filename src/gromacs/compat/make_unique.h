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
#ifndef GMX_COMPAT_MAKE_UNIQUE_H
#define GMX_COMPAT_MAKE_UNIQUE_H
/*! \libinternal
 * \file
 * \brief Provides template gmx::compat::make_unique
 *
 * The implementation of gmx::compat::make_unique is copied with little
 * modification from C++ standardization doc N3656
 * at http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3656.htm
 * though additional wrapping has been added for use in Gromacs.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup module_compatibility
 */
/*! \addtogroup module_compatibility
 * ### gmx::compat::make_unique
 *
 * * For C++ standard greater than or equal to C++14, `gmx::compat::make_unique`
 *   is an alias to `std::make_unique`.
 * * For C++ standard less than C++14,
 *   borrow `std::make_unique` implementation from C++ standardization doc N3656
 *   at http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3656.htm
 *
 */
#ifdef __cplusplus
#include <memory>

namespace gmx
{
namespace compat
{

///\cond
#if __cplusplus < 201402L

#include <cstddef>
#include <type_traits>
#include <utility>

// All gmx::compat code should use std::unique_ptr
using std::unique_ptr;

template<class T> struct _Unique_if {
    typedef unique_ptr<T> _Single_object;
};

template<class T> struct _Unique_if<T[]> {
    typedef unique_ptr<T[]> _Unknown_bound;
};

template<class T, size_t N> struct _Unique_if<T[N]> {
    typedef void _Known_bound;
};

template<class T, class ... Args>
typename _Unique_if<T>::_Single_object
make_unique(Args && ... args)
{
    return unique_ptr<T>(new T(std::forward<Args>(args) ...));
}

template<class T>
typename _Unique_if<T>::_Unknown_bound
make_unique(size_t n)
{
    typedef typename std::remove_extent<T>::type U;
    return unique_ptr<T>(new U[n]());
}

template<class T, class ... Args>
typename _Unique_if<T>::_Known_bound
make_unique(Args && ...) = delete;

#else

// Forward gmx::compat::make_unique<Explicit...>(Args...) to implementation in STL.
template<typename ... Explicit, typename ... Args>
auto make_unique(Args ... args)->decltype(std::make_unique<Explicit...>(std::forward<Args>(args) ...))
{
    return std::make_unique<Explicit...>(std::forward<Args>(args) ...);
}

#endif  // __cplusplus < 201402L
///\endcond

}      // namespace gmx::compat
}      // namespace gmx

#endif // defined __cplusplus
#endif // header guard
