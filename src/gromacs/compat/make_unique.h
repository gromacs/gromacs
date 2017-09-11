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
 * \defgroup module_compat C++ standard library compatibility helpers.
 * \brief Provide uniform interface to selected C++ standard library features.
 *
 * For some features not available on all platforms supported by
 * \Gromacs, provide back-ports or mappings to available standard
 * library implementations as appropriate.
 */

/*! \libinternal
 * \file
 * \brief Provides template gmx::compat::make_unique
 *
 * The implementation of gmx::compat::make_unique is copied with little
 * modification from C++ standardization doc N3656
 * at http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3656.htm
 * though additional wrapping has been added for use in \Gromacs.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup module_compat
 * \inlibraryapi
 */
/*! \addtogroup module_compat
 * ### gmx::compat::make_unique
 *
 * Provide a trivially adapted implementation of the C++ standard library `make_unique` function template.
 * When All supported \Gromacs build platforms provide `std::make_unique`, this should be removed.
 *
 */
#include <cstddef>

#include <memory>
#include <type_traits>
#include <utility>

namespace gmx
{
namespace compat
{

///\cond

// All gmx::compat code should use std::unique_ptr
using ::std::unique_ptr;

template<class T> struct Unique_if {
    typedef unique_ptr<T> Single_object;
};

template<class T> struct Unique_if<T[]> {
    typedef unique_ptr<T[]> Unknown_bound;
};

template<class T, size_t N> struct Unique_if<T[N]> {
    typedef void Known_bound;
};

template<class T, class ... Args>
typename Unique_if<T>::Single_object
make_unique(Args && ... args)
{
    return unique_ptr<T>(new T(::std::forward<Args>(args) ...));
}

template<class T>
typename Unique_if<T>::Unknown_bound
make_unique(size_t n)
{
    typedef typename ::std::remove_extent<T>::type U;
    return unique_ptr<T>(new U[n]());
}

template<class T, class ... Args>
typename Unique_if<T>::Known_bound
make_unique(Args && ...) = delete;

///\endcond

}      // namespace gmx::compat
}      // namespace gmx

#endif // header guard
