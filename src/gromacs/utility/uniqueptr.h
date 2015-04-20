/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares gmx::gmx_unique_ptr and supporting functionality.
 *
 * \author Roland Schulz <roland@utk.edu>
 * \author John Eblen <jeblen@acm.org>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_UNIQUEPTR_H
#define GMX_UTILITY_UNIQUEPTR_H

#include "config.h"

#ifdef GMX_CXX11   // C++11 Compiler
#include <memory>  // IWYU pragma: export
#include <utility> // IWYU pragma: export
#else              // C++03 Compiler
#include <boost/shared_ptr.hpp>
#endif

namespace gmx
{

//! \cond libapi
/*! \libinternal \class gmx_unique_ptr
 * \brief
 * Smart pointer for unique ownership.
 *
 * The \a type member typedef declares the actual smart pointer type.
 * If std::unique_ptr from C++11 is available, it is used, otherwise maps to
 * boost::shared_ptr. Because of this, there are some limitations to usage.
 * gmx::move() should be used to move the pointer.
 *
 * Avoid using directly as a type, use a typedef instead. Typical usage:
 * \code
   typedef gmx_unique_ptr<ExampleClass>::type ExampleClassPointer;
   \endcode
 *
 * \ingroup module_utility
 * \inlibraryapi
 */
/*! \libinternal \typedef gmx_unique_ptr::type
 * \brief The smart pointer type.
 * Work-around for the non-existence of template typedefs in C++03.
 */
#ifdef GMX_CXX11 // C++11 Compiler
using std::move;
template<typename T>
struct gmx_unique_ptr
{
    typedef std::unique_ptr<T> type;
};
#else // C++03 Compiler
/*! \fn boost::shared_ptr<T> &move(boost::shared_ptr<T> &ptr)
 * \brief
 * Moves gmx::gmx_unique_ptr type pointers
 *
 * For C++11 gmx::move is the std::move, for non-C++11 compilers, the
 * move operation is a no-op.
 *
 * \inlibraryapi
 */
template<typename T>
boost::shared_ptr<T> &move(boost::shared_ptr<T> &ptr)
{
    return ptr;
}
template<typename T>
struct gmx_unique_ptr
{
    typedef boost::shared_ptr<T> type;
};
#endif
//! \endcond

} // namespace gmx

#endif
