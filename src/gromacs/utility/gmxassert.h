/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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
 * Defines assert macros customized for Gromacs.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_GMXASSERT_H
#define GMX_UTILITY_GMXASSERT_H

#include <boost/current_function.hpp>

#include "gromacs/utility/basedefinitions.h"

//! \addtogroup module_utility
//! \{

/*! \def GMX_RELEASE_ASSERT
 * \brief
 * Macro for asserts that should also be present in the release version.
 *
 * Regardless of NDEBUG, this macro checks \p condition, and if it is not true,
 * it calls the assert handler.
 *
 * Although this macro currently calls abort() if the assertion fails, it
 * should only be used in a context where it is safe to throw an exception to
 * keep the option open.
 */
#ifdef GMX_DISABLE_ASSERTS
#define GMX_RELEASE_ASSERT(condition, msg)
#else
#define GMX_RELEASE_ASSERT(condition, msg) \
    ((void) ((condition) ? (void)0 : \
             ::gmx::internal::assertHandler(#condition, msg, \
                                            BOOST_CURRENT_FUNCTION, __FILE__, __LINE__)))
#endif
/*! \def GMX_ASSERT
 * \brief
 * Macro for debug asserts.
 *
 * If NDEBUG is defined, this macro expands to nothing.
 * If it is not defined, it will work exactly like ::GMX_RELEASE_ASSERT.
 *
 * \see ::GMX_RELEASE_ASSERT
 */
#ifdef NDEBUG
#define GMX_ASSERT(condition, msg)
#else
#define GMX_ASSERT(condition, msg) GMX_RELEASE_ASSERT(condition, msg)
#endif

//! \}

namespace gmx
{

/*! \cond internal */
namespace internal
{

/*! \brief
 * Called when an assert fails.
 *
 * Should not be called directly, but instead through ::GMX_ASSERT or
 * ::GMX_RELEASE_ASSERT.
 *
 * \ingroup module_utility
 */
gmx_noreturn
void assertHandler(const char *condition, const char *msg,
                   const char *func, const char *file, int line);

}   // namespace internal
//! \endcond

} // namespace gmx

#endif
