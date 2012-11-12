/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \file
 * \brief
 * Defines assert macros customized for Gromacs.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_GMXASSERT_H
#define GMX_UTILITY_GMXASSERT_H

#include <boost/current_function.hpp>
#include <boost/exception/detail/attribute_noreturn.hpp>

/*! \addtopublicapi
 * \{
 */

/*! \brief
 * Macro for asserts that should also be present in the release version.
 *
 * Regardless of NDEBUG, this macro checks \p condition, and if it is not true,
 * it calls the assert handler.
 *
 * Although this macro currently calls abort() if the assertion fails, it
 * should only be used in a context where it is safe to throw an exception to
 * keep the option open.
 */
#define GMX_RELEASE_ASSERT(condition, msg) \
    ((void) ((condition) ? (void)0 : \
                 ::gmx::internal::assertHandler(# condition, msg, \
                                                BOOST_CURRENT_FUNCTION, __FILE__, __LINE__)))
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

/*!\}*/

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
BOOST_ATTRIBUTE_NORETURN
void assertHandler(const char *condition, const char *msg,
                   const char *func, const char *file, int line);

}   // namespace internal
//! \endcond

} // namespace gmx

#endif
