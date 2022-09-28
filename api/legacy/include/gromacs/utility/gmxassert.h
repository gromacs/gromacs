/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2011- The GROMACS Authors
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
 * Defines assert macros customized for Gromacs.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_GMXASSERT_H
#define GMX_UTILITY_GMXASSERT_H

#include <filesystem>

#include "current_function.h"

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
#    define GMX_RELEASE_ASSERT(condition, msg)
#else
#    ifdef _MSC_VER
#        define GMX_RELEASE_ASSERT(condition, msg)                \
            ((void)((condition) ? (void)0                         \
                                : ::gmx::internal::assertHandler( \
                                        #condition, msg, GMX_CURRENT_FUNCTION, __FILE__, __LINE__)))
#    else
// Use an "immediately invoked function expression" to allow being
// used in constexpr context with older GCC versions
// https://akrzemi1.wordpress.com/2017/05/18/asserts-in-constexpr-functions/
#        define GMX_RELEASE_ASSERT(condition, msg)                                                         \
            ((void)((condition) ? (void)0 : [&]() {                                                        \
                ::gmx::internal::assertHandler(#condition, msg, GMX_CURRENT_FUNCTION, __FILE__, __LINE__); \
            }()))
#    endif
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
#    define GMX_ASSERT(condition, msg) ((void)0)
#else
#    define GMX_ASSERT(condition, msg) GMX_RELEASE_ASSERT(condition, msg)
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
[[noreturn]] void assertHandler(const char*                  condition,
                                const char*                  msg,
                                const char*                  func,
                                const std::filesystem::path& file,
                                int                          line);

} // namespace internal
//! \endcond

} // namespace gmx

#endif
