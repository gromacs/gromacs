/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
 * Basic types and macros used throughout \Gromacs.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_BASEDEFINITIONS_H
#define GMX_UTILITY_BASEDEFINITIONS_H

#include <stdint.h>

#include <cinttypes>
#include <cstddef>

//! Identical to bool
typedef bool gmx_bool;

#ifndef FALSE
/** False value for ::gmx_bool. */
#    define FALSE false
#endif
#ifndef TRUE
/** True value for ::gmx_bool. */
#    define TRUE true
#endif
/** Number of gmx_bool values. */
#define BOOL_NR 2

namespace gmx
{
/*! \brief Integer type for indexing into arrays or vectors
 *
 * Same as ptrdiff_t.
 */
using index = std::ptrdiff_t;

//! Return signed size of container
template<typename T>
index ssize(const T& t)
{
    return t.size();
}
} // namespace gmx

/* ICC, GCC, MSVC, Pathscale, PGI, XLC support __restrict.
 * Any other compiler can be added here. */
/*! \brief
 * Keyword to use in instead of C99 `restrict`.
 *
 * We cannot use `restrict` because it is only in C99, but not in C++.
 * This macro should instead be used to allow easily supporting different
 * compilers.
 */
#define gmx_restrict __restrict

/*! \def gmx_unused
 * \brief
 * Attribute to suppress compiler warnings about unused function parameters.
 *
 * This attribute suppresses compiler warnings about unused function arguments
 * by marking them as possibly unused.  Some arguments are unused but
 * have to be retained to preserve a function signature
 * that must match that of another function.
 * Some arguments are only used in *some* conditional compilation code paths
 * (e.g. MPI).
 */
#ifndef gmx_unused
#    ifdef __GNUC__
/* GCC, clang, and some ICC pretending to be GCC */
#        define gmx_unused __attribute__((unused))
#    elif (defined(__INTEL_COMPILER) || defined(__ECC)) && !defined(_MSC_VER)
/* ICC on *nix */
#        define gmx_unused __attribute__((unused))
#    elif defined(__PGI)
/* Portland group compilers */
#        define gmx_unused __attribute__((unused))
#    elif defined _MSC_VER
/* MSVC */
#        define gmx_unused /*@unused@*/
#    elif defined(__xlC__)
/* IBM */
#        define gmx_unused __attribute__((unused))
#    else
#        define gmx_unused
#    endif
#endif

/*! \brief Attribute to explicitly indicate that a parameter or
 * locally scoped variable is used just in debug mode.
 *
 * \ingroup module_utility
 */
#ifdef NDEBUG
#    define gmx_used_in_debug gmx_unused
#else
#    define gmx_used_in_debug
#endif

#ifndef __has_feature
/** For compatibility with non-clang compilers. */
#    define __has_feature(x) 0
#endif

/*! \brief
 * Macro to explicitly ignore an unused value.
 *
 * \ingroup module_utility
 *
 * \todo Deprecated - use gmx_unused
 */
#define GMX_UNUSED_VALUE(value) (void)value

#ifdef __clang__
#    define DO_PRAGMA(x) _Pragma(#    x)
#    define CLANG_DIAGNOSTIC_IGNORE(warning) \
        _Pragma("clang diagnostic push") DO_PRAGMA(clang diagnostic ignored #warning)
#    define DIAGNOSTIC_RESET _Pragma("clang diagnostic pop")
#else
//! Ignore specified clang warning until DIAGNOSTIC_RESET
#    define CLANG_DIAGNOSTIC_IGNORE(warning)
//! Reset all diagnostics to default
#    define DIAGNOSTIC_RESET
#endif

namespace gmx
{
namespace internal
{
/*! \cond internal */
/*! \internal \brief
 * Helper for ignoring values in macros.
 *
 * \ingroup module_utility
 */
template<typename T>
static inline void ignoreValueHelper(const T& /*unused*/)
{
}
//! \endcond
} // namespace internal
} // namespace gmx

/*! \brief
 * Macro to explicitly ignore a return value of a call.
 *
 * Mainly meant for ignoring values of functions declared with
 * `__attribute__((warn_unused_return))`.  Makes it easy to find those places if
 * they need to be fixed, and document the intent in cases where the return
 * value really can be ignored.  It also makes it easy to adapt the approach so
 * that they don't produce warnings.  A cast to void doesn't remove the warning
 * in gcc, while adding a dummy variable can cause warnings about an unused
 * variable.
 *
 * \ingroup module_utility
 */
#define GMX_IGNORE_RETURN_VALUE(call) ::gmx::internal::ignoreValueHelper(call)

#endif
