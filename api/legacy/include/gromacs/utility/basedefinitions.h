/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
 * Basic types and macros used throughout \Gromacs.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_BASEDEFINITIONS_H
#define GMX_UTILITY_BASEDEFINITIONS_H

#include <cinttypes>
#include <cstddef>
#include <cstdint>

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
using Index = std::ptrdiff_t;

//! Return signed size of container
template<typename T>
Index ssize(const T& t)
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
/* GCC, clang, and any pretending to be or based on them */
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

#if defined(__GNUC__) && !defined(__clang__)
#    define DO_PRAGMA(x) _Pragma(#    x)
#    define GCC_DIAGNOSTIC_IGNORE(warning) \
        _Pragma("GCC diagnostic push") DO_PRAGMA(GCC diagnostic ignored warning)
#    define GCC_DIAGNOSTIC_RESET _Pragma("GCC diagnostic pop")
#else
//! Ignore specified clang warning until GCC_DIAGNOSTIC_RESET
#    define GCC_DIAGNOSTIC_IGNORE(warning)
//! Reset all diagnostics to default
#    define GCC_DIAGNOSTIC_RESET
#endif

#if defined(__clang__) && !defined(DO_PRAGMA)
#    define DO_PRAGMA(x) _Pragma(#    x)
#    define CLANG_DIAGNOSTIC_IGNORE(warning) \
        _Pragma("clang diagnostic push") DO_PRAGMA(clang diagnostic ignored warning)
#    define CLANG_DIAGNOSTIC_RESET _Pragma("clang diagnostic pop")
#else
//! Ignore specified clang warning until CLANG_DIAGNOSTIC_RESET
#    define CLANG_DIAGNOSTIC_IGNORE(warning)
//! Reset all diagnostics to default
#    define CLANG_DIAGNOSTIC_RESET
#endif

#ifdef _MSC_VER
#    define MSVC_DIAGNOSTIC_IGNORE(id) __pragma(warning(push)) __pragma(warning(disable : id))
#    define MSVC_DIAGNOSTIC_RESET __pragma(warning(pop))
#else
//! Ignore specified MSVC warning until MSVC_DIAGNOSTIC_RESET
#    define MSVC_DIAGNOSTIC_IGNORE(warning)
//! Reset all diagnostics to default
#    define MSVC_DIAGNOSTIC_RESET
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
