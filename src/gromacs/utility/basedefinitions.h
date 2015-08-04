/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#ifndef _MSC_VER
#include <inttypes.h>
#endif

/*! \brief
 * Boolean type for use in \Gromacs C code.
 *
 * There is no standard size for 'bool' in C++, so when
 * we previously defined it to int for C code the data types
 * (and structs) would have different size depending on your compiler,
 * both at \Gromacs build time and when you use the library.
 * The only way around this is to NOT assume anything about the C++ type,
 * so we cannot use the name 'bool' in our C code anymore.
 */
typedef int gmx_bool;

#ifndef FALSE
/** False value for ::gmx_bool. */
#  define FALSE   0
#endif
#ifndef TRUE
/** True value for ::gmx_bool. */
#  define TRUE    1
#endif
/** Number of gmx_bool values. */
#define BOOL_NR 2

/*! \name Fixed-width integer types
 *
 * These types and macros provide the equivalent of 32- and 64-bit integer
 * types from C99 headers `stdint.h` and `inttypes.h`.  These headers are also
 * there in C++11.  The types and macros from here should be used instead of
 * `int32_t` etc.
 * MSVC doesn't support these before Visual Studio 2013.
 */
/*! \{ */
#ifdef _MSC_VER
typedef __int32 gmx_int32_t;
#define GMX_PRId32 "I32d"
#define GMX_SCNd32 "I32d"

typedef __int64 gmx_int64_t;
#define GMX_PRId64 "I64d"
#define GMX_SCNd64 "I64d"

typedef unsigned __int32 gmx_uint32_t;
#define GMX_PRIu32 "I32u"
#define GMX_SCNu32 "I32u"

typedef unsigned __int64 gmx_uint64_t;
#define GMX_PRIu64 "I64u"
#define GMX_SCNu64 "I64u"
#else
typedef int32_t gmx_int32_t;
#define GMX_PRId32 PRId32
#define GMX_SCNd32 SCNd32

typedef int64_t gmx_int64_t;
#define GMX_PRId64 PRId64
#define GMX_SCNd64 SCNd64

typedef uint32_t gmx_uint32_t;
#define GMX_PRIu32 PRIu32
#define GMX_SCNu32 SCNu32

typedef uint64_t gmx_uint64_t;
#define GMX_PRIu64 PRIu64
#define GMX_SCNu64 SCNu64
#endif

#define GMX_INT32_MAX INT32_MAX
#define GMX_INT32_MIN INT32_MIN

#define GMX_INT64_MAX INT64_MAX
#define GMX_INT64_MIN INT64_MIN

#define GMX_UINT32_MAX UINT32_MAX
#define GMX_UINT32_MIN UINT32_MIN

#define GMX_UINT64_MAX UINT64_MAX
#define GMX_UINT64_MIN UINT64_MIN
/*! \} */

/*! \def gmx_inline
 * \brief
 * Keyword to use in C code instead of C99 `inline`.
 *
 * Some of the C compilers we support do not recognize the C99 keyword
 * `inline`.  This macro should be used in C code and in shared C/C++ headers
 * to indicate a function is inlined.
 * C++ code should use plain `inline`, as that is already in C++98.
 */
#if !defined __cplusplus && defined _MSC_VER
#define gmx_inline __inline
#else
/* C++ or C99 */
#define gmx_inline inline
#endif

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

/*! \def gmx_cxx_const
 * \brief
 * Keyword to work around C/C++ differences in possible const keyword usage.
 *
 * Some functions that do not modify their input parameters cannot declare
 * those parameters as `const` and compile warning/error-free on both C and C++
 * compilers because of differences in `const` semantics.  This macro can be
 * used in cases where C++ allows `const`, but C does not like it, to make the
 * same declaration work for both.
 */
#ifdef __cplusplus
#define gmx_cxx_const const
#else
#define gmx_cxx_const
#endif

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
#ifdef __GNUC__
/* GCC, clang, and some ICC pretending to be GCC */
#  define gmx_unused __attribute__ ((unused))
#elif (defined(__INTEL_COMPILER) || defined(__ECC)) && !defined(_MSC_VER)
/* ICC on *nix */
#  define gmx_unused __attribute__ ((unused))
#elif defined(__PGI)
/* Portland group compilers */
#  define gmx_unused __attribute__ ((unused))
#elif defined _MSC_VER
/* MSVC */
#  define gmx_unused /*@unused@*/
#elif defined(__xlC__)
/* IBM */
#  define gmx_unused __attribute__ ((unused))
#else
#  define gmx_unused
#endif
#endif

#ifndef __has_feature
/** For compatibility with non-clang compilers. */
#define __has_feature(x) 0
#endif

/*! \def gmx_noreturn
 * \brief
 * Indicate that a function is not expected to return.
 */
#ifndef gmx_noreturn
#if defined(__GNUC__) || __has_feature(attribute_analyzer_noreturn)
#define gmx_noreturn __attribute__((noreturn))
#elif defined (_MSC_VER)
#define gmx_noreturn __declspec(noreturn)
#else
#define gmx_noreturn
#endif
#endif

/*! \def GMX_ALIGNMENT
 * \brief
 * Supports aligned variables */

/*! \def GMX_ALIGNED(type, alignment)
 * \brief
 * Declare variable with data alignment
 *
 * \param[in] type       Type of variable
 * \param[in] alignment  Alignment in multiples of type
 *
 * Typical usage:
 * \code
   GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH) buf[...];
   \endcode
 */
/* alignas(x) is not used even with GMX-CXX11 because it isn't in the list of
   tested features and thus might not be supported.
   MSVC before 2015 has align but doesn't support sizeof inside. */
#if defined(_MSC_VER) && (_MSC_VER > 1800 || defined(__ICL))
#  define GMX_ALIGNMENT 1
#  define GMX_ALIGNED(type, alignment) __declspec(align(alignment*sizeof(type))) type
#elif defined(__GNUC__) || defined(__clang__)
#  define GMX_ALIGNMENT 1
#  define GMX_ALIGNED(type, alignment) __attribute__ ((__aligned__(alignment*sizeof(type)))) type
#else
#  define GMX_ALIGNMENT 0
#  define GMX_ALIGNED(type, alignment)
#endif

/*! \brief
 * Macro to explicitly ignore an unused value.
 *
 * \ingroup module_utility
 */
#define GMX_UNUSED_VALUE(value) (void)value

#ifdef __cplusplus
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
template <typename T>
static inline void ignoreValueHelper(const T &)
{
}
//! \endcond
}   // namespace internal
}   // namespace gmx

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
#define GMX_IGNORE_RETURN_VALUE(call) \
        ::gmx::internal::ignoreValueHelper(call)
#endif

#endif
