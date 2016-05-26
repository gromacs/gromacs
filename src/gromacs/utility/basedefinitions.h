/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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
 *
 * (MSVC 2015 still doesn't support the format strings.)
 */
/*! \{ */
typedef int32_t gmx_int32_t;
typedef int64_t gmx_int64_t;
typedef uint32_t gmx_uint32_t;
typedef uint64_t gmx_uint64_t;

#ifdef _MSC_VER
#define GMX_PRId32 "I32d"
#define GMX_SCNd32 "I32d"

#define GMX_PRId64 "I64d"
#define GMX_SCNd64 "I64d"

#define GMX_PRIu32 "I32u"
#define GMX_SCNu32 "I32u"

#define GMX_PRIu64 "I64u"
#define GMX_SCNu64 "I64u"
#else
#define GMX_PRId32 PRId32
#define GMX_SCNd32 SCNd32

#define GMX_PRId64 PRId64
#define GMX_SCNd64 SCNd64

#define GMX_PRIu32 PRIu32
#define GMX_SCNu32 SCNu32

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

/*! \def GMX_CXX11_COMPILATION
 * \brief
 * Defined to 1 when compiling as C++11.
 *
 * While \Gromacs only supports C++11 compilation, there are some parts of the
 * code that are compiled with other tools than the actual C++ compiler, and
 * these may not support C++11.  Most notable such case is all of CUDA code
 * (with CUDA versions older than 6.5), but other types of kernels might also
 * have similar limitations in the future.
 *
 * The define is intended for conditional compilation in low-level headers that
 * need to support inclusion from such non-C++11 files, but get significant
 * benefit (e.g., for correctness checking or more convenient use) from C++11.
 * It should only be used for features that do not influence the ABI of the
 * header; e.g., static_asserts or additional helper methods.
 */
#if defined __cplusplus && __cplusplus >= 201103L
#    define GMX_CXX11_COMPILATION 1
#else
#    define GMX_CXX11_COMPILATION 0
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

/*! \def gmx_constexpr
 * \brief C++11 constexpr everywhere except MSVC 2013, where it is empty.
 *
 * Support for constexpr was not added until MSVC 2015, and it still
 * seems to be unreliable. Since interacting with parts of libc++ and
 * libstdc++ depend on it (for instance the min/max calls in our
 * random engines), we need to specify it for other compilers.
 */
#ifndef gmx_constexpr
#if !defined(_MSC_VER)
#    define gmx_constexpr constexpr
#else
#    define gmx_constexpr
#endif
#endif

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

#if (defined(__GNUC__) && !defined(__clang__)) || defined(__ibmxl__) || defined(__xlC__) || defined(__PATHCC__)
// Gcc-4.6.4 does not support alignas, but both gcc, pathscale and xlc
// support the standard GNU alignment attributes. PGI also sets __GNUC__ now,
// and mostly supports it. clang 3.2 does not support the GCC alignment attribute.
#    define GMX_ALIGNED(type, alignment) __attribute__ ((aligned(alignment*sizeof(type)))) type
#else
// If nothing else works we rely on C++11. This will for instance work for MSVC2015 and later.
// If you get an error here, find out what attribute to use to get your compiler to align
// data properly and add it as a case.
#    define GMX_ALIGNED(type, alignment) alignas(alignment*alignof(type)) type
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
