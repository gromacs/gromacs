/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#ifndef _simple_h
#define _simple_h

/* Information about integer data type sizes */
#include <limits.h>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#ifndef _MSC_VER
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif


#define XX      0           /* Defines for indexing in	*/
#define YY      1           /* vectors			*/
#define ZZ      2
#define DIM     3           /* Dimension of vectors		*/
#define XXXX    0           /* defines to index matrices */
#define XXYY    1
#define XXZZ    2
#define YYXX    3
#define YYYY    4
#define YYZZ    5
#define ZZXX    6
#define ZZYY    7
#define ZZZZ    8

/* There is no standard size for 'bool' in C++, so when
 * we previously defined it to int for C code the data types
 * (and structs) would have different size depending on your compiler,
 * both at gromacs build time and when you use the library.
 * The only way around this is to NOT assume anything about the C++ type,
 * so we cannot use the name 'bool' in our C code anymore.
 */

typedef int gmx_bool;

#ifndef FALSE
#  define FALSE   0
#endif
#ifndef TRUE
#  define TRUE    1
#endif
#define BOOL_NR 2


typedef int         atom_id;      /* To indicate an atoms id         */
#define NO_ATID     (atom_id)(~0) /* Use this to indicate invalid atid */

/*! \brief Double precision accuracy */
#define GMX_DOUBLE_EPS   2.2204460492503131e-16

/*! \brief Maximum double precision value - reduced 1 unit in last digit for MSVC */
#define GMX_DOUBLE_MAX   1.7976931348623157e+308

/*! \brief Minimum double precision value */
#define GMX_DOUBLE_MIN   2.2250738585072014e-308

/*! \brief Single precision accuracy */
#define GMX_FLOAT_EPS    1.19209290e-07F

/*! \brief Maximum single precision value - reduced 1 unit in last digit for MSVC */
#define GMX_FLOAT_MAX    3.40282346E+38F

/*! \brief Minimum single precision value */
#define GMX_FLOAT_MIN    1.175494351E-38F

#ifdef __PGI
/* The portland group x86 C/C++ compilers do not treat negative zero initializers
 * correctly, but "optimizes" them to positive zero, so we implement it explicitly.
 * These constructs are optimized to simple loads at compile time. If you want to
 * use them on other compilers those have to support gcc preprocessor extensions.
 * Note: These initializers might be sensitive to the endianness (which can
 * be different for byte and word order), so check that it works for your platform
 * and add a separate section if necessary before adding to the ifdef above.
 */
#    define GMX_DOUBLE_NEGZERO  ({ const union { int  di[2]; double d; } _gmx_dzero = {0, -2147483648}; _gmx_dzero.d; })
#    define GMX_FLOAT_NEGZERO   ({ const union { int  fi; float f; } _gmx_fzero = {-2147483648}; _gmx_fzero.f; })
#else
/*! \brief Negative zero in double */
#    define GMX_DOUBLE_NEGZERO  (-0.0)

/*! \brief Negative zero in float */
#    define GMX_FLOAT_NEGZERO   (-0.0f)
#endif


/* Check whether we already have a real type! */
#ifdef GMX_DOUBLE

#ifndef HAVE_REAL
typedef double      real;
#define HAVE_REAL
#endif

#define GMX_MPI_REAL      MPI_DOUBLE
#define GMX_REAL_EPS      GMX_DOUBLE_EPS
#define GMX_REAL_MIN      GMX_DOUBLE_MIN
#define GMX_REAL_MAX      GMX_DOUBLE_MAX
#define GMX_REAL_NEGZERO  GMX_DOUBLE_NEGZERO
#define gmx_real_fullprecision_pfmt "%21.14e"
#else

#ifndef HAVE_REAL
typedef float           real;
#define HAVE_REAL
#endif

#define GMX_MPI_REAL      MPI_FLOAT
#define GMX_REAL_EPS      GMX_FLOAT_EPS
#define GMX_REAL_MIN      GMX_FLOAT_MIN
#define GMX_REAL_MAX      GMX_FLOAT_MAX
#define GMX_REAL_NEGZERO  GMX_FLOAT_NEGZERO
#define gmx_real_fullprecision_pfmt "%14.7e"
#endif

typedef real            rvec[DIM];

typedef double          dvec[DIM];

typedef real            matrix[DIM][DIM];

typedef real            tensor[DIM][DIM];

typedef int             ivec[DIM];

typedef int             imatrix[DIM][DIM];

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

#if !defined __cplusplus && _MSC_VER
#define gmx_inline __inline
#else
/* C++ or C99 */
#define gmx_inline inline
#endif

/* ICC, GCC, MSVC, Pathscale, PGI, XLC support __restrict.
 * Any other compiler can be added here. We cannot
 * use restrict because it is in C99 but not in C++ */
#define gmx_restrict __restrict

/*
 * These attributes suppress compiler warnings about unused function arguments
 * by marking them as possibly unused. Some arguments are unused but
 * have to be retained to preserve a function signature
 * that must match that of another function.
 * Some arguments are only used in *some* code paths (e.g. MPI)
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

/* Standard sizes for char* string buffers */
#define STRLEN 4096
#define BIG_STRLEN 1048576


#ifdef __cplusplus
}
#endif

#endif
