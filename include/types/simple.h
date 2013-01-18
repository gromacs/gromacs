/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

/* Dont remove this instance of HAVE_CONFIG_H!!!
 *
 * We dont _require_ config.h here, but IF one is
 * available it might contain valuable information about simple types
 * that helps us automate things better and avoid bailing out.
 *
 * Note that this does not have to be the gromacs config.h - several
 * package setups define these simple types.
 */
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

/* Information about integer data type sizes */
#include <limits.h>

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
#define GMX_DOUBLE_EPS   1.11022302E-16

/*! \brief Maximum double precision value - reduced 1 unit in last digit for MSVC */
#define GMX_DOUBLE_MAX   1.79769312E+308

/*! \brief Minimum double precision value */
#define GMX_DOUBLE_MIN   2.22507386E-308

/*! \brief Single precision accuracy */
#define GMX_FLOAT_EPS    5.96046448E-08

/*! \brief Maximum single precision value - reduced 1 unit in last digit for MSVC */
#define GMX_FLOAT_MAX    3.40282346E+38

/*! \brief Minimum single precision value */
#define GMX_FLOAT_MIN    1.17549435E-38


/* Check whether we already have a real type! */
#ifdef GMX_DOUBLE

#ifndef HAVE_REAL
typedef double      real;
#define HAVE_REAL
#endif

#define GMX_MPI_REAL    MPI_DOUBLE
#define GMX_REAL_EPS    GMX_DOUBLE_EPS
#define GMX_REAL_MIN    GMX_DOUBLE_MIN
#define GMX_REAL_MAX    GMX_DOUBLE_MAX
#define gmx_real_fullprecision_pfmt "%21.14e"
#else

#ifndef HAVE_REAL
typedef float           real;
#define HAVE_REAL
#endif

#define GMX_MPI_REAL    MPI_FLOAT
#define GMX_REAL_EPS    GMX_FLOAT_EPS
#define GMX_REAL_MIN    GMX_FLOAT_MIN
#define GMX_REAL_MAX    GMX_FLOAT_MAX
#define gmx_real_fullprecision_pfmt "%14.7e"
#endif

typedef real            rvec[DIM];

typedef double          dvec[DIM];

typedef real            matrix[DIM][DIM];

typedef real            tensor[DIM][DIM];

typedef int             ivec[DIM];

typedef int             imatrix[DIM][DIM];


/* For the step count type gmx_large_int_t we aim for 8 bytes (64bit),
 * but we might only be able to get 4 bytes (32bit).
 *
 * We first try to find a type without reyling on any SIZEOF_XXX defines.
 *
 * Avoid using "long int" if we can. This type is really dangerous,
 * since the width frequently depends on compiler options, and they
 * might not be set correctly when (buggy) Cmake is detecting things.
 * Instead, start by looking for "long long", and just go down if we
 * have to (rarely on new systems). /EL 20100810
 */
#if ( (defined SIZEOF_LONG_LONG_INT && SIZEOF_LONG_LONG_INT == 8) || (defined LLONG_MAX && LLONG_MAX == 9223372036854775807LL) )

/* Long long int is 64 bit */
typedef long long int gmx_large_int_t;
#define gmx_large_int_fmt   "lld"
#define gmx_large_int_pfmt "%lld"
#define SIZEOF_GMX_LARGE_INT  8
#define GMX_LARGE_INT_MAX     9223372036854775807LL
#define GMX_LARGE_INT_MIN     (-GMX_LARGE_INT_MAX - 1LL)
#define GMX_MPI_LARGE_INT MPI_LONG_LONG_INT

#elif ( (defined SIZEOF_LONG_INT && SIZEOF_LONG_INT == 8) || (defined LONG_MAX && LONG_MAX == 9223372036854775807L) )

/* Long int is 64 bit */
typedef long int gmx_large_int_t;
#define gmx_large_int_fmt   "ld"
#define gmx_large_int_pfmt "%ld"
#define SIZEOF_GMX_LARGE_INT  8
#define GMX_LARGE_INT_MAX     9223372036854775807LL
#define GMX_LARGE_INT_MIN     (-GMX_LARGE_INT_MAX - 1LL)
#define GMX_MPI_LARGE_INT MPI_LONG_INT

#elif ( (defined SIZEOF_INT && SIZEOF_INT == 8) || (defined INT_MAX && INT_MAX == 9223372036854775807L) )

/* int is 64 bit */
typedef int gmx_large_int_t;
#define gmx_large_int_fmt   "d"
#define gmx_large_int_pfmt  "%d"
#define SIZEOF_GMX_LARGE_INT  8
#define GMX_LARGE_INT_MAX     9223372036854775807LL
#define GMX_LARGE_INT_MIN     (-GMX_LARGE_INT_MAX - 1LL)
#define GMX_MPI_LARGE_INT MPI_INT

#elif ( (defined INT_MAX && INT_MAX == 2147483647) || (defined SIZEOF_INT && SIZEOF_INT == 4) )

/* None of the above worked, try a 32 bit integer */
typedef int gmx_large_int_t;
#define gmx_large_int_fmt   "d"
#define gmx_large_int_pfmt "%d"
#define SIZEOF_GMX_LARGE_INT  4
#define GMX_LARGE_INT_MAX     2147483647
#define GMX_LARGE_INT_MIN     (-GMX_LARGE_INT_MAX - 1)
#define GMX_MPI_LARGE_INT MPI_INT

#else

#error "Cannot find any 32 or 64 bit integer data type. Please extend the gromacs simple.h file!"

#endif


#ifndef gmx_inline
/* config.h tests for inline definitions and should work on a much wider range
 * of compilers, but does not work with installed headers. These compiler checks
 * still enable a real inline keyword for the most common compilers.
 */

/* Try to define suitable inline keyword for gmx_inline.
 * Set it to empty if we cannot find one (and dont complain to the user)
 */
#ifndef __cplusplus

#ifdef __GNUC__
/* GCC */
#  define gmx_inline   __inline__
#elif (defined(__INTEL_COMPILER) || defined(__ECC)) && defined(__ia64__)
/* ICC */
#  define gmx_inline __inline__
#elif defined(__PATHSCALE__)
/* Pathscale */
#  define gmx_inline __inline__
#elif defined(__PGIC__)
/* Portland */
#  define gmx_inline __inline
#elif defined _MSC_VER
/* MSVC */
#  define gmx_inline __inline
#elif defined(__xlC__)
/* IBM */
#  define gmx_inline __inline
#else
#  define gmx_inline
#endif

#else
/* C++ */
#  define gmx_inline inline
#endif

#endif /* ifndef gmx_inline */


/* Restrict keywords. Note that this has to be done for C++ too, unless
 * it was set from the more general checks if we had config.h (gmx internal)
 */
#ifndef gmx_restrict

#ifdef __GNUC__
/* GCC */
#  define gmx_restrict   __restrict__
#elif (defined(__INTEL_COMPILER) || defined(__ECC)) && defined(__ia64__)
/* ICC */
#  define gmx_restrict __restrict__
#elif defined(__PATHSCALE__)
/* Pathscale */
#  define gmx_restrict __restrict
#elif defined(__PGIC__)
/* Portland */
#  define gmx_restrict __restrict
#elif defined _MSC_VER
/* MSVC */
#  define gmx_restrict __restrict
#elif defined(__xlC__)
/* IBM */
#  define gmx_restrict __restrict
#else
#  define gmx_restrict
#endif

#endif

/* Standard sizes for char* string buffers */
#define STRLEN 4096
#define BIG_STRLEN 1048576


#ifdef __cplusplus
}
#endif

#endif
