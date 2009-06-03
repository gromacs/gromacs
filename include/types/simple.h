/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * 
 * And Hey:
 * GRoups of Organic Molecules in ACtion for Science
 */

#ifndef _simple_h
#define _simple_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef CPLUSPLUS
extern "C" {
#endif

#ifndef FALSE
#define FALSE   0
#endif
#ifndef TRUE
#define TRUE    1
#endif
#define BOOL_NR 2

#define XX	0			/* Defines for indexing in	*/
#define	YY	1			/* vectors			*/
#define ZZ	2
#define DIM   	3			/* Dimension of vectors		*/
#define XXXX    0                       /* defines to index matrices */
#define XXYY    1
#define XXZZ    2
#define YYXX    3
#define YYYY    4
#define YYZZ    5
#define ZZXX    6
#define ZZYY    7
#define ZZZZ    8
#ifndef HAVE_BOOL
#define bool int
  /* typedef int         	bool; */
#endif


typedef int     	atom_id;	/* To indicate an atoms id         */
#define NO_ATID		(atom_id)(~0)	/* Use this to indicate invalid atid */

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
typedef double   	real;
#define HAVE_REAL
#endif
#define GMX_MPI_REAL    MPI_DOUBLE
#define GMX_REAL_EPS    GMX_DOUBLE_EPS
#define GMX_REAL_MIN    GMX_DOUBLE_MIN
#define GMX_REAL_MAX    GMX_DOUBLE_MAX
#else
#ifndef HAVE_REAL
typedef float           real;
#define HAVE_REAL
#endif
#define GMX_MPI_REAL    MPI_FLOAT
#define GMX_REAL_EPS    GMX_FLOAT_EPS
#define GMX_REAL_MIN    GMX_FLOAT_MIN
#define GMX_REAL_MAX    GMX_FLOAT_MAX
#endif

#ifndef VECTORIZATION_BUFLENGTH
#define VECTORIZATION_BUFLENGTH 1000
  /* The total memory size of the vectorization buffers will
   * be 5*sizeof(real)*VECTORIZATION_BUFLENGTH
   */
#endif  
typedef real        	rvec[DIM];

typedef double       	dvec[DIM];

typedef real	    	matrix[DIM][DIM];

typedef real        	tensor[DIM][DIM];

typedef int             ivec[DIM];

typedef int             imatrix[DIM][DIM];

/* For the step count type gmx_step_t we aim for 8 bytes (64bit),
 * but we might only be able to get 4 bytes (32bit).
 */
#if (!(defined SIZEOF_LONG_LONG_INT) || SIZEOF_LONG_INT == 8)
typedef long int gmx_step_t;
#define gmx_step_fmt   "ld"
#define gmx_step_pfmt "%ld"
#define SIZEOF_GMX_STEP_T SIZEOF_LONG_INT
#else
typedef long long int gmx_step_t;
#define gmx_step_fmt   "lld"
#define gmx_step_pfmt "%lld"
#define SIZEOF_GMX_STEP_T SIZEOF_LONG_LONG_INT
#endif

#ifdef CPLUSPLUS
}
#endif

#endif

