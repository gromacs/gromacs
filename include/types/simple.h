/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
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
/* Max number of nodes 	*/  
#define MAXNODES	256	
#ifndef HAVE_BOOL
#define bool int
  /* typedef int         	bool; */
#endif


typedef int     	atom_id;	/* To indicate an atoms id         */
#define NO_ATID		(atom_id)(~0)	/* Use this to indicate invalid atid */

  /* Check whether we already have a real type! */
#ifdef DOUBLE
#ifndef HAVE_REAL
typedef double   	real;
#define HAVE_REAL
#endif
#define GMX_MPI_REAL    MPI_DOUBLE
#define GMX_REAL_EPS    2.2e-16
#else
#ifndef HAVE_REAL
typedef float           real;
#define HAVE_REAL
#endif
#define GMX_MPI_REAL    MPI_FLOAT
#define GMX_REAL_EPS    1.2e-07
#endif

#ifndef VECTORIZATION_BUFLENGTH
#define VECTORIZATION_BUFLENGTH 1000
  /* The total memory size of the vectorization buffers will
   * be 5*sizeof(real)*VECTORIZATION_BUFLENGTH
   */
#endif  
typedef real        	rvec[DIM];

typedef real	    	matrix[DIM][DIM];

typedef real        	tensor[DIM][DIM];

typedef int             ivec[DIM];

typedef int             imatrix[DIM][DIM];

#ifdef CPLUSPLUS
}
#endif

#endif

