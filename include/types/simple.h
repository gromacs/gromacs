/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

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
#define MAXPROC 	256		/* Max number of processors 	*/

#ifndef HAVE_BOOL
#define bool int
  /* typedef int         	bool; */
#endif


typedef unsigned int 	atom_id;	/* To indicate an atoms id         */
typedef unsigned int 	atom_nr;	/* To indicate a number of atoms   */
/* Notice that atom_id cannot use its full range when the types of atom_id */
/* and atom_nr are equal.                                                  */
#define NO_ATID		(atom_id)(~0)	/* Use this to indicate invalid atid */

#ifdef DOUBLE
typedef double   	real;
#define GMX_MPI_REAL    MPI_DOUBLE
#define GMX_REAL_EPS    2.2e-16
#else
typedef float           real;
#define GMX_MPI_REAL    MPI_FLOAT
#define GMX_REAL_EPS    1.2e-07
#endif

typedef real        	rvec[DIM];

typedef real	    	matrix[DIM][DIM];

typedef real        	tensor[DIM][DIM];

typedef int             ivec[DIM];

typedef int             imatrix[DIM][DIM];

#ifndef HAVE_ULONG
typedef unsigned long   ulong;
#endif

#ifdef NEED_USHORT
typedef unsigned short  ushort;
#endif

#ifdef CPLUSPLUS
}
#endif


