/*
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
#ifndef _block_h
#define _block_h


#include "idef.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the block structure points into an array (usually of atom_ids).
   It is a list of starting indices for objects of consecutive ids, such
   as molecules. 
   For example, if this block denotes molecules, then the first molecule
   ranges from index[0] to index[1]-1 in the atom list.

   This makes the mapping from atoms to molecules O(Nmolecules) instead
   of O(Natoms) in size.  */
typedef struct {
  int nr;			/* The number of blocks			*/
  atom_id *index;		/* Array of indices (dim: nr+1) 	*/
  int nalloc_index;             /* The allocation size for index        */
} t_block;

typedef struct {
  int nr;			/* The number of blocks			*/
  atom_id *index;		/* Array of indices in a (dim: nr+1)	*/
  int nra;			/* The number of atoms 			*/
  atom_id *a;			/* Array of atom numbers in each group 	*/
				/* (dim: nra)				*/
				/* Block i (0<=i<nr) runs from		*/
				/* index[i] to index[i+1]-1. There will */
				/* allways be an extra entry in index	*/
				/* to terminate the table		*/
  int nalloc_index;             /* The allocation size for index        */
  int nalloc_a;                 /* The allocation size for a            */
} t_blocka;


#ifdef __cplusplus
}
#endif

#endif
