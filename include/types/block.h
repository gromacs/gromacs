/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct {
  int multinr[MAXNODES];       	/* The indices for the multinode
                                 * version. For n=0, the blocks run from 0
                                 * upto multinr[index[0]]. The blocks for 
                                 * node n (n>0) run from 
                                 * index[multinr[n-1]] to index[multinr[n]].
                                 */
  int nr;			/* The number of blocks			*/
  atom_id *index;		/* Array of indices in a (dim: nr+1)	*/
  int nra;			/* The number of atoms 			*/
  atom_id *a;			/* Array of atom numbers in each group 	*/
				/* (dim: nra)				*/
				/* Block i (0<=i<nr) runs from		*/
				/* index[i] to index[i+1]-1. There will */
				/* allways be an extra entry in index	*/
				/* to terminate the table		*/
} t_block;

