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

