/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */
typedef struct {
  int multinr[MAXPROC];		/* The indices for the multiprocessor 
                                 * version. For n=0, the blocks run from 0
                                 * upto multinr[index[0]]. The blocks for 
                                 * processor n (n>0) run from 
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

