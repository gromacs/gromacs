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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct {
  int	 nr;		/* Total number of charge groups	*/
  int	 nrx,nry,nrz;	/* The dimension of the grid		*/
  int    ncells;	/* Total number of cells		*/
  int    maxcells;	/* Max number of cells (when box grows)	*/
  int	 delta;		/* The distance in cells to be searched */
  int    gmax;		/* The size of the largest grid cell	*/
  int	 *cell_index;	/* The cell number of each cg		*/
  int    *index;	/* The index into a for each cell	*/
			/* The location of the cell in the index*/
			/* array can be found by calling xyz2ci	*/
  int    *nra;		/* The number of entries in a cell	*/
  int    *a;		/* The grid of cgs			*/
} t_grid;

