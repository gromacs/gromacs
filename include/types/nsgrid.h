/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
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

