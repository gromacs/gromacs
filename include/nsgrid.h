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

#ifndef _nsgrid_h
#define _nsgrid_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"

extern void init_grid(FILE *log,t_grid *grid,
		      int delta,matrix box,real rlong,int ncg);

extern void done_grid(t_grid *grid);

extern void grid_first(FILE *log,t_grid *grid,matrix box,real rlong);

extern void fill_grid(FILE *log,bool bDD,int cg_index[],
		      t_grid *grid,matrix box,
		      int ncg,int cg0,int cg1,rvec cg_cm[]);

extern void calc_elemnr(FILE *log,bool bDD,int cg_index[],
			t_grid *grid,int cg0,int cg1,int ncg);

extern void calc_ptrs(t_grid *grid);

extern void grid_last(FILE *log,bool bDD,int cg_index[],
		      t_grid *grid,int cg0,int cg1,int ncg);

extern int xyz2ci_(int nry,int nrz,int x,int y,int z);
#define xyz2ci(nry,nrz,x,y,z) ((nry)*(nrz)*(x)+(nrz)*(y)+(z))
/* Return the cell index */

extern void ci2xyz(t_grid *grid,int i,int *x,int *y,int *z);

extern void check_grid(FILE *log,t_grid *grid);

extern void print_grid(FILE *log,t_grid *grid,bool bDD,int cg_index[]);

extern void mv_grid(t_commrec *cr,bool bDD,int cg_index[],
		    t_grid *grid,int cgload[]);
/* Move the grid over processors */

#endif	/* ns_grid_h */


