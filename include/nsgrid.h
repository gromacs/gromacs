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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"

extern t_grid *init_grid(FILE *fplog,t_forcerec *fr);

extern void done_grid(t_grid *grid);

extern void get_nsgrid_boundaries(t_grid *grid,matrix box,
				  int ncg,rvec *cgcm,
				  rvec grid_x0,rvec grid_x1);

extern void set_grid_ncg(t_grid *grid,int ncg);

extern void grid_first(FILE *log,t_grid *grid,gmx_domdec_t *dd,
		       int ePBC,matrix box,rvec izones_x0,rvec izones_x1,
		       real rlong,int ncg,rvec dens_x0,rvec dens_x1);

extern void fill_grid(FILE *log,
		      gmx_domdec_zones_t *dd_zones,
		      t_grid *grid,matrix box,
		      int cg0,int cg1,rvec cg_cm[]);

extern void calc_elemnr(FILE *log,t_grid *grid,int cg0,int cg1,int ncg);

extern void calc_ptrs(t_grid *grid);

extern void grid_last(FILE *log,t_grid *grid,int cg0,int cg1,int ncg);

extern int xyz2ci_(int nry,int nrz,int x,int y,int z);
#define xyz2ci(nry,nrz,x,y,z) ((nry)*(nrz)*(x)+(nrz)*(y)+(z))
/* Return the cell index */

extern void ci2xyz(t_grid *grid,int i,int *x,int *y,int *z);

extern void check_grid(FILE *log,t_grid *grid);

extern void print_grid(FILE *log,t_grid *grid);

extern void mv_grid(t_commrec *cr,t_grid *grid);
/* Move the grid over processors */


