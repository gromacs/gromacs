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
 * Gromacs Runs On Most of All Computer Systems
 */

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

#define   GRID_STDDEV_FAC  sqrt(3)
#define NSGRID_STDDEV_FAC  2.0
/*
 * GRID_STDDEV_FAC * stddev is used to estimate the interaction density.
 * sqrt(3) gives a uniform load for a rectangular block of cg's.
 * For a sphere it is not a bad approximation for 4x1x1 up to 4x2x2.
 *
 * The extent of the neighborsearch grid is a bit larger than sqrt(3)
 * to account for less dense regions at the edges of the system.
 */

t_grid *init_grid(FILE *fplog,t_forcerec *fr);

void done_grid(t_grid *grid);

void get_nsgrid_boundaries(t_grid *grid,
				  gmx_domdec_t *dd,
				  matrix box,gmx_ddbox_t *ddbox,
				  rvec *gr0,rvec *gr1,
				  int ncg,rvec *cgcm,
				  rvec grid_x0,rvec grid_x1,
				  real *grid_density);
/* Return the ns grid boundaries grid_x0 and grid_x1
 * and the estimate for the grid density.
 * For non-bounded dimensions the boundaries are determined
 * from the average and std.dev. of cgcm.
 * The are determined from box, unless gr0!=NULL or gr1!=NULL,
 * then they are taken from gr0 or gr1.
 * With dd and unbounded dimensions, the proper grid borders for cells
 * on the edges are determined from cgcm.
 */

void grid_first(FILE *log,t_grid *grid,
		       gmx_domdec_t *dd,const gmx_ddbox_t *ddbox,
		       int ePBC,matrix box,rvec izones_x0,rvec izones_x1,
		       real rlong,real grid_density);

void fill_grid(FILE *log,
		      gmx_domdec_zones_t *dd_zones,
		      t_grid *grid,int ncg_tot,
		      int cg0,int cg1,rvec cg_cm[]);
/* Allocates space on the grid for ncg_tot cg's.
 * Fills the grid with cg's from cg0 to cg1.
 * When cg0 is -1, contiues filling from grid->nr to cg1.
 */

void calc_elemnr(FILE *log,t_grid *grid,int cg0,int cg1,int ncg);

void calc_ptrs(t_grid *grid);

void grid_last(FILE *log,t_grid *grid,int cg0,int cg1,int ncg);

int xyz2ci_(int nry,int nrz,int x,int y,int z);
#define xyz2ci(nry,nrz,x,y,z) ((nry)*(nrz)*(x)+(nrz)*(y)+(z))
/* Return the cell index */

void ci2xyz(t_grid *grid,int i,int *x,int *y,int *z);

void check_grid(FILE *log,t_grid *grid);

void print_grid(FILE *log,t_grid *grid);

void mv_grid(t_commrec *cr,t_grid *grid);
/* Move the grid over processors */

#ifdef __cplusplus
}
#endif


