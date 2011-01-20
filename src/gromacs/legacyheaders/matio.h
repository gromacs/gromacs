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

#ifndef _matio_h
#define _matio_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

gmx_bool matelmt_cmp(t_xpmelmt e1, t_xpmelmt e2);

t_matelmt searchcmap(int n,t_mapping map[],t_xpmelmt c);
/* Seach in the map for code 'c' and return entry number. 
 * return -1 if not found
 */

int getcmap(FILE *in,const char *fn,t_mapping **map);
/* Read the mapping table from in, return number of entries */

int readcmap(const char *fn,t_mapping **map);
/* Read the mapping table from fn, return number of entries */

void printcmap(FILE *out,int n,t_mapping map[]);
/* print mapping table to out */

void writecmap(const char *fn,int n,t_mapping map[]);
/* print mapping table to fn */

int read_xpm_matrix(const char *fnm, t_matrix **matrix);
/* Reads a number of matrices from .xpm file fnm and returns this number */

real **matrix2real(t_matrix *matrix,real **mat);
/* Converts an matrix in a t_matrix struct to a matrix of reals
 * When mat==NULL memory will be allocated 
 * Returns NULL when something went wrong
 */

void write_xpm_m(FILE *out, t_matrix m);
/* Writes a t_matrix struct to .xpm file */ 

void write_xpm3(FILE *out,unsigned int flags,
		       const char *title,const char *legend,
		       const char *label_x,const char *label_y,
		       int n_x,int n_y,real axis_x[],real axis_y[],
		       real *matrix[],real lo,real mid,real hi,
		       t_rgb rlo,t_rgb rmid,t_rgb rhi,int *nlevels);
/* See write_xpm.
 * Writes a colormap varying as rlo -> rmid -> rhi.
 */
void write_xpm_split(FILE *out,unsigned int flags,
			    const char *title,const char *legend,
			    const char *label_x,const char *label_y,
			    int n_x,int n_y,real axis_x[],real axis_y[],
			    real *matrix[],
			    real lo_top,real hi_top,int *nlevel_top,
			    t_rgb rlo_top,t_rgb rhi_top,
			    real lo_bot,real hi_bot,int *nlevel_bot,
			    gmx_bool bDiscreteColor,
			    t_rgb rlo_bot,t_rgb rhi_bot);
/* See write_xpm.
 * Writes a colormap with separate above and below diagonal colormaps.
 * If bDiscrete then a colormap with 16 fixed colors is used, first  of
 * which is white.
 */

void write_xpm(FILE *out,unsigned int flags,
		      const char *title,const char *legend,
		      const char *label_x,const char *label_y,
		      int n_x,int n_y,real t_x[],real t_y[],
		      real *matrix[],real lo,real hi,
		      t_rgb rlo,t_rgb rhi,int *nlevels);
/* out        xpm file
 * flags      flags, defined types/matrix.h
 *            MAT_SPATIAL_X
 *            MAT_SPATIAL_Y
 *            Defines if x and y are spatial dimensions,
 *            when not, there are n axis ticks at the middle of the elements,
 *            when set, there are n+1 axis ticks at the edges of the elements.
 * title      matrix title
 * legend     label for the continuous legend
 * label_x    label for the x-axis
 * label_y    label for the y-axis
 * n_x, n_y   size of the matrix
 * axis_x[]   the x-ticklabels (n_x or n_x+1)
 * axis_y[]   the y-ticklables (n_y or n_y+1)
 * *matrix[]  element x,y is matrix[x][y]
 * lo         output lower than lo is set to lo
 * hi         output higher than hi is set to hi
 * rlo        rgb value for level lo
 * rhi        rgb value for level hi
 * nlevels    number of color levels for the output
 */

real **mk_matrix(int nx, int ny, gmx_bool b1D);

void done_matrix(int nx, real ***m);

void clear_matrix(int nx, int ny, real **m);

#ifdef __cplusplus
}
#endif

#endif	/* _matio_h */
