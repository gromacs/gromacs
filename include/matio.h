/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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

#ifndef _matio_h
#define _matio_h

static char *SRCID_matio_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) matio.h 1.11 5/20/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"

extern bool matelmt_cmp(t_xpmelmt e1, t_xpmelmt e2);

extern t_matelmt searchcmap(int n,t_mapping map[],t_xpmelmt c);
/* Seach in the map for code 'c' and return entry number. 
 * return -1 if not found
 */

extern int getcmap(FILE *in,char *fn,t_mapping **map);
/* Read the mapping table from in, return number of entries */

extern int readcmap(char *fn,t_mapping **map);
/* Read the mapping table from fn, return number of entries */

extern void printcmap(FILE *out,int n,t_mapping map[]);
/* print mapping table to out */

extern void writecmap(char *fn,int n,t_mapping map[]);
/* print mapping table to fn */

extern int read_xpm_matrix(char *fnm, t_matrix **matrix);
/* Reads a number of matrices from .xpm file fnm and returns this number */

extern int matrix2real(t_matrix *matrix, real ***mat);
/* Converts an integer matrix in a t_matrix struct to a matrix of reals
 * mat is snewed by matrix2real
 * Returns 0 when something went wrong
 */

extern void write_xpm_m(FILE *out, t_matrix m);
/* Writes a t_matrix struct to .xpm file */ 

extern void write_xpm3(FILE *out,
		       char *title,char *legend,char *label_x,char *label_y,
		       int n_x,int n_y,real axis_x[],real axis_y[],
		       real *matrix[],real lo,real mid,real hi,
		       t_rgb rlo,t_rgb rmid,t_rgb rhi,int *nlevels);
/* See write_xpm.
 * Writes a colormap varying as rlo -> rmid -> rhi.
 */

extern void write_xpm(FILE *out,
		      char *title,char *legend,char *label_x,char *label_y,
		      int n_x,int n_y,real t_x[],real t_y[],
		      real *matrix[],real lo,real hi,
		      t_rgb rlo,t_rgb rhi,int *nlevels);
/* out        xpm file
 * title      matrix title
 * legend     label for the continuous legend
 * label_x    label for the x-axis
 * label_y    label for the y-axis
 * n_x, n_y   size of the matrix
 * axis_x[]   the x-ticklabels
 * axis_y[]   the y-ticklables
 * *matrix[]  element x,y is matrix[x][y]
 * lo         output lower than lo is set to lo
 * hi         output higher than hi is set to hi
 * rlo        rgb value for level lo
 * rhi        rgb value for level hi
 * nlevels    number of levels for the output minus one
 */

#endif	/* _matio_h */
