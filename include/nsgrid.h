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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _nsgrid_h
#define _nsgrid_h

static char *SRCID_nsgrid_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) nsgrid.h 1.7 11/23/92"
#endif /* HAVE_IDENT */

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


