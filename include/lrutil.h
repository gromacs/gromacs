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
 * Great Red Owns Many ACres of Sand 
 */

#ifndef _lrutil_h
#define _lrutil_h

static char *SRCID_lrutil_h = "$Id$";

#include <math.h>
#include "typedefs.h"
#include "complex.h"

extern void set_LRconsts(FILE *log,real r1,real rc,rvec box,t_forcerec *fr);
/* Set constants necessary for Long Range electrostatics calculations */

extern real gk(real k,real rc,real r1);
/* Compute the Ghat function for a single k-value */

extern real calc_dx2(rvec xi,rvec xj,rvec box);

extern void calc_dx(rvec xi,rvec xj,rvec box,rvec dx);

extern real phi_sr(int nj,rvec x[],real charge[],real rc,real r1,
		   rvec box,real phi[],t_block *excl);

extern real shiftfunction(real r1,real rc,real R);

extern real spreadfunction(real r1,real rc,real R);

extern real potential(real r1,real rc,real R);

extern void calc_ener(FILE *fp,char *title,bool bHeader,
		      int nmol,int natoms,
		      real phi[],real charge[],t_block *excl);

extern real calc_selfenergy(FILE *fp,int natoms,real charge[],t_block *excl);
/* Calculate the self energy when using long range electrostatics methods.
 * Since this is a constant, it is computed only once and stored in
 * a local variable. On the next call this variable is returned straight 
 * away. No forces are computed so r1 must be > the distance between
 * excluded pairs.
 */

static void calc_lll(rvec box,rvec lll)
{
  lll[XX] = 2.0*M_PI/box[XX];
  lll[YY] = 2.0*M_PI/box[YY];
  lll[ZZ] = 2.0*M_PI/box[ZZ];
}

static void calc_k(rvec lll,int ix,int iy,int iz,int nx,int ny,int nz,rvec k)
{
#define IDX(i,n,x)  (i<=n/2) ? (i*x) : ((i-n)*x)
  k[XX] = IDX(ix,nx,lll[XX]);
  k[YY] = IDX(iy,ny,lll[YY]);
  k[ZZ] = IDX(iz,nz,lll[ZZ]);
#undef IDX
}

/******************************************************************
 *
 *   PLOTTING ROUTINES FOR DEBUGGING
 *
 ******************************************************************/
 
extern void plot_phi(char *fn,rvec box,int natoms,rvec x[],real phi[]);
/* Plot potential (or whatever) in a postscript matrix */

extern void print_phi(char *fn,int natoms,rvec x[],real phi[]);
/* Print to a text file in x y phi format */

extern void plot_qtab(char *fn,int nx,int ny,int nz,real ***qtab);
/* Plot a charge table to a postscript matrix */

extern void write_pqr(char *fn,t_atoms *atoms,rvec x[],real phi[],real dx);
/* Write a pdb file where the potential phi is printed as B-factor (for
 * viewing with rasmol). All atoms are moved over a distance dx in the X 
 * direction, to enable viewing of two data sets simultaneously with rasmol
 */
#endif
