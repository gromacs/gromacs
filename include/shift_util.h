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

#ifndef _shift_util_h
#define _shift_util_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "typedefs.h"
#include "gmxcomplex.h"
#include "physics.h"

/* Routines to set global constants for speeding up the calculation
 * of potentials, forces and gk-values.
 */
extern void set_shift_consts(FILE *log,real r1,real rc,rvec box,t_forcerec *fr);

extern real gk(real k,real rc,real r1);
/* Compute the Ghat function for a single k-value */

extern real gknew(real k,real rc,real r1);
/* Compute the (new!) Ghat function for a single k-value */

extern void pr_scalar_gk(char *fn,int nx,int ny,int nz,rvec box,real ***ghat);

extern real calc_dx2(rvec xi,rvec xj,rvec box);

extern void calc_dx(rvec xi,rvec xj,rvec box,rvec dx);

extern real phi_sr(FILE *log,int nj,rvec x[],real charge[],real rc,real r1,
		   rvec box,real phi[],t_block *excl,rvec f_sr[],bool bOld);

extern real shiftfunction(real r1,real rc,real R);

extern real spreadfunction(real r1,real rc,real R);

extern real potential(real r1,real rc,real R);

extern void calc_ener(FILE *fp,char *title,bool bHeader,
		      int nmol,int natoms,
		      real phi[],real charge[],t_block *excl);

extern real shift_LRcorrection(FILE *fp,t_nsborder *nsb,
			       t_commrec *cr,t_forcerec *fr,
			       real charge[],t_block *excl,rvec x[],
			       bool bOld,matrix box,matrix lrvir);
/* Calculate the self energy and forces
 * when using long range electrostatics methods.
 * Part of this is a constant, it is computed only once and stored in
 * a local variable. The remainder is computed every step.
 * PBC is taken into account. (Erik L.) 
 */

extern void calc_weights(int iatom,int nx,int ny,int nz,
			 rvec x,rvec box,rvec invh,ivec ixyz,real WXYZ[]);

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

extern void write_grid_pqr(char *fn,int nx,int ny,int nz,real ***phi);
extern void write_pqr(char *fn,t_atoms *atoms,rvec x[],real phi[],real dx);
/* Write a pdb file where the potential phi is printed as B-factor (for
 * viewing with rasmol). All atoms are moved over a distance dx in the X 
 * direction, to enable viewing of two data sets simultaneously with rasmol
 */

/******************************************************************
 *
 *   ROUTINES FOR GHAT MANIPULATION
 *
 ******************************************************************/
 
extern void symmetrize_ghat(int nx,int ny,int nz,real ***ghat);
/* Symmetrize the Ghat function. It is assumed that the 
 * first octant of the Ghat function is either read or generated
 * (all k-vectors from 0..nx/2 0..ny/2 0..nz/2).
 * Since Gk depends on the absolute value of k only, 
 * symmetry operations may shorten the time to generate it.
 */
 
extern void mk_ghat(FILE *fp,int nx,int ny,int nz,real ***ghat,
		    rvec box,real r1,real rc,bool bSym,bool bOld);
/* Generate a Ghat function from scratch. The ghat grid should
 * be allocated using the mk_rgrid function. When bSym, only
 * the first octant of the function is generated by direct calculation
 * and the above mentioned function is called for computing the rest.
 * When !bOld a new experimental function form will be used.
 */

extern real ***rd_ghat(FILE *log,char *fn,ivec igrid,rvec gridspacing,
		       rvec beta,int *porder,real *rshort,real *rlong);
/* Read a Ghat function from a file as generated by the program
 * mk_ghat. The grid size (number of grid points) is returned in
 * igrid, the space between grid points in gridspacing,
 * beta is a constant that determines the contribution of first
 * and second neighbours in the grid to the force
 * (See Luty et al. JCP 103 (1995) 3014)
 * porder determines whether 8 (when porder = 1) or 27 (when
 * porder = 2) neighbouring grid points are used for spreading
 * the charge.
 * rshort and rlong are the lengths used for generating the Ghat
 * function.
 */
		  
extern void wr_ghat(char *fn,int n1max,int n2max,int n3max,real h1,
		    real h2,real h3,real ***ghat,int nalias,
		    int porder,int niter,bool bSym,rvec beta,
		    real r1,real rc,real pval,real zval,real eref,real qopt);
/* Write a ghat file. (see above) */

extern void pr_scalar_gk(char *fn,int nx,int ny,int nz,rvec box,real ***ghat);

extern real analyse_diff(FILE *log,char *label,
			 int natom,rvec ffour[],rvec fpppm[],
			 real phi_f[],real phi_p[],real phi_sr[],
			 char *fcorr,char *pcorr,
			 char *ftotcorr,char *ptotcorr);
/* Analyse difference between forces from fourier (_f) and other (_p)
 * LR solvers (and potential also).
 * If the filenames are given, xvgr files are written.
 * returns the root mean square error in the force.
 */

#endif
