/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gyas ROwers Mature At Cryogenic Speed
 */

#ifndef _poisson_h
#define _poisson_h

static char *SRCID_poisson_h = "$Id$";
#include "typedefs.h"

#define llim2 (-3)
#define ulim2  (3)

/* typedef for poisson solver */
typedef struct {
  int  nx,ny,nz;
  real ***ptr;
} t_PSgrid;

extern void unpack_PSgrid(t_PSgrid *grid,int *nx,int *ny,int *nz,real ****ptr);

extern void symmetrize_PSgrid(FILE *fp,t_PSgrid *grid,real sum);

extern void calc_nxyz(int nx,int ny,int nz,
		      int **nnx,int **nny,int **nnz);
/* Calculate tables to comput modulo (instead of function call) */
		      
extern real ps_gather_f(FILE *log,bool bVerbose,
			int natoms,rvec x[],rvec f[],real charge[],rvec box,
			real pot[],t_PSgrid *grid,rvec beta,t_nrnb *nrnb);

extern void spread_q_poisson(FILE *log,bool bVerbose,bool bCoulomb,
			     int natoms,rvec x[],real prop[],rvec box,
			     real rc,t_PSgrid *grid,t_nrnb *nrnb,
			     bool bOld,real r1);
/* Spreading charges (if BCoulomb)  or C6 (or any other property)
 * by convolution in real space of the spread function with   
 * the neighbouring grid points, (that is within the cut-off rc).
 * bOld and r1 are for backwards compatibility and testing.
 */

extern int solve_poisson(FILE *log,t_PSgrid *pot,t_PSgrid *rho,
			 bool bVerbose,t_nrnb *nrnb,int maxnit,real tol,
			 rvec box);
/* Solves a Poisson equation: Nabla^2 Phi = Rho, using Gauss-Seidel relaxation
 * returns the number of iterations.
 */

static void calc_invh_h(rvec box,int nx,int ny,int nz,rvec invh,rvec h)
{
  invh[XX] = nx/box[XX];
  invh[YY] = ny/box[YY];
  invh[ZZ] = nz/box[ZZ];
  h[XX]    = 1.0/invh[XX];
  h[YY]    = 1.0/invh[YY];
  h[ZZ]    = 1.0/invh[ZZ];
}

extern real do_poisson(FILE *log,       bool bVerbose,
		       t_inputrec *ir,  int natoms,
		       rvec x[],        rvec f[],
		       real charge[],   rvec box,
		       real phi[],      t_commrec *cr,
		       t_nrnb *nrnb,    int *nit,
		       bool bOld);
/* Calculate potentials etc. using a poisson solver */

extern real do_optimize_poisson(FILE *log,       bool bVerbose,
				t_inputrec *ir,  int natoms,
				rvec x[],        rvec f[],
				real charge[],   rvec box,
				real phi[],      t_commrec *cr,
				t_nrnb *nrnb,    rvec f_ref[],
				real phi_ref[],  rvec beta,
				bool bOld);
/* Does the same, but optimizes beta by comparison with something else */

#endif
