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
static char *SRCID_pssolve_c = "$Id$";
#include "poisson.h"
#include "physics.h"
#include "vec.h"
	
int solve_poisson(FILE *log,t_PSgrid *pot,t_PSgrid *rho,
		   bool bVerbose,t_nrnb *nrnb,int maxnit,real tol,
		   rvec box)
{
  /* A simple Gauss-Seidel relaxation method for solving the Poisson
   * equation: Nabla^2 Pot = -Rho
   * (epsilon0 should be in Rho in the case of coulomb potential)
   * Successive overrelaxation is applied, which speeds up considerably.
   *
   */
  static bool bFirst = TRUE;
  static int  *nnx,*nny,*nnz;
  static real fac,dx_2,dy_2,dz_2,fac_1,omega;
  real   deviation,val_ijk,epsrho;
  real   dx2,dy2,dz2,residual,sum,xi;
  int    nit;
  int    i,j,k,i_1,i1,j_1,j1,k_1,k1;
  int    nx,ny,nz;
  real   ***pot_ptr,***rho_ptr;
  
  unpack_PSgrid(pot,&nx,&ny,&nz,&pot_ptr);
  unpack_PSgrid(rho,&nx,&ny,&nz,&rho_ptr);
  
  if (bFirst) {
    fprintf(log,"Solving Poisson Equation on %dx%dx%d grid\n",nx,ny,nz);
    calc_nxyz(nx,ny,nz,&nnx,&nny,&nnz);
    
    /* Prefactors corresponding to grid spacing squared for 
     * finite differencing. See Num. Res. Chap. 19.
     * Note that grid spacing may not be identical in all directions.
     */
    dx2   = sqr(box[XX]/nx);
    dy2   = sqr(box[YY]/ny);
    dz2   = sqr(box[ZZ]/nz);
    dx_2  = 1.0/dx2;
    dy_2  = 1.0/dy2;
    dz_2  = 1.0/dz2;
    fac_1 = 2.0*(dx_2+dy_2+dz_2);
    fac   = 1.0/fac_1;

    /* SOR factor from Atkins: 
     * An Introduction to Numerical Analysis (p. 486) 
     */
    xi    = 1-2*sqr(sin(M_PI/(2*nx)));
    omega = 2/(1+sqrt(1-xi*xi));
    fprintf(log,"xi = %g, omega = %g\n",xi,omega);

    bFirst = FALSE;
  }

  /* Solve by simple averaging */
  nit  = 0;
  do {
    deviation = 0.0;
    sum = 0.0;
    for(i=0; (i<nx); i++) {
      /* Indices left and right in the grid */
      i_1 = nnx[i-1+nx];
      i1  = nnx[i+1+nx];
      for(j=0; (j<ny); j++) {
	/* Indices up and down in the grid */
	j_1 = nny[j-1+ny];
	j1  = nny[j+1+ny];
	for(k=0; (k<nz); k++) {
	  /* Indices fore and back in the grid */
	  k_1 = nnz[k-1+nz];
	  k1  = nnz[k+1+nz];
	  
	  /* Get the new value by averaging surrounding grid points */
	  epsrho  = rho_ptr[i][j][k];
	  val_ijk = pot_ptr[i][j][k];
	      
	  /* Calculate the error in the current potential */
	  residual= (dx_2*(pot_ptr[i_1][j][k] + pot_ptr[i1][j][k]) +
		     dy_2*(pot_ptr[i][j_1][k] + pot_ptr[i][j1][k]) +
		     dz_2*(pot_ptr[i][j][k_1] + pot_ptr[i][j][k1]) -
		     fac_1*val_ijk + epsrho);
	  
	  deviation        += fabs(residual);
	  pot_ptr[i][j][k] += omega*fac*residual;
	  sum              += pot_ptr[i][j][k];
	}
      }
    }
    
    /* This call may or may not be necessary. For a 32 cube grid and coulomb
     * interaction it does NOT make a difference in the forces nor in the
     * potential. It is quite expensive however...
     */
    symmetrize_PSgrid(NULL,pot,sum);
    
    deviation = (deviation/(nx*ny*nz));

    if (bVerbose)
      fprintf(stderr,"\rnit: %5d  dev: %8.3f  sum: %8.3f\n",
	      nit,deviation,sum);
    
    nit ++;
  } while ((nit < maxnit) && (deviation > tol));
  if (bVerbose)
    fprintf(stderr,"\n");
  if (nit == maxnit)
    fatal_error(0,"Poisson Solver did *not* converge in %d iterations\n",nit);
  
  return nit;
}









