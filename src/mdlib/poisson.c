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
static char *SRCID_poisson_c = "$Id$";
#include <math.h>
#include "smalloc.h"
#include "typedefs.h"
#include "shift_util.h"
#include "macros.h"
#include "physics.h"
#include "nrnb.h"
#include "vec.h"
#include "pppm.h"
#include "poisson.h"

t_PSgrid *mk_PSgrid(int nx,int ny,int nz)
{
  t_PSgrid *ps;
  int      i,j;
  
  snew(ps,1);
  ps->nx=nx;
  ps->ny=ny;
  ps->nz=nz;
  snew(ps->ptr,nx);
  for(i=0; (i<nx); i++) {
    snew(ps->ptr[i],ny);
    for(j=0; (j<ny); j++)
      snew(ps->ptr[i][j],nz);
  }
  return ps;
}

void unpack_PSgrid(t_PSgrid *grid,int *nx,int *ny,int *nz,real ****ptr)
{
  *nx  = grid->nx;
  *ny  = grid->ny;
  *nz  = grid->nz;
  *ptr = grid->ptr;
}

void copy_PSgrid(t_PSgrid *dest,t_PSgrid *src)
{
  int  i,j,k;
  int  nx,ny,nz;
  real ***src_ptr,***dst_ptr;
  
  unpack_PSgrid(dest,&nx,&ny,&nz,&dst_ptr);
  unpack_PSgrid(src,&nx,&ny,&nz,&src_ptr);
  
  for(i=0; (i<nx); i++)
    for(j=0; (j<ny); j++)
      for(k=0; (k<nz); k++)
	dst_ptr[i][j][k] = src_ptr[i][j][k];
}

void clear_PSgrid(t_PSgrid *grid)
{
  int  i,j,k;
  int  nx,ny,nz;
  real ***ptr;
  
  unpack_PSgrid(grid,&nx,&ny,&nz,&ptr);
  
  for(i=0; (i<nx); i++)
    for(j=0; (j<ny); j++)
      for(k=0; (k<nz); k++)
	ptr[i][j][k] = 0.0;
}

void symmetrize_PSgrid(FILE *fp,t_PSgrid *grid,real sum)
{
  int  i,j,k;
  int  nx,ny,nz;
  real ***ptr;
  real ming=0,maxg=0;
  ivec imin={-1,-1,-1},imax={-1,-1,-1};
  
  unpack_PSgrid(grid,&nx,&ny,&nz,&ptr);

  ming = maxg = 0;  
  if (sum == 0.0) {
    sum  = 0;
    ming = maxg = ptr[0][0][0];
    for(i=0; (i<nx); i++)
      for(j=0; (j<ny); j++)
	for(k=0; (k<nz); k++) {
	  sum += ptr[i][j][k];
	  if (ptr[i][j][k] < ming) {
	    ming     = ptr[i][j][k];
	    imin[XX] = i;
	    imin[YY] = j;
	    imin[ZZ] = k;
	  }
	  else if (ptr[i][j][k] > maxg) {
	    maxg     = ptr[i][j][k];
	    imax[XX] = i;
	    imax[YY] = j;
	    imax[ZZ] = k;
	  }
	}
  }
  sum /= (nx*ny*nz);
  for(i=0; (i<nx); i++)
    for(j=0; (j<ny); j++)
      for(k=0; (k<nz); k++)
	ptr[i][j][k] -= sum;
  if (fp)
    fprintf(fp,"Symmetrize_PSgrid: sum = %g, ming = %g(%d,%d,%d), maxg = %g(%d,%d,%d)\n",
	    sum,
	    ming,imin[XX],imin[YY],imin[ZZ],
	    maxg,imax[XX],imax[YY],imax[ZZ]);
}

void calc_nxyz(int nx,int ny,int nz,
	       int **nnx,int **nny,int **nnz)
{
  int i;
  
  snew(*nnx,3*nx);
  snew(*nny,3*ny);
  snew(*nnz,3*nz);
  for(i=0; (i<3*nx); i++)
    (*nnx)[i] = i % nx;
  for(i=0; (i<3*ny); i++)
    (*nny)[i] = i % ny;
  for(i=0; (i<3*nz); i++)
    (*nnz)[i] = i % nz;
}
	
real do_poisson(FILE *log,       bool bVerbose,
		t_inputrec *ir,  int natoms,
		rvec x[],        rvec f[],
		real charge[],   rvec box,
		real phi[],      t_commrec *cr,
		t_nrnb *nrnb,    int *nit,   
		bool bOld)
{
  static  bool bFirst  = TRUE;
  static  t_PSgrid *pot,*rho;
  static  int       maxnit;
  static  real      r1,rc;
  static  rvec      beta;
  
  const     real tol = 1e-2;
  int       m;
  real      ener;
  
  ener = 0.0;
  
  if (bFirst) {
    maxnit = ir->userint1;

    fprintf(log,"Will use Poisson Solver for long-range electrostatics\n");
    fprintf(log,"Grid size is %d x %d x %d\n",ir->nkx,ir->nky,ir->nkz);

    if ((ir->nkx < 4) || (ir->nky < 4) || (ir->nkz < 4)) 
      fatal_error(0,"Grid must be at least 4 points in all directions");
    
    pot = mk_PSgrid(ir->nkx,ir->nky,ir->nkz);
    rho = mk_PSgrid(ir->nkx,ir->nky,ir->nkz);
    
    r1 = ir->rcoulomb_switch;
    rc = ir->rcoulomb;
    for(m=0; (m<DIM); m++)
      beta[m] = 1.85;
      
  }

  /* Make the grid empty and spread the charges */
  clear_PSgrid(rho);
  spread_q_poisson(log,bVerbose,TRUE,natoms,x,charge,box,rc,rho,nrnb,bOld,r1);
  
  /* Subtract average charge distribution from each grid point. This does not
   * influence the forces, but it may influence the potential and the energy.
   * On the other hand the charge distribution should sum up to zero if the
   * system is neutral!
   */
  symmetrize_PSgrid(debug,rho,0.0);
  
  /* For the first time we solve the potential we copy the charge distribution to
   * the potential as a start for the solver. Later we use the previous step's
   * solution.
   */
  if (bFirst)
    copy_PSgrid(pot,rho);
  
  /* Second step: solving the poisson equation in real space */
  *nit = solve_poisson(log,pot,rho,bVerbose,nrnb,maxnit,tol,box);
  
  symmetrize_PSgrid(debug,pot,0.0);
  
  /* Third and last step: gather the forces, energies and potential
   * from the grid.
   */
  ener = ps_gather_f(log,bVerbose,natoms,x,f,charge,box,phi,pot,beta,nrnb);
  
  bFirst = FALSE;
  
  return ener;
}

real do_optimize_poisson(FILE *log,       bool bVerbose,
			 t_inputrec *ir,  int natoms,
			 rvec x[],        rvec f[],
			 real charge[],   rvec box,
			 real phi[],      t_commrec *cr,
			 t_nrnb *nrnb,    rvec f_ref[],
			 real phi_ref[],  rvec beta,
			 bool bOld)
{
#define BMIN 1.6
#define DB   0.025
#define NB   (1+(int)((2.1-BMIN)/DB))
  static  bool bFirst  = TRUE;
  static  bool bSecond = TRUE;
  static  t_PSgrid *pot,*rho;
  static  int       maxnit;
  static  real      r1,rc;
  
  real      rmsf[NB][NB][NB],rmsf_min,rrmsf;
  ivec      minimum;
  const     real tol = 1e-2;
  int       i,m,bx,by,bz;
  char      buf[128];
  real      ener;
  
  ener = 0.0;
  
  if (bFirst) {
    maxnit = ir->userint1;

    fprintf(log,"Will use Poisson Solver for long-range electrostatics\n");
    fprintf(log,"Grid size is %d x %d x %d\n",ir->nkx,ir->nky,ir->nkz);

    if ((ir->nkx < 4) || (ir->nky < 4) || (ir->nkz < 4)) 
      fatal_error(0,"Grid must be at least 4 points in all directions");
    
    pot = mk_PSgrid(ir->nkx,ir->nky,ir->nkz);
    rho = mk_PSgrid(ir->nkx,ir->nky,ir->nkz);
    
    r1 = ir->rcoulomb_switch;
    rc = ir->rcoulomb;
    for(m=0; (m<DIM); m++)
      beta[m] = 4.0/3.0;
      
    bFirst = FALSE;
  }

  /* Make the grid empty */
  clear_PSgrid(rho);
  spread_q_poisson(log,bVerbose,TRUE,natoms,x,charge,box,rc,rho,nrnb,
		   bOld,r1);
  
  symmetrize_PSgrid(log,rho,0.0);
  if (bSecond) 
    copy_PSgrid(pot,rho);
  
  /* Second step: solving the poisson equation in real space */
  (void) solve_poisson(log,pot,rho,bVerbose,nrnb,maxnit,tol,box);
  
  symmetrize_PSgrid(log,pot,0.0);
  /* Third and last step: gather the forces, energies and potential
   * from the grid.
   */
#define BETA(n) (BMIN+n*DB)
  /* Optimization of beta in progress */
  for(bx=0; (bx<NB); bx++) {
    beta[XX] = BETA(bx);
    for(by=0; (by<NB); by++) {
      beta[YY] = BETA(by);
      for(bz=0; (bz<NB); bz++) {
	beta[ZZ] = BETA(bz);
	  
	for(i=0; (i<natoms); i++) {
	  phi[i] = 0.0;
	  clear_rvec(f[i]);
	}
	ener = ps_gather_f(log,bVerbose,natoms,x,f,charge,box,
			   phi,pot,beta,nrnb);
	sprintf(buf,"Poisson, beta = %g\n",beta[XX]);
	rmsf[bx][by][bz] = 
	  analyse_diff(log,buf,natoms,f_ref,f,phi_ref,phi,NULL,
		       /*"fcorr.xvg","pcorr.xvg"*/NULL,NULL,NULL,NULL);
      }
    }
  }
  rmsf_min = rmsf[0][0][0];
  minimum[XX] = minimum[YY] = minimum[ZZ] = 0;
  for(bx=0; (bx<NB); bx++) {
    beta[XX] = BETA(bx);
    for(by=0; (by<NB); by++) {
      beta[YY] = BETA(by);
      for(bz=0; (bz<NB); bz++) {
	beta[ZZ] = BETA(bz);
	rrmsf    = rmsf[bx][by][bz];
	
	fprintf(log,"Beta: %6.3f  %6.3f  %6.3f  RMSF: %8.3f\n",
		beta[XX],beta[YY],beta[ZZ],rrmsf);
	if (rrmsf < rmsf_min) {
	  rmsf_min = rrmsf;
	  minimum[XX] = bx;
	  minimum[YY] = by;
	  minimum[ZZ] = bz;
	}
      }
    }
  }
  beta[XX] = BETA(minimum[XX]);
  beta[YY] = BETA(minimum[YY]);
  beta[ZZ] = BETA(minimum[ZZ]);
  fprintf(log,"Minimum RMSF %8.3f at Beta = %6.3f  %6.3f  %6.3f\n",
	  rmsf_min,beta[XX],beta[YY],beta[ZZ]);
  /* Computing optimum once more... */
  for(i=0; (i<natoms); i++) {
    phi[i] = 0.0;
    clear_rvec(f[i]);
  }
  ener = ps_gather_f(log,bVerbose,natoms,x,f,charge,box,phi,pot,beta,nrnb);

  return ener;
}
