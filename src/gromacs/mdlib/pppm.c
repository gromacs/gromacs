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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include "physics.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "xvgr.h"
#include "gmx_fatal.h"
#include "txtdump.h"
#include "network.h"
#include "nrnb.h"
#include "pppm.h"
#include "coulomb.h"
#include "mdrun.h"
#include "gmx_fft.h"
#include "pme.h"

#ifdef GMX_MPI
#include "gmx_parallel_3dfft.h"
#endif

#define llim  (-1)
#define ulim   (1)
#define llim2 (-3)
#define ulim2  (3)



/* PPPM temporarily disabled while working on 2D PME */
#define DISABLE_PPPM



#ifndef DISABLE_PPPM

/* TODO: fix thread-safety */

static void calc_invh(rvec box,int nx,int ny,int nz,rvec invh)
{
  invh[XX] = nx/box[XX];
  invh[YY] = ny/box[YY];
  invh[ZZ] = nz/box[ZZ];
}

static void calc_weights(int iatom,int nx,int ny,int nz,
			 rvec x,rvec box,rvec invh,ivec ixyz,real WXYZ[])
{
  const  real half=0.5;
  tensor wxyz;
  real   abc,ttt,fact;
#ifdef DEBUG
  real   wtot;
#endif
  ivec   nxyz;
  int    it,j,k,m,nm;
  real   Wx,Wy,Wzx,Wzy,Wzz;
  
  fact = 0.125;
  
  nxyz[XX] = nx;
  nxyz[YY] = ny;
  nxyz[ZZ] = nz;
  for(m=0; (m<DIM); m++) {
    /* Put particle in the box... */  
    ttt = x[m]*invh[m];
    
    /* Round to nearest grid point. Do the math in integer! */
    it  = (ttt+0.5);
    nm  = nxyz[m];
    if (it < 0) {
      ttt+=nm;
      it +=nm;
    }
    else if (it >= nm) {
      ttt-=nm;
      it -=nm;
    }
    if ((it < 0) || (it >= nxyz[m]))
      gmx_fatal(FARGS,"iatom = %d, it = %d, x=%f, ttt=%f",iatom,it,x[m],ttt);
    ixyz[m]    = it;
    
    /* Fraction (offset) from grid point */
    abc        = ttt - (real)ixyz[m];
    
    wxyz[m][0] = sqr((real)(half  - abc));
    wxyz[m][1] = 1.5 - 2.0*sqr(abc);
    wxyz[m][2] = sqr((real)(half  + abc));
  }
  Wzx=wxyz[ZZ][XX];
  Wzy=wxyz[ZZ][YY];
  Wzz=wxyz[ZZ][ZZ];
  for(j=m=0; (j<DIM); j++) {
    Wx = wxyz[XX][j]*fact;
    for(k=0; (k<DIM); k++,m+=3) {
      Wy        = Wx*wxyz[YY][k];
      WXYZ[m]   = Wy*Wzx;
      WXYZ[m+1] = Wy*Wzy;
      WXYZ[m+2] = Wy*Wzz;
    }
  }
#ifdef DEBUG
  wtot = 0;
  for(j=0; (j<27); j++)
    wtot+=WXYZ[j];
  fprintf(stderr,"wtot = %g\n",wtot);
#endif
#ifdef HACK
  for(j=0; (j<27); j++)
    WXYZ[j] = 0;
  WXYZ[13] = 1.0;
#endif
}

static void calc_nxyz(int nx,int ny,int nz,
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
	
static void spread_q(FILE *log,gmx_bool bVerbose,
		     int start,int nr,
		     rvec x[],real charge[],rvec box,
		     t_fftgrid *grid,t_nrnb *nrnb)
{
  static gmx_bool bFirst = TRUE;
  static int  *nnx,*nny,*nnz;
  rvec   invh;
  real   qi,qwt;
#ifdef DEBUG
  real   qt;
#endif
  real   WXYZ[27];
  ivec   ixyz;
  int    i,iX,iY,iZ,index;
  int    jx,jy,jz,jcx,jcy,jcz;
  int    nxyz;
  int    nx,ny,nz,nx2,ny2,nz2,la2,la12;
  real *ptr;
  
  unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
  
  calc_invh(box,nx,ny,nz,invh);

  if (bFirst) {
    if (log) {
      fprintf(log,
	      "Spreading Charges using Triangle Shaped on %dx%dx%d grid\n",
	      nx,ny,nz);
      fprintf(log,"invh = %10g,%10g,%10g\n",invh[XX],invh[YY],invh[ZZ]);
    }
  
    calc_nxyz(nx,ny,nz,&nnx,&nny,&nnz);
    
    bFirst = FALSE;
  }

  for(i=start; (i<start+nr); i++) {
    qi=charge[i];

    /* Each charge is spread over the nearest 27 grid cells,
     * thus we loop over -1..1 in X,Y and Z direction
     * We apply the TSC (triangle shaped charge)
     * see Luty et. al, JCP 103 (1995) 3014
     */
    
    if (qi != 0.0) {
      calc_weights(i,nx,ny,nz,x[i],box,invh,ixyz,WXYZ);
      iX  = ixyz[XX] + nx;
      iY  = ixyz[YY] + ny;
      iZ  = ixyz[ZZ] + nz;

#ifdef DEBUG
      qt=0;
#endif
      nxyz = 0;
      for(jx=-1; (jx<=1); jx++) {
	jcx = nnx[iX + jx];
	for(jy=-1; (jy<=1); jy++) {
	  jcy = nny[iY + jy];
	  for(jz=-1; (jz<=1); jz++,nxyz++) {
	    jcz   = nnz[iZ + jz];
	    index = INDEX(jcx,jcy,jcz);
	    qwt   = qi*WXYZ[nxyz];
	    grid->ptr[index]+=qwt;
#ifdef DEBUG
	    qt   += qwt;
	    if (bVerbose && log)
	      fprintf(log,"spread %4d %2d %2d %2d  %10.3e, weight=%10.3e\n",
		      index,jcx,jcy,jcz,grid->ptr[index],WXYZ[nxyz]);
#endif
	  }
	}
      }
#ifdef DEBUG
      if (log)
	fprintf(log,"q[%3d] = %6.3f, qwt = %6.3f\n",i,qi,qt);
#endif
    }
  }
  inc_nrnb(nrnb,eNR_SPREADQ,9*nr);
  inc_nrnb(nrnb,eNR_WEIGHTS,3*nr);
}

static real gather_inner(int JCXYZ[],real WXYZ[],int ixw[],int iyw[],int izw[],
			 int la2,int la12,
		  real c1x,real c1y,real c1z,real c2x,real c2y,real c2z,
		  real qi,rvec f,real ptr[])
{
  real pi,fX,fY,fZ,weight;
  int  jxyz,m,jcx,jcy,jcz;
  int  jcx0,jcy0,jcz0;
  
  pi = 0.0;
  fX = 0.0;
  fY = 0.0;
  fZ = 0.0;
  
  /* Now loop over 27 surrounding vectors */      
  for(jxyz=m=0; (jxyz < 27); jxyz++,m+=3) {
    jcx    = JCXYZ[m];
    jcy    = JCXYZ[m+1];
    jcz    = JCXYZ[m+2];
    weight = WXYZ[jxyz];
    
    jcx0   = ixw[jcx];
    jcy0   = iyw[jcy];
    jcz0   = izw[jcz];

    /* Electrostatic Potential ! */
    pi += weight * ptr[INDEX(jcx0,jcy0,jcz0)];

    /* Forces */
    fX += weight * ((c1x*(ptr[INDEX(ixw[jcx-1],jcy0,jcz0)] - 
			  ptr[INDEX(ixw[jcx+1],jcy0,jcz0)] )) +
		    (c2x*(ptr[INDEX(ixw[jcx-2],jcy0,jcz0)] - 
			  ptr[INDEX(ixw[jcx+2],jcy0,jcz0)] )));
    fY += weight * ((c1y*(ptr[INDEX(jcx0,iyw[jcy-1],jcz0)] -
			  ptr[INDEX(jcx0,iyw[jcy+1],jcz0)] ))  +
		    (c2y*(ptr[INDEX(jcx0,iyw[jcy-2],jcz0)] -
			  ptr[INDEX(jcx0,iyw[jcy+2],jcz0)] )));
    fZ += weight * ((c1z*(ptr[INDEX(jcx0,jcy0,izw[jcz-1])] -
			  ptr[INDEX(jcx0,jcy0,izw[jcz+1])] ))  +
		    (c2z*(ptr[INDEX(jcx0,jcy0,izw[jcz-2])] -
			  ptr[INDEX(jcx0,jcy0,izw[jcz+2])] )));
  }
  f[XX] += qi*fX;
  f[YY] += qi*fY;
  f[ZZ] += qi*fZ;
  
  return pi;
}

static real gather_f(FILE *log,gmx_bool bVerbose,
		     int start,int nr,rvec x[],rvec f[],real charge[],rvec box,
		     real pot[],t_fftgrid *grid,rvec beta,t_nrnb *nrnb)
{
  static gmx_bool bFirst=TRUE;
  static int  *nnx,*nny,*nnz;
  static int  JCXYZ[81];
  int    i,m;
  real   energy;
  real   qi,pi;
  ivec   ixyz;
  rvec   invh,c1,c2;
  real   WXYZ[27];
  real   c1x,c1y,c1z,c2x,c2y,c2z;
  int    ixw[7],iyw[7],izw[7];
  int    ll;
  int    nx,ny,nz,nx2,ny2,nz2,la2,la12;
  real *ptr;
  
  unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,&la2,&la12,TRUE,&ptr);
  
  calc_invh(box,nx,ny,nz,invh);
  
  for(m=0; (m<DIM); m++) {
    c1[m] = (beta[m]/2.0)*invh[m];
    c2[m] = ((1.0-beta[m])/4.0)*invh[m];
  }
  c1x = c1[XX];
  c1y = c1[YY];
  c1z = c1[ZZ];
  c2x = c2[XX];
  c2y = c2[YY];
  c2z = c2[ZZ];

  if (bFirst) {
    if (log) {
      fprintf(log,"Gathering Forces using Triangle Shaped on %dx%dx%d grid\n",
	      nx,ny,nz);
      fprintf(log,"beta = %10g,%10g,%10g\n",beta[XX],beta[YY],beta[ZZ]);
      fprintf(log,"c1   = %10g,%10g,%10g\n",c1[XX],c1[YY],c1[ZZ]);
      fprintf(log,"c2   = %10g,%10g,%10g\n",c2[XX],c2[YY],c2[ZZ]);
      fprintf(log,"invh = %10g,%10g,%10g\n",invh[XX],invh[YY],invh[ZZ]);
    }
    calc_nxyz(nx,ny,nz,&nnx,&nny,&nnz);

    for(i=0; (i<27); i++) {
      JCXYZ[3*i]   = 2 + (i/9);
      JCXYZ[3*i+1] = 2 + (i/3) % 3;
      JCXYZ[3*i+2] = 2 + (i % 3); 
    }
    
    bFirst = FALSE;
  }

  energy=0.0;  	  
  for(i=start; (i<start+nr); i++) {
    /* Each charge is spread over the nearest 27 grid cells,
     * thus we loop over -1..1 in X,Y and Z direction
     * We apply the TSC (triangle shaped charge)
     * see Luty et. al, JCP 103 (1995) 3014
     */
     
    calc_weights(i,nx,ny,nz,x[i],box,invh,ixyz,WXYZ);

    for(ll=llim2; (ll<=ulim2); ll++) {
      ixw[ll-llim2] = nnx[ixyz[XX]+ll+nx];
      iyw[ll-llim2] = nny[ixyz[YY]+ll+ny];
      izw[ll-llim2] = nnz[ixyz[ZZ]+ll+nz];
    }
    
    qi      = charge[i];
    pi      = gather_inner(JCXYZ,WXYZ,ixw,iyw,izw,la2,la12,
			   c1x,c1y,c1z,c2x,c2y,c2z,
			   qi,f[i],ptr);
    
    energy += pi*qi;
    pot[i]  = pi;
  }
  
  inc_nrnb(nrnb,eNR_GATHERF,27*nr);
  inc_nrnb(nrnb,eNR_WEIGHTS,3*nr);
  
  return energy*0.5;
}

static void convolution(FILE *fp,gmx_bool bVerbose,t_fftgrid *grid,real ***ghat,
			t_commrec *cr)
{
  int      i,j,k,index;
  real     gk;
  int      nx,ny,nz,nx2,ny2,nz2,la2,la12;
  t_complex  *ptr;
  int      *nTest;
  int jstart,jend;
  
  unpack_fftgrid(grid,&nx,&ny,&nz,&nx2,&ny2,&nz2,
		 &la2,&la12,FALSE,(real **)&ptr);
  snew(nTest,grid->nptr);
  
  if(PAR(cr)) {
#if (defined GMX_MPI && !defined GMX_WITHOUT_FFTW)
    jstart=grid->pfft.local_y_start_after_transpose;
    jend=jstart+grid->pfft.local_ny_after_transpose;

    for(j=jstart; (j<jend); j++) { /* local cells */
	for(i=0; (i<nx); i++) {
	    for(k=0;k<(nz/2+1); k++) {
		gk    = ghat[i][j][k];
		index = INDEX(j,i,k);
		ptr[index].re *= gk;
		ptr[index].im *= gk;
		nTest[index]++;
	    }
	}
    }
#ifdef DEBUG
    for(j=jstart; (j<jend); j++) {
	for(i=0; (i<nx); i++) {
	    for(k=0; k<(nz/2+1); k++) {
		index = INDEX(j,i,k);
		if (nTest[index] != 1)
		    fprintf(fp,"Index %d sucks, set %d times\n",index,nTest[index]);
	    }
	}
    }
#endif /* DEBUG */
#endif /* GMX_MPI */
  } else { /* if not running in parallel */
      for(i=0; (i<nx); i++) {
	  for(j=0; (j<ny); j++) {
	      for(k=0;k<(nz/2+1); k++) {
		  gk    = ghat[i][j][k];
		  index = INDEX(i,j,k);
		  ptr[index].re *= gk;
		  ptr[index].im *= gk;
		  nTest[index]++;
	      }
	  }
      }
#ifdef DEBUG
      for(i=0; (i<nx); i++) {
	  for(j=0; (j<ny); j++) {
	      for(k=0; k<(nz/2+1); k++) {
		  index = INDEX(i,j,k);
		  if (nTest[index] != 1)
		      fprintf(fp,"Index %d sucks, set %d times\n",index,nTest[index]);
	      }
	  }
      }
#endif	
  }
  sfree(nTest);
}


void solve_pppm(FILE *fp,t_commrec *cr,
		t_fftgrid *grid,real ***ghat,rvec box,
		gmx_bool bVerbose,t_nrnb *nrnb)
{
  int  ntot,npppm;
  
  /*  if (bVerbose) 
      print_fftgrid(fp,"Q-Real",grid,grid->nxyz,"qreal.pdb",box,TRUE);*/
  
  gmxfft3D(grid,GMX_FFT_REAL_TO_COMPLEX,cr);
  
  /*  if (bVerbose) {
      print_fftgrid(fp,"Q-k",grid,1.0,"qk-re.pdb",box,TRUE);
      print_fftgrid(fp,"Q-k",grid,1.0,"qk-im.pdb",box,FALSE);
      fprintf(stderr,"Doing convolution\n");
      }*/
  
  convolution(fp,bVerbose,grid,ghat,cr); 
  
  /*  if (bVerbose) 
      print_fftgrid(fp,"Convolution",grid,1.0,
      "convolute.pdb",box,TRUE);*/
  
  gmxfft3D(grid,GMX_FFT_COMPLEX_TO_REAL,cr);
  
  /*  if (bVerbose) 
      print_fftgrid(fp,"Potential",grid,1.0,"pot.pdb",box,TRUE);*/
  
  ntot  = grid->nxyz;  
  npppm = ntot*log((real)ntot)/log(2.0);
  inc_nrnb(nrnb,eNR_FFT,2*npppm);
  inc_nrnb(nrnb,eNR_CONV,ntot);
}


static rvec      beta;
static real      ***ghat=NULL;
static t_fftgrid *grid=NULL;

#endif


int gmx_pppm_init(FILE *log,      t_commrec *cr,
                  const output_env_t oenv, gmx_bool bVerbose,
                  gmx_bool bOld,      matrix box,
                  char *ghatfn,   t_inputrec *ir,
                  gmx_bool bReproducible)
{
  int   nx,ny,nz,m,porder;
  ivec  grids;
  real  r1,rc;
  const real tol = 1e-5;
  rvec  box_diag,spacing;

#ifdef DISABLE_PPPM
    gmx_fatal(FARGS,"PPPM is not functional in the current version, we plan to implement PPPM through a small modification of the PME code.");
    return -1;
#else
    
#ifdef GMX_WITHOUT_FFTW
  gmx_fatal(FARGS,"PPPM used, but GROMACS was compiled without FFTW support!\n");
#endif

  if (log) {
    if (cr != NULL) {
      if (cr->nnodes > 1)
	fprintf(log,"Initializing parallel PPPM.\n");
    }
    fprintf(log,
	    "Will use the PPPM algorithm for long-range electrostatics\n");
  }
 
  for(m=0; m<DIM; m++)
    box_diag[m] = box[m][m];

  if (!gmx_fexist(ghatfn)) {    
    beta[XX]=beta[YY]=beta[ZZ]= 1.85;
    nx     = ir->nkx;
    ny     = ir->nky;
    nz     = ir->nkz;
   
    if (log) {
      fprintf(log,"Generating Ghat function\n");
      fprintf(log,"Grid size is %d x %d x %d\n",nx,ny,nz);
    }

    if ((nx < 4) || (ny < 4) || (nz < 4)) 
      gmx_fatal(FARGS,"Grid must be at least 4 points in all directions");
      
    ghat   = mk_rgrid(nx,ny,nz);
    mk_ghat(NULL,nx,ny,nz,ghat,box_diag,
	    ir->rcoulomb_switch,ir->rcoulomb,TRUE,bOld);
    
    if (bVerbose)
      pr_scalar_gk("generghat.xvg",oenv,nx,ny,nz,box_diag,ghat);
  }
  else {
    fprintf(stderr,"Reading Ghat function from %s\n",ghatfn);
    ghat = rd_ghat(log,oenv,ghatfn,grids,spacing,beta,&porder,&r1,&rc);
    
    /* Check whether cut-offs correspond */
    if ((fabs(r1-ir->rcoulomb_switch)>tol) || (fabs(rc-ir->rcoulomb)>tol)) {
      if (log) {
	fprintf(log,"rcoulomb_switch = %10.3e  rcoulomb = %10.3e"
		"  r1 = %10.3e  rc = %10.3e\n",
		ir->rcoulomb_switch,ir->rcoulomb,r1,rc);
	fflush(log);
      }
      gmx_fatal(FARGS,"Cut-off lengths in tpb file and Ghat file %s "
		  "do not match\nCheck your log file!",ghatfn);
    }
      
    /* Check whether boxes correspond */
    for(m=0; (m<DIM); m++)
      if (fabs(box_diag[m]-grids[m]*spacing[m]) > tol) {
	if (log) {
	  pr_rvec(log,0,"box",box_diag,DIM,TRUE);
	  pr_rvec(log,0,"grid-spacing",spacing,DIM,TRUE);
	  pr_ivec(log,0,"grid size",grids,DIM,TRUE);
	  fflush(log);
	}
	gmx_fatal(FARGS,"Box sizes in tpb file and Ghat file %s do not match\n"
		    "Check your log file!",ghatfn);
      }

    if (porder != 2)
      gmx_fatal(FARGS,"porder = %d, should be 2 in %s",porder,ghatfn);
      
    nx = grids[XX];
    ny = grids[YY];
    nz = grids[ZZ];
    
    if (bVerbose)
      pr_scalar_gk("optimghat.xvg",oenv,nx,ny,nz,box_diag,ghat);
  }
  /* Now setup the FFT things */
#ifdef GMX_MPI
  cr->mpi_comm_mygroup=cr->mpi_comm_mysim;
#endif
  grid = mk_fftgrid(nx,ny,nz,NULL,NULL,cr,bReproducible);
  
  return 0;
#endif
}

int gmx_pppm_do(FILE *log,       gmx_pme_t pme,
		gmx_bool bVerbose,
		rvec x[],        rvec f[],
		real charge[],   rvec box,
		real phi[],      t_commrec *cr,
		int start,       int nr,
		t_nrnb *nrnb,
		int pme_order,   real *energy)
{
#ifdef DISABLE_PPPM
    gmx_fatal(FARGS,"PPPM temporarily disabled while working on 2DPME\n");
    return -1;
#else

    /* Make the grid empty */
  clear_fftgrid(grid);
  
  /* First step: spreading the charges over the grid. */
  spread_q(log,bVerbose,start,nr,x,charge,box,grid,nrnb);
  
  /* In the parallel code we have to sum the grids from neighbouring nodes */
  if (PAR(cr))
    gmx_sum_qgrid(pme,cr,grid,GMX_SUM_QGRID_FORWARD);
  
  /* Second step: solving the poisson equation in Fourier space */
  solve_pppm(log,cr,grid,ghat,box,bVerbose,nrnb);
  
  /* In the parallel code we have to sum once again... */
  if (PAR(cr))
    gmx_sum_qgrid(pme,cr,grid,GMX_SUM_QGRID_BACKWARD);
  
  /* Third and last step: gather the forces, energies and potential
   * from the grid.
   */
  *energy = gather_f(log,bVerbose,start,nr,x,f,charge,box,
		     phi,grid,beta,nrnb);
  
  return 0;
#endif
}

#ifndef DISABLE_PPPM
static int gmx_pppm_opt_do(FILE *log,       gmx_pme_t pme,
			   t_inputrec *ir,  gmx_bool bVerbose,
			   int natoms,
			   rvec x[],        rvec f[],
			   real charge[],   rvec box,
			   real phi[],      t_commrec *cr,
			   t_nrnb *nrnb,    rvec beta,
			   t_fftgrid *grid, gmx_bool bOld,
			   real *energy)
{
  real      ***ghat;
  int       nx,ny,nz;
  
  if (log)
    fprintf(log,"Generating Ghat function\n");
  nx     = ir->nkx;
  ny     = ir->nky;
  nz     = ir->nkz;
  ghat   = mk_rgrid(nx,ny,nz);
  mk_ghat(NULL,nx,ny,nz,ghat,box,ir->rcoulomb_switch,ir->rcoulomb,TRUE,bOld);
  
  /* pr_scalar_gk("generghat.xvg",nx,ny,nz,box,ghat); */
  
  /* Now start the actual PPPM procedure.
   * First step: spreading the charges over the grid.
   */
  /* Make the grid empty */
  clear_fftgrid(grid);
  
  spread_q(log,bVerbose,0,natoms,x,charge,box,grid,nrnb);
  
  /* Second step: solving the poisson equation in Fourier space */
  solve_pppm(log,cr,grid,ghat,box,bVerbose,nrnb);
  
  /* Third and last step: gather the forces, energies and potential
   * from the grid.
   */
  *energy = gather_f(log,bVerbose,0,natoms,x,f,charge,box,phi,grid,beta,nrnb);

  free_rgrid(ghat,nx,ny);
    
  return 0;
}

#endif
