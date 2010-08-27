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
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "futil.h"
#include "vec.h"
#include "physics.h"
#include "coulomb.h"
#include "pppm.h"
#include "xvgr.h"
#include "gmxfio.h"
#include "pppm.h"
#include "smalloc.h"

static void calc_k(rvec lll,int ix,int iy,int iz,int nx,int ny,int nz,rvec k)
{
#define IDX(i,n,x)  (i<=n/2) ? (i*x) : ((i-n)*x)
  k[XX] = IDX(ix,nx,lll[XX]);
  k[YY] = IDX(iy,ny,lll[YY]);
  k[ZZ] = IDX(iz,nz,lll[ZZ]);
#undef IDX
}

void symmetrize_ghat(int nx,int ny,int nz,real ***ghat)
/* Symmetrize the Ghat function. It is assumed that the 
 * first octant of the Ghat function is either read or generated
 * (all k-vectors from 0..nx/2 0..ny/2 0..nz/2).
 * Since Gk depends on the absolute value of k only, 
 * symmetry operations may shorten the time to generate it.
 */
 
{
  int  i,j,k;
  int  iip,jjp,kkp;
  real ggg;

  /* fprintf(stderr,"Symmetrizing Ghat function\n");   */
  /* Only the lower octant of the rectangle has been saved,
   * so we must construct the other 7 octants by symmetry operations.
   */
  for(i=0; (i<=nx/2); i++) {
    iip = (nx-i) % nx;
    for(j=0; (j<=ny/2); j++) {
      jjp = (ny-j) % ny;
      for(k=0; (k<=nz/2); k++) {
	kkp = (nz-k) % nz;
	ggg                 = ghat[i][j][k];
	ghat[i]  [jjp][k]   = ggg;
	ghat[i]  [j]  [kkp] = ggg;
	ghat[i]  [jjp][kkp] = ggg;
	ghat[iip][j]  [k]   = ggg;
	ghat[iip][jjp][k]   = ggg;
	ghat[iip][j]  [kkp] = ggg;
	ghat[iip][jjp][kkp] = ggg;
      }
    }
  }
}

void mk_ghat(FILE *fp,int nx,int ny,int nz,real ***ghat,
	     rvec box,real r1,real rc,gmx_bool bSym,gmx_bool bOld)
{
  int  ix,iy,iz;
  int  ixmax,iymax,izmax;
  real k2,ggg;
  rvec k,lll;
  
  calc_lll(box,lll);
    
  if (bSym) {
    ixmax=nx/2+1;
    iymax=ny/2+1;
    izmax=nz/2+1;
  }
  else {
    ixmax=nx;
    iymax=ny;
    izmax=nz;
  }
  /* Loop over lattice vectors in fourier space */    
  for(ix=0; (ix < ixmax); ix++) {
    for(iy=0; (iy < iymax); iy++) {
      for(iz=0; (iz < izmax); iz++) {
	calc_k(lll,ix,iy,iz,nx,ny,nz,k);
	k2 = iprod(k,k);
	if ((ix == 0) && (iy == 0) && (iz == 0))
	  ggg = 0.0;
	else {
	  if (bOld)
	    ggg = gk(sqrt(k2),rc,r1)/(k2*EPSILON0);
	  else
	    ggg = gknew(sqrt(k2),rc,r1)/(k2*EPSILON0);
	}
	ghat[ix][iy][iz]=ggg;
      }
    }
  }
  if (bSym)
    symmetrize_ghat(nx,ny,nz,ghat);
}

void wr_ghat(const char *fn,const output_env_t oenv,  
	     int n1max,int n2max,int n3max,real h1,real h2,real h3,
	     real ***ghat,int nalias,int porder,int niter,gmx_bool bSym,rvec beta,
	     real r1,real rc,real pval,real zval,real eref,real qopt)
{
  FILE *fp;
  int  N1MAX,N2MAX,N3MAX;
  gmx_bool bNL=FALSE;
  real rx,ry,rz;
  int  ii,jj,kk,nn;
  
  fp = gmx_fio_fopen(fn,"w");
  fprintf(fp,"%8d  %8d  %8d  %15.10e  %15.10e %15.10e\n",
	  n1max,n2max,n3max,h1,h2,h3);
  fprintf(fp,"%8d  %8d  %8d  %8d  %15.10e  %15.10e  %15.10e\n",
	  nalias,porder,niter,bSym,beta[XX],beta[YY],beta[ZZ]);
  fprintf(fp,"%10g  %10g  %10g  %10g  %10g  %10g\n",
	  rc,r1,pval,zval,eref,qopt);
  
  if (bSym) {
    N1MAX = n1max/2+1;
    N2MAX = n2max/2+1;
    N3MAX = n3max/2+1;
  }
  else {
    N1MAX = n1max;
    N2MAX = n2max;
    N3MAX = n3max;
  }
  for(ii=0; (ii<N1MAX); ii++) {
    for(jj=0; (jj<N2MAX); jj++) {
      for(kk=0,nn=1; (kk<N3MAX); kk++,nn++) { 
	bNL=FALSE;
	fprintf(fp,"  %12.5e",ghat[ii][jj][kk]);
	if ((nn % 6) == 0) {
	  fprintf(fp,"\n");
	  bNL=TRUE;
	}
      }
      if (!bNL)
	fprintf(fp,"\n");
    }
  }
  gmx_fio_fclose(fp);
  
  fp=xvgropen("ghat.xvg","G-Hat","k","gk",oenv);
  for(ii=0; (ii<N1MAX); ii++) {
    rx=sqr((real)(ii*h1));
    for(jj=0; (jj<N2MAX); jj++) {
      ry=rx+sqr((real)(jj*h2));
      for(kk=0; (kk<N3MAX); kk++) {
	rz=ry+sqr((real)(kk*h3));
	fprintf(fp,"%10g  %10g\n",sqrt(rz),ghat[ii][jj][kk]);
      }
    }
  }
  gmx_fio_fclose(fp);
}

void pr_scalar_gk(const char *fn,const output_env_t oenv,int nx,int ny,int nz,
                  rvec box,real ***ghat)
{
  FILE *fp;
  int  ii,jj,kk;
  real k1;
  rvec k,lll;
  
  calc_lll(box,lll);
  
  fp=xvgropen(fn,"G-Hat","k","gk",oenv);
  for(ii=0; (ii<nx); ii++) {
    for(jj=0; (jj<ny); jj++) {
      for(kk=0; (kk<nz); kk++) {
	calc_k(lll,ii,jj,kk,nx,ny,nz,k);
	k1 = norm(k);
	fprintf(fp,"%10g  %10g\n",k1,ghat[ii][jj][kk]);
      }
    }
  }
  gmx_fio_fclose(fp);
}

static real ***mk_rgrid(int nx,int ny,int nz)
{
  real *ptr1;
  real **ptr2;
  real ***ptr3;
  int  i,j,n2,n3;

  snew(ptr1,nx*ny*nz);
  snew(ptr2,nx*ny);
  snew(ptr3,nx);

  n2=n3=0;
  for(i=0; (i<nx); i++) {
    ptr3[i]=&(ptr2[n2]);
    for(j=0; (j<ny); j++,n2++) {
      ptr2[n2] = &(ptr1[n3]);
      n3 += nz;
    }
  }
  return ptr3;
}

real ***rd_ghat(FILE *log,const output_env_t oenv,char *fn,ivec igrid,
                rvec gridspace, rvec beta,int *porder,real *r1,real *rc)
{
  FILE   *in;
  real   ***gh;
  double gx,gy,gz,alX,alY,alZ,ddd;
  double acut,pval,zval,eref,qopt,r11;
  int    nalias,niter,bSym;
  int    ix,iy,iz,ixmax,iymax,izmax;
  
  in=gmx_fio_fopen(fn,"r");
  if(6 != fscanf(in,"%d%d%d%lf%lf%lf",&ix,&iy,&iz,&gx,&gy,&gz))
  {
      gmx_fatal(FARGS,"Error reading from file %s",fn);
  }


  igrid[XX]=ix, igrid[YY]=iy, igrid[ZZ]=iz;
  gridspace[XX]=gx,  gridspace[YY]=gy,  gridspace[ZZ]=gz;
  if(7 != fscanf(in,"%d%d%d%d%lf%lf%lf",&nalias,porder,&niter,&bSym,&alX,&alY,&alZ))
  {
      gmx_fatal(FARGS,"Error reading from file %s",fn);
  }

  if(6 != fscanf(in,"%lf%lf%lf%lf%lf%lf",&acut,&r11,&pval,&zval,&eref,&qopt))
  {
    gmx_fatal(FARGS,"Error reading from file %s",fn);
  }

  *r1 = r11;
  *rc = acut;
  
  fprintf(log,"\nOpened %s for reading ghat function\n",fn);
  fprintf(log,"gridsize: %10d %10d %10d\n",ix,iy,iz);
  fprintf(log,"spacing:  %10g %10g %10g\n",gx,gy,gz);
  fprintf(log,"    nalias    porder     niter      bSym      beta[X-Z]\n"
	  "%10d%10d%10d%10d%10g%10g%10g\n",
	  nalias,*porder,niter,bSym,alX,alY,alZ);
  fprintf(log,"      acut        r1      pval      zval      eref      qopt\n"
	  "%10g%10g%10g%10g%10g%10g\n",acut,*r1,pval,zval,eref,qopt);
  fflush(log);
  
  beta[XX] = alX;
  beta[YY] = alY;
  beta[ZZ] = alZ;
  
  gh      = mk_rgrid(ix,iy,iz);
  if (bSym) {
    ixmax=igrid[XX]/2+1;
    iymax=igrid[YY]/2+1;
    izmax=igrid[ZZ]/2+1;
  }
  else {
    ixmax=igrid[XX];
    iymax=igrid[YY];
    izmax=igrid[ZZ];
  }
  fprintf(log,"Reading ghat of %d %d %d\n",ixmax,iymax,izmax);
  for(ix=0; (ix<ixmax); ix++)
    for(iy=0; (iy<iymax); iy++)
      for(iz=0; (iz<izmax); iz++) {
	if( 1 != fscanf(in,"%lf",&ddd))
        {
	    gmx_fatal(FARGS,"Error reading from file %s",fn);
	}

	gh[ix][iy][iz] = ddd;
      }
  gmx_fio_fclose(in);

  wr_ghat("output.hat",oenv,igrid[XX],igrid[YY],igrid[ZZ],gx,gy,gz,gh,
	  nalias,*porder,niter,bSym,beta,
	  *r1,*rc,pval,zval,eref,qopt);
    
  if (bSym) 
    symmetrize_ghat(igrid[XX],igrid[YY],igrid[ZZ],gh);
  
  fprintf(log,"\nSuccessfully read ghat function!\n");
  
  
  return gh;
}

