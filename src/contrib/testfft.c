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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#include <math.h>
#include <stdio.h>
#include "typedefs.h"
#include "macros.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/xvgr.h"
#include "complex.h"
#include "fftgrid.h"
#include "mdrun.h"

void testfft(FILE *fp,t_complex ***grid,int nx,int ny,int nz,gmx_bool bFirst)
{
#ifdef USE_SGI_FFT
#ifdef GMX_DOUBLE
  static    zomplex   *coeff;
#else
  static    complex   *coeff;
#endif  
  static    int la1,la2;
#endif
  t_complex *cptr;
  real      *gptr,*fqqq,fg,fac;
  int       ntot,i,j,k,m,n,ndim[4];
  int       npppm;
  
  ndim[0] = 0;
  ndim[1] = nx;
  ndim[2] = ny;
  ndim[3] = nz;
  
  ntot    = nx*ny*nz;
  cptr    = grid[0][0];
  fqqq    = &(grid[0][0][0].re);
  
#ifdef USE_SGI_FFT
  if (bFirst) {
    fprintf(fp,"Going to use SGI optimized FFT routines.\n");
#ifdef GMX_DOUBLE
    coeff  = zfft3di(nx,ny,nz,NULL);
#else
    coeff  = cfft3di(nx,ny,nz,NULL);
#endif
    bFirst = FALSE;
  }
  la1 = nx;
  la2 = ny;
#ifdef GMX_DOUBLE
  zfft3d(1,nx,ny,nz,(zomplex *)cptr,la1,la2,coeff);
#else
  cfft3d(1,nx,ny,nz,(complex *)cptr,la1,la2,coeff);
#endif
#else
  fourn(fqqq-1,ndim,3,1);
#endif
  
#ifdef USE_SGI_FFT
#ifdef GMX_DOUBLE
  zfft3d(-1,nx,ny,nz,(zomplex *)cptr,la1,la2,coeff);
#else
  cfft3d(-1,nx,ny,nz,(complex *)cptr,la1,la2,coeff);
#endif
#else
  fourn(fqqq-1,ndim,3,-1);
#endif
}

void testrft(FILE *fp,real ***grid,int nx,int ny,int nz,gmx_bool bFirst)
{
#ifdef USE_SGI_FFT
#ifdef GMX_DOUBLE
  static    double *coeff;
#else
  static    float *coeff;
#endif
  static    int la1,la2;
#endif
  real      *cptr;
  real      *gptr,*fqqq,fg,fac;
  int       ntot,i,j,k,m,n,ndim[4];
  int       npppm,job;
  
  ndim[0] = 0;
  ndim[1] = nx;
  ndim[2] = ny;
  ndim[3] = nz;
  
  ntot    = nx*ny*nz;
  cptr    = grid[0][0];
  fqqq    = &(grid[0][0][0]);
  
#ifdef USE_SGI_FFT
  if (bFirst) {
    fprintf(fp,"Going to use SGI optimized FFT routines.\n");
#ifdef GMX_DOUBLE
    coeff  = dfft3di(nx,ny,nz,NULL);
#else
    coeff  = sfft3di(nx,ny,nz,NULL);
#endif
    bFirst = FALSE;
  }
  job = 1;
  la1 = nx+2;
  la2 = ny;
#ifdef GMX_DOUBLE
  dzfft3d(job,nx,ny,nz,cptr,la1,la2,coeff);
#else
  scfft3d(job,nx,ny,nz,cptr,la1,la2,coeff);
#endif
#else
  fourn(fqqq-1,ndim,3,1);
#endif

  job = -1;
  
#ifdef USE_SGI_FFT
#ifdef GMX_DOUBLE
  zdfft3d(job,nx,ny,nz,cptr,la1,la2,coeff);
#else
  csfft3d(job,nx,ny,nz,cptr,la1,la2,coeff);
#endif
#else
  fourn(fqqq-1,ndim,3,-1);
#endif
}

int main(int argc,char *argv[])
{
  FILE      *fp;
  int       nnn[] = { 8, 10, 12, 15, 16, 18, 20, 24, 25, 27, 30, 32, 36, 40,
		      45, 48, 50, 54, 60, 64, 72, 75, 80, 81, 90, 100 };
#define NNN asize(nnn)
  int       *niter;
  int       i,j,n,nit,ntot,n3;
  double    t,nflop;
  double    *rt,*ct;
  t_complex ***g;
  real      ***h;
  
  snew(rt,NNN);
  snew(ct,NNN);
  snew(niter,NNN);
  
  for(i=0; (i<NNN); i++) {
    n = nnn[i];
    fprintf(stderr,"\rReal %d     ",n);
    if (n < 16)
      niter[i] = 100;
    else if (n < 26)
      niter[i] = 50;
    else if (n < 51)
      niter[i] = 10;
    else
      niter[i] = 5;
    nit = niter[i];
      
    h   = mk_rgrid(n+2,n,n);
    start_time();
    for(j=0; (j<nit); j++) {
      testrft(stdout,h,n,n,n,(j==0));
    }
    update_time();
    rt[i] = node_time();
    free_rgrid(h,n,n);
    
    fprintf(stderr,"\rComplex %d     ",n);
    g   = mk_cgrid(n,n,n);
    start_time();
    for(j=0; (j<nit); j++) {
      testfft(stdout,g,n,n,n,(j==0));
    }
    update_time();
    ct[i] = node_time();
    free_cgrid(g,n,n);
  }
  fprintf(stderr,"\n");
  fp=xvgropen("timing.xvg","FFT timings per grid point","n","t (s)");
  for(i=0; (i<NNN); i++) {
    n3 = 2*niter[i]*nnn[i]*nnn[i]*nnn[i];
    fprintf(fp,"%10d  %10g  %10g\n",nnn[i],rt[i]/n3,ct[i]/n3);
  }
  gmx_fio_fclose(fp);
  
  return 0;
}
