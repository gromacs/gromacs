/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_testfft_c = "$Id$";

#include <math.h>
#include <stdio.h>
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "xvgr.h"
#include "complex.h"
#include "grids.h"
#include "mdrun.h"
#ifdef USE_SGI_FFT
#include <fft.h>
#endif

void testfft(FILE *fp,t_complex ***grid,int nx,int ny,int nz,bool bFirst)
{
#ifdef USE_SGI_FFT
#ifdef DOUBLE
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
#ifdef DOUBLE
    coeff  = zfft3di(nx,ny,nz,NULL);
#else
    coeff  = cfft3di(nx,ny,nz,NULL);
#endif
    bFirst = FALSE;
  }
  la1 = nx;
  la2 = ny;
#ifdef DOUBLE
  zfft3d(1,nx,ny,nz,(zomplex *)cptr,la1,la2,coeff);
#else
  cfft3d(1,nx,ny,nz,(complex *)cptr,la1,la2,coeff);
#endif
#else
  fourn(fqqq-1,ndim,3,1);
#endif
  
#ifdef USE_SGI_FFT
#ifdef DOUBLE
  zfft3d(-1,nx,ny,nz,(zomplex *)cptr,la1,la2,coeff);
#else
  cfft3d(-1,nx,ny,nz,(complex *)cptr,la1,la2,coeff);
#endif
#else
  fourn(fqqq-1,ndim,3,-1);
#endif
}

void testrft(FILE *fp,real ***grid,int nx,int ny,int nz,bool bFirst)
{
#ifdef USE_SGI_FFT
#ifdef DOUBLE
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
#ifdef DOUBLE
    coeff  = dfft3di(nx,ny,nz,NULL);
#else
    coeff  = sfft3di(nx,ny,nz,NULL);
#endif
    bFirst = FALSE;
  }
  job = 1;
  la1 = nx+2;
  la2 = ny;
#ifdef DOUBLE
  dzfft3d(job,nx,ny,nz,cptr,la1,la2,coeff);
#else
  scfft3d(job,nx,ny,nz,cptr,la1,la2,coeff);
#endif
#else
  fourn(fqqq-1,ndim,3,1);
#endif

  job = -1;
  
#ifdef USE_SGI_FFT
#ifdef DOUBLE
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
    rt[i] = cpu_time();
    free_rgrid(h,n,n);
    
    fprintf(stderr,"\rComplex %d     ",n);
    g   = mk_cgrid(n,n,n);
    start_time();
    for(j=0; (j<nit); j++) {
      testfft(stdout,g,n,n,n,(j==0));
    }
    update_time();
    ct[i] = cpu_time();
    free_cgrid(g,n,n);
  }
  fprintf(stderr,"\n");
  fp=xvgropen("timing.xvg","FFT timings per grid point","n","t (s)");
  for(i=0; (i<NNN); i++) {
    n3 = 2*niter[i]*nnn[i]*nnn[i]*nnn[i];
    fprintf(fp,"%10d  %10g  %10g\n",nnn[i],rt[i]/n3,ct[i]/n3);
  }
  fclose(fp);
  
  return 0;
}
