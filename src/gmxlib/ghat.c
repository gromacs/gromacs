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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_ghat_c = "$Id$";

#include <stdio.h>
#include "typedefs.h"
#include "futil.h"
#include "vec.h"
#include "physics.h"
#include "lrutil.h"
#include "grids.h"
#include "ghat.h"

void symmetrize_ghat(int nx,int ny,int nz,real ***ghat)
{
  int  i,j,k;
  int  iip,jjp,kkp;
  real ggg;

  fprintf(stderr,"Symmetrizing Ghat function\n");  
  /* Only the lower octant of the rectangle has been saved,
   * so we must construct the other 7 octants by symmetry operations.
   */
  for(k=0; (k<=nz/2); k++) {
    for(j=0; (j<=ny/2); j++) {
      kkp = (nz-k) % nz;
      jjp = (ny-j) % ny;
      for(i=0; (i<=nx/2); i++) {
	iip = (nx-i) % nx;
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
	     rvec box,real r1,real rc,bool bSym)
{
  int  ix,iy,iz;
  int  ixmax,iymax,izmax;
  int  m;
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
	else 
	  ggg = gk(sqrt(k2),rc,r1)/(k2*EPSILON0);
	
	ghat[ix][iy][iz]=ggg;
      }
    }
  }
  if (bSym)
    symmetrize_ghat(nx,ny,ny,ghat);
}

real ***rd_ghat(FILE *log,char *fn,ivec igrid,rvec gridspace,
		real *beta,int *porder,real *rshort,real *rlong)
{
  FILE   *in;
  real   ***gh;
  double gx,gy,gz,alpha,ddd;
  double acut,r1,pval,zval,eref,qopt;
  int    nalias,niter,bSym;
  int    ix,iy,iz,jx,jy,jz,m,ixmax,iymax,izmax;
  
  in=ffopen(fn,"r");
  fscanf(in,"%d%d%d%lf%lf%lf",&ix,&iy,&iz,&gx,&gy,&gz);
  igrid[XX]=ix, igrid[YY]=iy, igrid[ZZ]=iz;
  gridspace[XX]=gx,  gridspace[YY]=gy,  gridspace[ZZ]=gz;
  fscanf(in,"%d%d%d%lf%d",&nalias,porder,&niter,&alpha,&bSym);
  fscanf(in,"%lf%lf%lf%lf%lf%lf",&acut,&r1,&pval,&zval,&eref,&qopt);
  
  fprintf(log,"\nOpening %s for reading ghat function\n",fn);
  fprintf(log,"gridsize: %10d %10d %10d\n",ix,iy,iz);
  fprintf(log,"spacing:  %10g %10g %10g\n",gx,gy,gz);
  fprintf(log,"    nalias    porder     niter      beta    bSym\n%10d%10d%10d%10g%10d\n",
	  nalias,*porder,niter,alpha,bSym);
  fprintf(log,"      acut        r1      pval      zval      eref      qopt\n"
	  "%10g%10g%10g%10g%10g%10g\n",acut,r1,pval,zval,eref,qopt);
  fflush(log);
  
  *beta   = alpha;
  *rshort = r1;
  *rlong  = acut;
  
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
	fscanf(in,"%lf",&ddd);
	gh[ix][iy][iz] = ddd;
      }
  ffclose(in);
  
  if (bSym) 
    symmetrize_ghat(igrid[XX],igrid[YY],igrid[ZZ],gh);
  
  fprintf(log,"\nSuccessfully read ghat function!\n");
  
  
  return gh;
}

