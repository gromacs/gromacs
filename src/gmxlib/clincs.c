/*
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
 * S  C  A  M  O  R  G
 */

#include <math.h>
#include "main.h"
#include "update.h"
#include "physics.h"
#include "vec.h"

void clincs(rvec *x,rvec *xp,int ncons,int ncm,int cmax,
	    int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,
	    real *blc,real *blcc,real *blm,
	    int nrec,real *invmass,rvec * r,
	    real *rhs1,real *rhs2,real *sol,real wangle,int *warn,
	    real *lambda)
{
  int     b,i,j,k,n,b4,rec,nr,n1,nc4;
  real    tmp0,tmp1,tmp2,im1,im2,mvb,rlen,len,wfac,lam;  
  real    u0,u1,u2,v0,v1,v2;
  real    *tmp;

  *warn=0;
  n1=ncons-ncm;
  nc4=(cmax-4)*n1;

  /* Compute normalized i-j vectors */
  for(b=0;b<ncons;b++) {
    i=bla1[b];
    j=bla2[b];
    tmp0=x[i][0]-x[j][0];
    tmp1=x[i][1]-x[j][1];
    tmp2=x[i][2]-x[j][2];
    rlen=invsqrt(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2);
    r[b][0]=rlen*tmp0;
    r[b][1]=rlen*tmp1;
    r[b][2]=rlen*tmp2;
  }
  
  for(b=0;b<n1;b++) {
    b4=4*b;
    tmp0=r[b][0];
    tmp1=r[b][1];
    tmp2=r[b][2];
    len=bllen[b];
    i=bla1[b];
    j=bla2[b];
    nr=blnr[b];
    for(n=0;n<nr;n++) {
      k=blbnb[b4+n];
      blm[b4+n]=blcc[b4+n]*(tmp0*r[k][0]+tmp1*r[k][1]+tmp2*r[k][2]);
    }
    mvb=blc[b]*(tmp0*(xp[i][0]-xp[j][0])+
		tmp1*(xp[i][1]-xp[j][1])+    
		tmp2*(xp[i][2]-xp[j][2])-len);
    rhs1[b]=mvb;
    sol[b]=mvb;
  }
  
  for(b=n1;b<ncons;b++) {
    b4=cmax*b-nc4;
    tmp0=r[b][0];
    tmp1=r[b][1];
    tmp2=r[b][2];
    len=bllen[b];
    i=bla1[b];
    j=bla2[b];
    nr=blnr[b];
    for(n=0;n<nr;n++) {
      k=blbnb[b4+n];
      blm[b4+n]=blcc[b4+n]*(tmp0*r[k][0]+tmp1*r[k][1]+tmp2*r[k][2]); 
    }
    mvb=blc[b]*(tmp0*(xp[i][0]-xp[j][0])+
		tmp1*(xp[i][1]-xp[j][1])+    
		tmp2*(xp[i][2]-xp[j][2])-len);
    rhs1[b]=mvb;
    sol[b]=mvb;
  }
  
  
  for(rec=0;rec<nrec;rec++) {
    for(b=0;b<n1;b++) {
      b4=4*b;
      mvb=0;
      for(n=0;n<4;n++) {
	j=blbnb[b4+n];
	mvb=mvb+blm[b4+n]*rhs1[j];
      }
      rhs2[b]=mvb;
      sol[b]=sol[b]+mvb;
    }
    for(b=n1;b<ncons;b++) {
      b4=cmax*b-nc4;
      mvb=0;
      nr=blnr[b];
      for(n=0;n<nr;n++) {
	j=blbnb[b4+n];
	mvb=mvb+blm[b4+n]*rhs1[j];
      }
      rhs2[b]=mvb;
      sol[b]=sol[b]+mvb;
    }
    tmp=rhs1;
    rhs1=rhs2;
    rhs2=tmp;
  }
  
  for(b=0;b<ncons;b++) {
    i=bla1[b];
    j=bla2[b];
    mvb=blc[b]*sol[b];
    lambda[b]=mvb;
    im1=invmass[i];
    im2=invmass[j];
    tmp0=r[b][0]*mvb;
    tmp1=r[b][1]*mvb;
    tmp2=r[b][2]*mvb;
    u0=xp[i][0]-tmp0*im1;
    u1=xp[i][1]-tmp1*im1;
    u2=xp[i][2]-tmp2*im1;
    v0=xp[j][0]+tmp0*im2;
    v1=xp[j][1]+tmp1*im2;
    v2=xp[j][2]+tmp2*im2;
    xp[i][0]=u0;
    xp[i][1]=u1;
    xp[i][2]=u2;
    xp[j][0]=v0;
    xp[j][1]=v1;
    xp[j][2]=v2;
  }
  
  
  
  /*     
  ********  Correction for centripetal effects  ********  
  */
  
  wfac=cos(DEG2RAD*wangle);
  wfac=wfac*wfac;
  
  for(b=0;b<ncons;b++) {
    len=bllen[b];
    i=bla1[b];
    j=bla2[b];
    tmp0=xp[i][0]-xp[j][0];
    tmp1=xp[i][1]-xp[j][1];
    tmp2=xp[i][2]-xp[j][2];
    u1=len*len;
    u0=2.*u1-(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2);
    if (u0 < wfac*u1) *warn=b;  
    if (u0 < 0) u0=0;
    mvb=blc[b]*(len-sqrt(u0));
    rhs1[b]=mvb;
    sol[b]=mvb;
  }
  
  for(rec=0;rec<nrec;rec++) {
    for(b=0;b<n1;b++) {
      b4=4*b;
      mvb=0;
      for(n=0;n<4;n++) {
	j=blbnb[b4+n];
	mvb=mvb+blm[b4+n]*rhs1[j];
      }
      rhs2[b]=mvb;
      sol[b]=sol[b]+mvb;
    }
    for(b=n1;b<ncons;b++) {
      b4=cmax*b-nc4;
      mvb=0;
      nr=blnr[b];
      for(n=0;n<nr;n++) {
	j=blbnb[b4+n];
	mvb=mvb+blm[b4+n]*rhs1[j];
      }
      rhs2[b]=mvb;
      sol[b]=sol[b]+mvb;
    }
    tmp=rhs1;
    rhs1=rhs2;
    rhs2=tmp;
  }
  
  for(b=0;b<ncons;b++) {
    i=bla1[b];
    j=bla2[b];
    lam=lambda[b];
    mvb=blc[b]*sol[b];
    lambda[b]=lam+mvb;
    im1=invmass[i];
    im2=invmass[j];
    tmp0=r[b][0]*mvb;
    tmp1=r[b][1]*mvb;
    tmp2=r[b][2]*mvb;
    u0=xp[i][0]-tmp0*im1;
    u1=xp[i][1]-tmp1*im1;
    u2=xp[i][2]-tmp2*im1;
    v0=xp[j][0]+tmp0*im2;
    v1=xp[j][1]+tmp1*im2;
    v2=xp[j][2]+tmp2*im2;
    xp[i][0]=u0;
    xp[i][1]=u1;
    xp[i][2]=u2;
    xp[j][0]=v0;
    xp[j][1]=v1;
    xp[j][2]=v2;
  }
}



