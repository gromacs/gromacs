/*
 * $Id$
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

#include <math.h>
#include "main.h"
#include "constr.h"
#include "physics.h"
#include "vec.h"
#include "pbc.h"

void clincsp(rvec *x,rvec *f,rvec *fp,t_pbc *pbc,int ncons,
	     int *bla1,int *bla2,int *blnr,int *blbnb,
	     real *blc,real *blcc,real *blm,
	     int nrec,real *invmass,rvec *r,
	     real *rhs1,real *rhs2,real *sol)
{
  int     b,i,j,k,n,it,rec;
  real    tmp0,tmp1,tmp2,im1,im2,mvb,rlen,len,wfac,lam;  
  real    u0,u1,u2,v0,v1,v2;
  rvec    dx;
  real    *tmp;

  /* Compute normalized i-j vectors */
  if (pbc) {
    for(b=0;b<ncons;b++) {
      pbc_dx(pbc,x[bla1[b]],x[bla2[b]],dx);
      unitv(dx,r[b]);
    }
  } else {
    for(b=0;b<ncons;b++) {
      rvec_sub(x[bla1[b]],x[bla2[b]],dx);
      unitv(dx,r[b]);
    } /* 16 ncons flops */
  }
  
  for(b=0;b<ncons;b++) {
    tmp0=r[b][0];
    tmp1=r[b][1];
    tmp2=r[b][2];
    i=bla1[b];
    j=bla2[b];
    for(n=blnr[b];n<blnr[b+1];n++) {
      k=blbnb[n];
      blm[n]=blcc[n]*(tmp0*r[k][0]+tmp1*r[k][1]+tmp2*r[k][2]); 
    } /* 6 nr flops */
    mvb=blc[b]*(tmp0*(f[i][0]-f[j][0])+
		tmp1*(f[i][1]-f[j][1])+    
		tmp2*(f[i][2]-f[j][2]));
    rhs1[b]=mvb;
    sol[b]=mvb;
    /* 7 flops */
  }
  /* Together: 23*ncons + 6*nrtot flops */
    
  for(rec=0;rec<nrec;rec++) {
    for(b=0;b<ncons;b++) {
      mvb=0;
      for(n=blnr[b];n<blnr[b+1];n++) {
	j=blbnb[n];
	mvb=mvb+blm[n]*rhs1[j];
      }
      rhs2[b]=mvb;
      sol[b]=sol[b]+mvb;
    }
    tmp=rhs1;
    rhs1=rhs2;
    rhs2=tmp;
  } /* nrec*(ncons+2*nrtot) flops */
  
  for(b=0;b<ncons;b++) {
    i=bla1[b];
    j=bla2[b];
    mvb=blc[b]*sol[b];
    im1=invmass[i];
    im2=invmass[j];
    tmp0=r[b][0]*mvb;
    tmp1=r[b][1]*mvb;
    tmp2=r[b][2]*mvb;
    u0=fp[i][0]-tmp0*im1;
    u1=fp[i][1]-tmp1*im1;
    u2=fp[i][2]-tmp2*im1;
    v0=fp[j][0]+tmp0*im2;
    v1=fp[j][1]+tmp1*im2;
    v2=fp[j][2]+tmp2*im2;
    fp[i][0]=u0;
    fp[i][1]=u1;
    fp[i][2]=u2;
    fp[j][0]=v0;
    fp[j][1]=v1;
    fp[j][2]=v2;
  } /* 16 ncons flops */
}

void clincs(rvec *x,rvec *xp,t_pbc *pbc,int ncons,
	    int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,
	    real *blc,real *blcc,real *blm,
	    int nit,int nrec,real *invmass,rvec *r,
	    real *rhs1,real *rhs2,real *sol,real wangle,int *warn,
	    real *lambda)
{
  int     b,i,j,k,n,it,rec;
  real    tmp0,tmp1,tmp2,im1,im2,mvb,rlen,len,wfac,lam;  
  real    u0,u1,u2,v0,v1,v2;
  rvec    dx;
  real    *tmp;

  *warn=0;

  if (pbc) {
    /* Compute normalized i-j vectors */
    for(b=0; b<ncons; b++) {
      pbc_dx(pbc,x[bla1[b]],x[bla2[b]],dx);
      unitv(dx,r[b]);
    }  
    for(b=0; b<ncons; b++) {
      for(n=blnr[b]; n<blnr[b+1]; n++) {
	blm[n] = blcc[n]*iprod(r[b],r[blbnb[n]]);
      }
      pbc_dx(pbc,xp[bla1[b]],xp[bla2[b]],dx);
      mvb = blc[b]*(iprod(r[b],dx) - bllen[b]);
      rhs1[b] = mvb;
      sol[b]  = mvb;
    }
  } else {
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
    } /* 16 ncons flops */
    
    for(b=0;b<ncons;b++) {
      tmp0=r[b][0];
      tmp1=r[b][1];
      tmp2=r[b][2];
      len=bllen[b];
      i=bla1[b];
      j=bla2[b];
      for(n=blnr[b];n<blnr[b+1];n++) {
	k=blbnb[n];
	blm[n]=blcc[n]*(tmp0*r[k][0]+tmp1*r[k][1]+tmp2*r[k][2]); 
      } /* 6 nr flops */
      mvb=blc[b]*(tmp0*(xp[i][0]-xp[j][0])+
		  tmp1*(xp[i][1]-xp[j][1])+    
		  tmp2*(xp[i][2]-xp[j][2])-len);
      rhs1[b]=mvb;
      sol[b]=mvb;
      /* 8 flops */
    }
    /* Together: 24*ncons + 6*nrtot flops */
  }
    
  for(rec=0;rec<nrec;rec++) {
    for(b=0;b<ncons;b++) {
      mvb=0;
      for(n=blnr[b];n<blnr[b+1];n++) {
	j=blbnb[n];
	mvb=mvb+blm[n]*rhs1[j];
      }
      rhs2[b]=mvb;
      sol[b]=sol[b]+mvb;
    }
    tmp=rhs1;
    rhs1=rhs2;
    rhs2=tmp;
  } /* nrec*(ncons+2*nrtot) flops */
  
  for(b=0;b<ncons;b++) {
    i=bla1[b];
    j=bla2[b];
    mvb=blc[b]*sol[b];
    lambda[b]=-mvb;
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
  } /* 16 ncons flops */
  
  
  
  /*     
  ********  Correction for centripetal effects  ********  
  */
  
  wfac=cos(DEG2RAD*wangle);
  wfac=wfac*wfac;

  for(it=0; it<nit; it++) {
  
    for(b=0;b<ncons;b++) {
      len = bllen[b];
      if (pbc) {
	pbc_dx(pbc,xp[bla1[b]],xp[bla2[b]],dx);
      } else {
	rvec_sub(xp[bla1[b]],xp[bla2[b]],dx);
      }
      u1 = len*len;
      u0 = 2*u1 - norm2(dx);
      if (u0 < wfac*u1) *warn = b;	
      if (u0 < 0) u0 = 0;
      mvb = blc[b]*(len - u0*invsqrt(u0));
      rhs1[b] = mvb;
      sol[b]  = mvb;
    } /* 18*ncons flops */
    
    for(rec=0;rec<nrec;rec++) {
      for(b=0;b<ncons;b++) {
	mvb=0;
	for(n=blnr[b];n<blnr[b+1];n++) {
	  j=blbnb[n];
	  mvb=mvb+blm[n]*rhs1[j];
	}
	rhs2[b]=mvb;
	sol[b]=sol[b]+mvb;
      }
      tmp=rhs1;
      rhs1=rhs2;
      rhs2=tmp;
    } /* nrec*(ncons+2*nrtot) flops */ 
    
    for(b=0;b<ncons;b++) {
      i=bla1[b];
      j=bla2[b];
      lam=lambda[b];
      mvb=blc[b]*sol[b];
      lambda[b]=lam-mvb;
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
    } /* 17 ncons flops */
  } /* nit*ncons*(35+9*nrec) flops */
  /* Total:
   * 24*ncons + 6*nrtot + nrec*(ncons+2*nrtot)
   * + nit * (18*ncons + nrec*(ncons+2*nrtot) + 17 ncons)
   *
   * (24+nrec)*ncons + (6+2*nrec)*nrtot
   * + nit * ((35+nrec)*ncons + 2*nrec*nrtot)
   * if nit=1
   * (59+nrec)*ncons + (6+4*nrec)*nrtot
   */
}
