/*
 *       $Id$
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
static char *SRCID_clincs_c = "$Id$";

#include <math.h>
#include "main.h"
#include "constr.h"
#include "physics.h"
#include "vec.h"

void clincs(rvec *x,rvec *xp,int ncons,
	    int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,
	    real *blc,real *blcc,real *blm,
	    int nit,int nrec,real *invmass,rvec * r,
	    real *rhs1,real *rhs2,real *sol,real wangle,int *warn,
	    real *lambda)
{
  int     b,i,j,k,n,it,rec;
  real    tmp0,tmp1,tmp2,im1,im2,mvb,rlen,len,wfac,lam;  
  real    u0,u1,u2,v0,v1,v2;
  real    *tmp;

  *warn=0;

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

void cconerr(real *max,real *rms,int *imax,rvec *xprime,
	     int ncons,int *bla1,int *bla2,real *bllen)
     
{
  real      len,d,ma,ms,tmp0,tmp1,tmp2;
  int       b,i,j,im;
  
  ma=0;
  ms=0;
  im=0;
  for(b=0;b<ncons;b++) {
    i=bla1[b];
    j=bla2[b];
    tmp0=xprime[i][0]-xprime[j][0];
    tmp1=xprime[i][1]-xprime[j][1];
    tmp2=xprime[i][2]-xprime[j][2];
    len=sqrt(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2);
    d=fabs(len/bllen[b]-1);
    if (d > ma) {
      ma=d;
      im=b;
    }
    ms=ms+d*d;
  }
  *max=ma;
  *rms=sqrt(ms/ncons);
  *imax=im;
}

void lincs_warning(rvec *x,rvec *xprime,
		   int ncons,int *bla1,int *bla2,real *bllen,real wangle)
{
  int b,i,j;
  rvec v0,v1;
  real wfac,d0,d1,cosine;
  char buf[STRLEN];
  
  wfac=cos(DEG2RAD*wangle);
  
  sprintf(buf,"bonds that rotated more than %g degrees:\n"
	  " atom 1 atom 2  angle  previous, current, constraint length\n",
	  wangle);
  fprintf(stderr,buf);
  fprintf(stdlog,buf); 

  for(b=0;b<ncons;b++) {
    i=bla1[b];
    j=bla2[b];
    rvec_sub(x[i],x[j],v0);
    rvec_sub(xprime[i],xprime[j],v1);
    d0=norm(v0);
    d1=norm(v1);
    cosine=iprod(v0,v1)/(d0*d1);
    if (cosine<wfac) {
      sprintf(buf," %6d %6d  %5.1f  %8.4f %8.4f    %8.4f\n",
	      i+1,j+1,RAD2DEG*acos(cosine),d0,d1,bllen[b]);
      fprintf(stderr,buf);
      fprintf(stdlog,buf);
    }
  }
}


