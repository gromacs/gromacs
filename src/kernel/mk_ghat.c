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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include "copyrite.h"
#include "macros.h"
#include "smalloc.h"
#include "typedefs.h"
#include "futil.h"
#include "xvgr.h"
#include "vec.h"
#include "maths.h"
#include "coulomb.h"
#include "physics.h"
#include "statutil.h"
#include "tpxio.h"
#include "copyrite.h"

const real tol = 1e-8;

real sinx_x(real x)
{
  if (x == 0)
    return 1.0;
  else
    return sin(x)/x;
}

real  uhat(int porder,real k1,real k2,real k3,real h1,real h2,real h3)
{
  real fac1,fac2,fac3;
  real f123,fff;
  int  i;
  
  fac1 = sinx_x(k1*h1*0.5);
  fac2 = sinx_x(k2*h2*0.5);
  fac3 = sinx_x(k3*h3*0.5);

  fff  = 1;
  f123 = fac1*fac2*fac3;
  for(i=1; (i<=porder+1); i++)
    fff *= f123;

  if (fabs(fff) < tol)
    return 0.0;
  else
    return fff;
}

real  uhat1D(int porder,real k1,real h1)
{
  real fac1;
  real fff;
  int  i;
  
  fac1 = sinx_x(k1*h1*0.5);

  fff  = 1;
  for(i=1; (i<=porder+1); i++)
    fff *= fac1;

  if (fabs(fff) < tol)
    return 0.0;
  else
    return fff;
}

real shat(real acut,real kmag,real r1)
{
  return gk(kmag,acut,r1);
}

real dhat(real alpha,real k,real h)
{
  real dh;

  dh = alpha*(sin(k*h)/h) + (1.0-alpha)*(sin(2.0*k*h)/(2.0*h) );
  
  if (fabs(dh) < tol)
    return 0.0;
  else
    return dh;
}

real cufrth(int porder,real k1,real k2,real k3,
	    real h1,real h2,real h3,int nalias)
{
  real kn1,kn2,kn3,tmp;
  int  n1,n2,n3;
  real ufrth;
  
  real twopi=2*M_PI;

  if (nalias == 0) {
    tmp   = uhat(porder,k1,k2,k3,h1,h2,h3);
    ufrth = tmp*tmp*tmp*tmp;
  }
  else {
    ufrth = 0.0;
    for( n1 = -nalias; (n1<= nalias); n1++) {
      for( n2 = -nalias; (n2<= nalias); n2++) {
	for( n3 = -nalias; (n3<= nalias); n3++) {
	  kn1 = k1 + n1*twopi/h1;
	  kn2 = k2 + n2*twopi/h2;
	  kn3 = k3 + n3*twopi/h3;
	  tmp = uhat(porder,kn1,kn2,kn3,h1,h2,h3);
	  ufrth += tmp*tmp*tmp*tmp;
	}
      }
    }
  }
  if (fabs(ufrth) < tol)
    return 0.0;
  else
    return ufrth;
}

real crsqal(real acut,real r1,real k1,real k2,real k3,
	    real h1,real h2,real h3,int nalias)
{
  real kn1,kn2,kn3;
  real h_1,h_2,h_3;
  real ksq,kmag,tmp,Rsqal;
  int  n1,n2,n3;

  real twopi=2*M_PI;

  h_1=twopi/h1;
  h_2=twopi/h2;
  h_3=twopi/h3;
  if (nalias==0) {
    ksq   = k1*k1 + k2*k2 + k3*k3;
    kmag  = sqrt(ksq);
    tmp   = shat(acut,kmag,r1);
    Rsqal = tmp*tmp/(ksq);
  }
  else {
    Rsqal = 0.0;
    for(n1 = -nalias; (n1<= nalias); n1++) {
      for( n2 = -nalias; (n2<= nalias); n2++) {
	for( n3 = -nalias; (n3<= nalias); n3++) {
	  kn1   = k1 + n1*h_1;
	  kn2   = k2 + n2*h_2;
	  kn3   = k3 + n3*h_3;
	  ksq   = kn1*kn1 + kn2*kn2 + kn3*kn3;
	  kmag  = sqrt(ksq);
	  tmp   = shat(acut,kmag,r1);
	  Rsqal = Rsqal + tmp*tmp/ksq;
	}
      }
    }
  }
  return Rsqal;
}


real usqsq(int porder,real k1,real k2,real k3,real h1,real h2,real h3)
{
  const real tt=2.0/3.0;
  const real tx=2.0/15.0;
  real t1,t2,t3,t12,t22,t32,tmp;
  
  t1 = sin(k1*h1*0.5);
  t2 = sin(k2*h2*0.5);
  t3 = sin(k3*h3*0.5);
  t12 = t1*t1;
  t22 = t2*t2;
  t32 = t3*t3;
  
  if (porder == 1) 
    tmp = (1.0-tt*t12)*(1-tt*t22)*(1-tt*t32);
  else if (porder == 2) 
    tmp = ( (1.0 - t12 + tx*t12*t12)*
	    (1.0 - t22 + tx*t22*t22)*
	    (1.0 - t32 + tx*t32*t32) );
  else 
    gmx_fatal(FARGS,"porder = %d in usqsq",porder);
  
  return tmp*tmp;
}

real usqsq1D(int porder,real k1,real h1)
{
  const real tt=2.0/3.0;
  const real tx=2.0/15.0;
  real t1,t12,tmp;
  
  t1 = sin(k1*h1*0.5);
  t12 = t1*t1;
  
  if (porder == 1) 
    tmp = (1.0-tt*t12);
  else if (porder == 2) 
    tmp = ( (1.0 - t12 + tx*t12*t12));
  else 
    gmx_fatal(FARGS,"porder = %d in usqsq",porder);
  
  return tmp*tmp;
}

real  ursum(int term,int porder,real acut,real r1,
	    real k1,real k2,real k3,real h1,real h2,real h3,int nalias)
{
  real kt,ksq,kmag;
  /*   real kcutsq; */
  real kn1,kn2,kn3,urs,tmp;
  real h_1,h_2,h_3;
  int  n1,n2,n3;

  real twopi=2*M_PI;
  h_1=twopi/h1;
  h_2=twopi/h2;
  h_3=twopi/h3;
  /*
    c
    c     for large enough values of k, the terms become negligable
    c     if shat(k) = exp(-k^2/4*acut) < eps
    c     kcutsq = 4*alpha* (-ln(eps))
    c     eps = 10^-6, -ln(eps) = 14
    c     eps = 10^-10, -ln(eps) = 23
    c     eps = 10^-20, -ln(eps) = 46
    c
    c     kcutsq = 4.0*acut*115;
  */

  if (nalias==0) {
    if (term==XX) kt = k1;
    if (term==YY) kt = k2;
    if (term==ZZ) kt = k3;
    ksq = k1*k1 + k2*k2 + k3*k3;
    kmag = sqrt(ksq);
    tmp = uhat(porder,k1,k2,k3,h1,h2,h3);
    urs = tmp*tmp*kt*shat(acut,kmag,r1)/(ksq);
  }
  else {
    urs = 0.0;
    for(n1 = -nalias; (n1<= nalias); n1++) {
      for( n2 = -nalias; (n2<= nalias); n2++) {
	for( n3 = -nalias; (n3<= nalias); n3++) {
	  kn1 = k1 + n1*h_1;
	  kn2 = k2 + n2*h_2;
	  kn3 = k3 + n3*h_3;
	  ksq = kn1*kn1 + kn2*kn2 + kn3*kn3;

	  if (term==XX) kt = kn1;
	  if (term==YY) kt = kn2;
	  if (term==ZZ) kt = kn3;
	  if (kt != 0.0) {
	    kmag = sqrt(ksq);
	    tmp = uhat(porder,kn1,kn2,kn3,h1,h2,h3);
	    if (tmp != 0.0)
	      urs = urs + tmp*tmp*kt*shat(acut,kmag,r1)/ksq;
	  }
	}
      }
    }
  }
  return urs;
}
 
real  ursum1D(int term,int porder,real acut,real r1,real k1,real h1,int nalias)
{
  real kt,ksq,kmag;
/*   real kcutsq; */
  real kn1,urs,tmp;
  real h_1;
  int  n1;

  real twopi=2*M_PI;
  h_1=twopi/h1;
  /*
    c
    c     for large enough values of k, the terms become negligable
    c     if shat(k) = exp(-k^2/4*acut) < eps
    c     kcutsq = 4*alpha* (-ln(eps))
    c     eps = 10^-6, -ln(eps) = 14
    c     eps = 10^-10, -ln(eps) = 23
    c     eps = 10^-20, -ln(eps) = 46
    c
    */
/*   kcutsq = 4.0*acut*115; */

  if (nalias==0) {
    if (term==1) kt = k1;
    ksq = k1*k1;
    kmag = sqrt(ksq);
    tmp = uhat1D(porder,k1,h1);
    urs = tmp*tmp*kt*shat(acut,kmag,r1)/(EPSILON0*ksq);
  }
  else {
    urs = 0.0;
    for(n1 = -nalias; (n1<= nalias); n1++) {
      kn1 = k1 + n1*h_1;
      ksq = kn1*kn1;
      /*c              if (ksq.lt.kcutsq) then*/
      if (term==XX) kt = kn1;
      if (kt != 0.0) {
	kmag = sqrt(ksq);
	tmp = uhat1D(porder,kn1,h1);
	if (tmp != 0.0)
	  urs = urs + tmp*tmp*kt*shat(acut,kmag,r1)/(EPSILON0*ksq);
      }
      /*c              endif*/
    }
  }
  return urs;
}
 
real sym(int indx,int maxind)
{
  if ( (indx == 0 ) || (indx == maxind/2) ) 
    return 1.0;
  else
    return 2.0;
}

void calc(gmx_bool bSym,gmx_bool bVerbose,
	  const int n1max,const int n2max,const int n3max,
	  const real h1,const real h2,const real h3,
	  int nalias,int porder,real acut,real r1,const real alpha,
	  const gmx_bool bSearch,
	  real ***ghat,real *ppval,real *zzval,real *eeref,real *qqopt)
{     
  real box1,box2,box3;
  real k1,k2,k3,ksq,kmag;
  real gnumer,dsq,gdenom,gsq;
  real ufrth,rsqal,rsq;
  real symfac;
  int  l1,l2,l3;
  real twopi=2*M_PI;
  real d1,d2,d3,u1,u2,u3,ss,gg;
  real pval,zval,eref,qopt;
  int  N1MAX,N2MAX,N3MAX;
  
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
    
  box1 = n1max*h1;
  box2 = n2max*h2;
  box3 = n3max*h3;

  pval = 0.0;
  zval = 0.0;
  eref = 0.0;
  qopt = 0.0;

  for(l1=0; (l1<N1MAX); l1++) {
    if (bVerbose)
      fprintf(stderr,"\rl1=%5d  qopt=%12.6e",l1,qopt);
      
    k1   = twopi*l1/box1;
    d1   = dhat(alpha,k1,h1);
    
    for(l2=0; (l2<N2MAX); l2++) {
      k2   = twopi*l2/box2;
      d2   = dhat(alpha,k2,h2);
      
      for(l3=0; (l3<N3MAX); l3++) {
	if (((l1+l2+l3) == 0) /*|| (l1 == n1max/2) || (l2 == n2max/2) || 
	    (l3 == n3max/2)*/)
	  ghat[0][0][0] = 0.0;
	else {
	  k3   = twopi*l3/box3;
	  ksq  = k1*k1 + k2*k2 + k3*k3;
	  kmag = sqrt(ksq);
	  
	  d3   = dhat(alpha,k3,h3);
	  u1   = ursum(XX,porder,acut,r1,k1,k2,k3,h1,h2,h3,nalias);
	  u2   = ursum(YY,porder,acut,r1,k1,k2,k3,h1,h2,h3,nalias);
	  u3   = ursum(ZZ,porder,acut,r1,k1,k2,k3,h1,h2,h3,nalias);
	  
	  gnumer = d1*u1+d2*u2+d3*u3;
	  dsq    = d1*d1+d2*d2+d3*d3;
	  gdenom = dsq*usqsq(porder,k1,k2,k3,h1,h2,h3);
	  if (bSym)
	    symfac = sym(l1,n1max)*sym(l2,n2max)*sym(l3,n3max);
	  else
	    symfac = 1.0;
	  
	  rsqal  = crsqal(acut,r1,k1,k2,k3,h1,h2,h3,nalias);

	  if (gdenom != 0)
	    qopt  += (symfac*(rsqal - sqr(gnumer)/gdenom));
	  if (debug)
	    fprintf(debug,"rsqal: %10.3e, gnumer: %10.3e, gdenom: %10.3e, ratio: %10.3e\n",
		    rsqal,gnumer,gdenom,gnumer/gdenom);
	
#ifdef DEBUG
	  if ((l1 == n1max/2) || (l2 == n2max/2) || (l3 == n3max/2))
	    printf("L(%2d,%2d,%2d)  D(%10.3e,%10.3e,%10.3e) U(%10.3e,%10.3e,%10.3e) gnumer=%10.3em dsq=%10.3e, gdenom=%10.3e, ghat=%10.3e\n",
		   l1,l2,l3,d1,d2,d3,u1,u2,u3,gnumer,dsq,gdenom,
		   (gdenom == 0.0) ? 0 : gnumer/gdenom);
#endif
	  if (!bSearch) {
	    if (gdenom != 0)
	      gg = gnumer/gdenom;
	    else
	      gg = 0.0;
	    ghat[l1][l2][l3] = gg/EPSILON0;
	    gsq = gg*gg;

	    ufrth = cufrth(porder,k1,k2,k3,h1,h2,h3,nalias);
	    pval  = pval + symfac*
	      (dsq*gsq*(usqsq(porder,k1,k2,k3,h1,h2,h3)-ufrth));
	    if (gdenom != 0)
	      zval  = zval + symfac*
		(dsq*gsq*ufrth - 2.0*gnumer*gnumer/gdenom + rsqal);
	    ss    = shat(acut,kmag,r1);
	    rsq   = ss*ss/ksq;
	    eref  = eref + symfac* (rsqal - rsq);
	  }
	}
      }
    }
  }
  if (bVerbose)
    fprintf(stderr,"\n");
  *ppval = pval/(box1*box2*box3);
  *zzval = zval/(box1*box2*box3);
  *eeref = eref/(box1*box2*box3);
  *qqopt = qopt/(EPSILON0*box1*box2*box3);
}

void calc1D(gmx_bool bSym,gmx_bool bVerbose,
	    const int n1max,const int n2max,const int n3max,
	    const real h1,const real h2,const real h3,
	    int nalias,int porder,real acut,real r1,const real alpha,
	    const gmx_bool bSearch,
	    real ***ghat,real *ppval,real *zzval,real *eeref,real *qqopt)
{     
  real box1,box2,box3;
  real k1,k2,k3;
  real gnumer,dsq,gdenom;
  real rsqal;
  real symfac;
  int  l1;
  real twopi=2*M_PI;
  real d1,u1;
  real pval,zval,eref,qopt;
  int  N1MAX;
/*   int  N2MAX,N3MAX; */
  
  if (bSym) {
    N1MAX = n1max/2+1;
/*     N2MAX = n2max/2+1; */
/*     N3MAX = n3max/2+1; */
  }
  else {
    N1MAX = n1max;
/*     N2MAX = n2max; */
/*     N3MAX = n3max; */
  }
    
  box1 = n1max*h1;
  box2 = n2max*h2;
  box3 = n3max*h3;

  pval = 0.0;
  zval = 0.0;
  eref = 0.0;
  qopt = 0.0;

  k2 = k3 = 0;
  
  for(l1=0; (l1<N1MAX); l1++) {
    if (bVerbose)
      fprintf(stderr,"\rl1=%5d  qopt=%12.6e",l1,qopt);
      
    k1   = twopi*l1/box1;
    d1   = dhat(alpha,k1,h1);
    
    if (l1 == 0) 
      ghat[0][0][0] = 0.0;
    else {
      u1   = ursum1D(XX,porder,acut,r1,k1,h1,nalias);
	  
      gnumer = d1*u1;
      dsq    = d1*d1;
      gdenom = dsq*usqsq(porder,k1,k2,k3,h1,h2,h3);
      if (bSym)
	symfac = sym(l1,n1max);
      else
	symfac = 1.0;
      
      rsqal  = crsqal(acut,r1,k1,k2,k3,h1,h2,h3,nalias);
      
      if (gdenom != 0)	  
	qopt  += symfac*(rsqal - (gnumer*gnumer)/gdenom);
    }
  }
  if (bVerbose)
    fprintf(stderr,"\n");
  *ppval = pval/(box1*box2*box3);
  *zzval = zval/(box1*box2*box3);
  *eeref = eref/(box1*box2*box3);
  *qqopt = qopt/(box1*box2*box3);
}

void read_params(char *fn,t_inputrec *ir,rvec boxs)
{
  real   t,lambda;
  int    step,natoms,m;
  matrix box;
  
  /* Read topology and coordinates */
  read_tpx(fn,&step,&t,&lambda,ir,box,&natoms,NULL,NULL,NULL,NULL);
  for(m=0; (m<DIM); m++)
    boxs[m] = box[m][m];
}
      
int main(int argc,char *argv[]) 
{
  /* FILE   *fp; */
  const  real gold=0.38197;
  const  int  mxiter=12;
  int    n1max,n2max,n3max;
  real   h1,h2,h3;
  real   box[3];
  int    nalias,porder;
  real   acut,alpha,r1;
  rvec   beta;
  gmx_bool   bSearch,bConv;
  /* gmx_bool   bNL=FALSE; */
  real   ***ghat;
  real   pval,zval,eref,qopt,norm;
  real   alpha0,alpha1,alpha2,alpha3,alptol;
  real   q0,q1,q2,q3;
  int    niter;
  /* int    ii,jj,kk,nn; */
  t_inputrec ir;
  t_filenm fnm[] = {
    { efTPX, NULL, NULL,   ffREAD },
    { efHAT, "-o", "ghat", ffWRITE }
  };
#define NFILE asize(fnm)
  static gmx_bool bVerbose=FALSE,bCubic=TRUE,bSym=TRUE;
  static t_pargs pa[] = {
    { "-v",     FALSE, etBOOL, &bVerbose, "Verbose on"},
    { "-cubic", FALSE, etBOOL, &bCubic,   "Force beta to be the same in all directions" },
    { "-sym",   FALSE, etBOOL, &bSym,     "HIDDENUse symmetry for the generation of ghat function (turn off for debugging only!)" }
  };
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,TRUE,NFILE,fnm,asize(pa),pa,0,NULL,0,NULL);
  
  read_params(ftp2fn(efTPX,NFILE,fnm),&ir,box);
  acut   = ir.rcoulomb;
  r1     = ir.rcoulomb_switch;
  n1max  = ir.nkx;
  n2max  = ir.nky;
  n3max  = ir.nkz;
  h1     = box[XX]/n1max;
  h2     = box[YY]/n2max;
  h3     = box[ZZ]/n3max;
  
  /* These are not parameters. They are fixed. 
   * Luty et al. determined their optimal values to be 2.
   */
  nalias = 2;
  porder = 2;
  /*bSym   = TRUE;*/
  
  set_LRconsts(stdout,r1,acut,box,NULL);  

  printf("Grid cell size is %8.4f x %8.4f x %8.4f nm\n",h1,h2,h3);
  
  ghat=mk_rgrid(n1max,n2max,n3max);
  
  bSearch = TRUE;

  niter = 0;
  if (!bSearch) {
    /* c          alpha = 1.33*/
    alpha = 1.0;
    calc(bSym,bVerbose,n1max,n2max,n3max,h1,h2,h3,
	 nalias,porder,acut,r1,alpha,bSearch,
	 ghat,&pval,&zval,&eref,&qopt);
  }
  else /* Elsje */ {
    alpha0 = 1.0;
    alpha1 = 4.0/3.0;
    alpha3 = 2.0;
    calc(bSym,bVerbose,n1max,n2max,n3max,h1,h2,h3,
	 nalias,porder,acut,r1,alpha0,bSearch,
	 ghat,&pval,&zval,&eref,&q0);
    calc(bSym,bVerbose,n1max,n2max,n3max,h1,h2,h3,
	 nalias,porder,acut,r1,alpha1,bSearch,
	 ghat,&pval,&zval,&eref,&q1);
    calc(bSym,bVerbose,n1max,n2max,n3max,h1,h2,h3,
	 nalias,porder,acut,r1,alpha3,bSearch,
	 ghat,&pval,&zval,&eref,&q3);
	 
    /* if ( (q1 > q0) || (q1 > q3) ) 
      gmx_fatal(FARGS,"oops q1=%f,q0=%f,q3=%f",q1,q0,q3); */
    alpha2 = alpha1 + gold*(alpha3-alpha1);
    calc(bSym,bVerbose,n1max,n2max,n3max,h1,h2,h3,
	 nalias,porder,acut,r1,alpha2,bSearch,
	 ghat,&pval,&zval,&eref,&q2);
    bConv  = FALSE;
    alptol = 0.05;
    
    fprintf(stderr,"q0 = %10g, q1= %10g, q2 = %10g, q3 = %10g\n",
	    q0,q1,q2,q3);
    
    while ((!bConv) && (niter < mxiter)) {
      fprintf(stderr,"%2d  %10.4f  %10.4f  %10.4f  %10.4f\n", 
	      niter,alpha0,alpha1,alpha2,alpha3);
      niter = niter + 1;
      if (q2 < q1) {
	alpha0 = alpha1;
	q0     = q1;
	alpha1 = alpha2;
	q1     = q2;
	alpha2 = alpha1 + gold*(alpha3-alpha1);
	calc(bSym,bVerbose,n1max,n2max,n3max,h1,h2,h3,
	     nalias,porder,acut,r1,alpha2,bSearch,
	     ghat,&pval,&zval,&eref,&q2);
      }
      else {
	alpha3 = alpha2;
	q3     = q2;
	alpha2 = alpha1;
	q2     = q1;
	alpha1 = alpha2 - gold*(alpha2-alpha0);
	calc(bSym,bVerbose,n1max,n2max,n3max,h1,h2,h3,
	     nalias,porder,acut,r1,alpha1,bSearch,
	     ghat,&pval,&zval,&eref,&q1);
      }
      if ((alpha3-alpha0) < alptol) 
	bConv = TRUE;
    }
    bSearch = FALSE;
    if (q1 < q2) {
      alpha = alpha1;
      calc(bSym,bVerbose,n1max,n2max,n3max,h1,h2,h3,
	   nalias,porder,acut,r1,alpha,bSearch,
	   ghat,&pval,&zval,&eref,&qopt);
    }
    else {
      alpha = alpha2;
      calc(bSym,bVerbose,n1max,n2max,n3max,h1,h2,h3,
	   nalias,porder,acut,r1,alpha,bSearch,
	   ghat,&pval,&zval,&eref,&qopt);
    }
    
    fprintf(stderr,"%10g  %10g  %10g  %10g  %10g  %10g\n",
	    acut,r1,pval,zval,eref,qopt);
    norm = sqr(1.0/(4.0*M_PI*h1*h1))*(4.0/3.0*M_PI*h1*h1*h1);
    if (norm != 0) {
      pval = sqrt(fabs(pval)/norm)*100.0;
      zval = sqrt(fabs(zval)/norm)*100.0;
      eref = sqrt(fabs(eref)/norm)*100.0;
      qopt = sqrt(fabs(qopt)/norm)*100.0;
    }
    beta[XX] = beta[YY] = beta[ZZ] = alpha;
    wr_ghat(ftp2fn(efHAT,NFILE,fnm),n1max,n2max,n3max,h1,h2,h3,
	    ghat,nalias,porder,niter,bSym,beta,
	    r1,acut,pval,zval,eref,qopt);
	    
    /*    fp=ftp2FILE(efHAT,NFILE,fnm,"w");
	  fprintf(fp,"%8d  %8d  %8d  %15.10e  %15.10e %15.10e\n",
	  n1max,n2max,n3max,h1,h2,h3);
	  fprintf(fp,"%8d  %8d  %8d  %8d  %15.10e  %15.10e  %15.10e\n",
	  nalias,porder,niter,bSym,alpha,alpha,alpha);
	  fprintf(fp,"%10g  %10g  %10g  %10g  %10g  %10g\n",
	  acut,r1,pval,zval,eref,qopt);
	  {
	  int  N1MAX,N2MAX,N3MAX;
	  
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
	  ffclose(fp);
	  }
    */
    pr_scalar_gk("ghat.xvg",n1max,n2max,n3max,box,ghat);
  }
  return 0;
}

