/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "pbc.h"
#include "txtdump.h"
#include "vec.h"
#include "nrnb.h"
#include "constr.h"
#include "copyrite.h"


static void pv(FILE *log,char *s,rvec x)
{
  int m;

  fprintf(log,"%5s:",s);
  for(m=0; (m<DIM); m++)
    fprintf(log,"  %10.3f",x[m]);
  fprintf(log,"\n");
  fflush(log);
}

void cshake(atom_id iatom[],int ncon,int *nnit,int maxnit,
	    real dist2[],real xp[],real rij[],real m2[],real omega,
	    real invmass[],real tt[],real lagr[],int *nerror)
{
  /*
   *     r.c. van schaik and w.f. van gunsteren
   *     eth zuerich
   *     june 1992
   *     Adapted for use with Gromacs by David van der Spoel november 92 and later.
   */
  const   real mytol=1e-6;
  
  int     ll,i,j,i3,j3,l3;
  int     ix,iy,iz,jx,jy,jz;
  real    toler,rpij2,rrpr,tx,ty,tz,diff,acor,im,jm;
  real    xh,yh,zh,rijx,rijy,rijz;
  real    tix,tiy,tiz;
  real    tjx,tjy,tjz;
  int     nit,error,iconv,nconv;
  
  error=0;
  nconv=1;
  for (nit=0; (nit<maxnit) && (nconv != 0) && (error == 0); nit++) {
    nconv=0;
    for(ll=0; (ll<ncon) && (error == 0); ll++) {
      l3    = 3*ll;
      rijx  = rij[l3+XX];
      rijy  = rij[l3+YY];
      rijz  = rij[l3+ZZ];
      i     = iatom[l3+1];
      j     = iatom[l3+2];
      i3    = 3*i;
      j3    = 3*j;
      ix    = i3+XX;
      iy    = i3+YY;
      iz    = i3+ZZ;
      jx    = j3+XX;
      jy    = j3+YY;
      jz    = j3+ZZ;
      
      tx      = xp[ix]-xp[jx];
      ty      = xp[iy]-xp[jy];
      tz      = xp[iz]-xp[jz];
      rpij2   = tx*tx+ty*ty+tz*tz;
      toler   = dist2[ll];
      diff    = toler-rpij2;
      
      /* iconv is zero when the error is smaller than a bound */
      iconv   = fabs(diff)*tt[ll];
      
      if (iconv != 0) {
	nconv   = nconv + iconv;
	rrpr    = rijx*tx+rijy*ty+rijz*tz;
	
	if (rrpr < toler*mytol) 
	  error=ll;
	else {
	  acor      = omega*diff*m2[ll]/rrpr;
	  lagr[ll] += acor;
	  xh        = rijx*acor;
	  yh        = rijy*acor;
	  zh        = rijz*acor;
	  im        = invmass[i];
	  jm        = invmass[j];
	  xp[ix] += xh*im;
	  xp[iy] += yh*im;
	  xp[iz] += zh*im;
	  xp[jx] -= xh*jm;
	  xp[jy] -= yh*jm;
	  xp[jz] -= zh*jm;
	}
      }
    }
  }
  *nnit=nit;
  *nerror=error;
}

int vec_shakef(FILE *log,
	       int natoms,real invmass[],int ncon,
	       t_iparams ip[],t_iatom *iatom,
	       real tol,rvec x[],rvec xp[],real omega,
	       bool bFEP,real lambda,real lagr[],
	       bool bCalcVir,tensor rmdr)
{
  static  rvec *rij=NULL;
  static  real *M2=NULL,*tt=NULL,*dist2=NULL;
  static  int  maxcon=0;
  int     maxnit=1000;
  int     nit,ll,i,j,type;
  t_iatom *ia;
  real    L1,tol2,toler;
  real    mm,tmp;
  int     error;
    
  if (ncon > maxcon) {
    srenew(rij,ncon);
    srenew(M2,ncon);
    srenew(tt,ncon);
    srenew(dist2,ncon);
    maxcon=ncon;
#ifdef DEBUG
    fprintf(log,"shake: maxcon = %d\n",maxcon);
#endif
  }

  L1=1.0-lambda;
  tol2=2.0*tol;
  ia=iatom;
  for(ll=0; (ll<ncon); ll++,ia+=3) {
    type  = ia[0];
    i=ia[1];
    j=ia[2];
    
    mm=2*(invmass[i]+invmass[j]);
    rij[ll][XX]=x[i][XX]-x[j][XX];
    rij[ll][YY]=x[i][YY]-x[j][YY];
    rij[ll][ZZ]=x[i][ZZ]-x[j][ZZ];
    M2[ll]=1.0/mm;
    if (bFEP) 
      toler = sqr(L1*ip[type].shake.dA + lambda*ip[type].shake.dB);
    else
      toler = sqr(ip[type].shake.dA);
    dist2[ll] = toler;
    tt[ll] = 1.0/(toler*tol2);
  }

  cshake(iatom,ncon,&nit,maxnit,dist2,xp[0],rij[0],M2,omega,invmass,tt,lagr,&error);

  if (nit >= maxnit) {
    fprintf(log,"Shake did not converge in %d steps\n",maxnit);
    nit=0;
  }
  else if (error != 0) {
    fprintf(log,"Inner product between old and new vector <= 0.0!\n"
	    "constraint #%d atoms %u and %u\n",
	    error-1,iatom[3*(error-1)+1]+1,iatom[3*(error-1)+2]+1);
    nit=0;
  }

  /* Constraint virial and correct the lagrange multipliers for the length */
  ia=iatom;
  for(ll=0; (ll<ncon); ll++,ia+=3) {
    if (bCalcVir) {
      mm = lagr[ll];
      for(i=0; i<DIM; i++) {
	tmp = mm*rij[ll][i];
	for(j=0; j<DIM; j++)
	  rmdr[i][j] -= tmp*rij[ll][j];
      }
      /* 21 flops */
    }

    type  = ia[0];
    if (bFEP) 
      toler = L1*ip[type].shake.dA + lambda*ip[type].shake.dB;
    else
      toler = ip[type].shake.dA;
    lagr[ll] *= toler;
  }
  
  return nit;
}

static void check_cons(FILE *log,int nc,rvec x[],rvec xp[],
		       t_iparams ip[],t_iatom *iatom,
		       real invmass[])
{
  t_iatom *ia;
  int     ai,aj;
  int     i;
  real    d,dp;
  rvec    dx;

  fprintf(log,
	  "    i     mi      j     mj      before       after   should be\n");
  ia=iatom;
  for(i=0; (i<nc); i++,ia+=3) {
    ai=ia[1];
    aj=ia[2];
    rvec_sub(x[ai],x[aj],dx);
    d=norm(dx);
    rvec_sub(xp[ai],xp[aj],dx);
    dp=norm(dx);
    fprintf(log,"%5d  %5.2f  %5d  %5.2f  %10.5f  %10.5f  %10.5f\n",
	    ai+1,1.0/invmass[ai],
	    aj+1,1.0/invmass[aj],d,dp,ip[ia[0]].shake.dA);
  }
}

bool bshakef(FILE *log,int natoms,real invmass[],int nblocks,int sblock[],
	     t_idef *idef,t_inputrec *ir,matrix box,rvec x_s[],rvec xp[],
	     t_nrnb *nrnb,real lambda,real *dvdlambda,
	     bool bCalcVir,tensor rmdr,bool bDumpOnError)
{
  static  bool bFirst=TRUE;
  static  real *lagr;
  /* Stuff for successive overrelaxation */
  static  real delta=0.1;
  static  real omega=1.0;
  static  int  gamma=1000000;
  
  t_iatom *iatoms;
  real    *lam,dt_2,dvdl;
  int     i,n0,ncons,blen,type;
  int     tnit=0,trij=0;
  
#ifdef DEBUG
  fprintf(log,"nblocks=%d, sblock[0]=%d\n",nblocks,sblock[0]);
#endif
  ncons=idef->il[F_SHAKE].nr/3;
  if (bFirst) {
    if (ir->bShakeSOR) 
      please_cite(log,"Barth95a");
    snew(lagr,ncons);
    bFirst=FALSE;
  }
  for(i=0; i<ncons; i++)
    lagr[i] =0;
  
  iatoms = &(idef->il[F_SHAKE].iatoms[sblock[0]]);
  lam    = lagr;
  for(i=0; (i<nblocks); ) {
    blen  = (sblock[i+1]-sblock[i]);
    blen /= 3;
    n0 = vec_shakef(log,natoms,invmass,blen,idef->iparams,
		    iatoms,ir->shake_tol,x_s,xp,omega,
		    ir->efep!=efepNO,lambda,lam,bCalcVir,rmdr);
#ifdef DEBUGSHAKE
    check_cons(log,blen,x_s,xp,idef->iparams,iatoms,invmass);
#endif
    
    if (n0 == 0) {
      if (bDumpOnError)
	check_cons(log,blen,x_s,xp,idef->iparams,iatoms,invmass);
      return FALSE;
    }
    tnit   += n0*blen;
    trij   += blen;
    iatoms += 3*blen;	/* Increment pointer! */
    lam    += blen;
    i++;
  }
  if (ir->efep != efepNO) {
    dt_2 = 1/sqr(ir->delta_t);
    dvdl = 0;
    for(i=0; i<ncons; i++) {
      type = idef->il[F_SHAKE].iatoms[3*i];
      dvdl += lagr[i]*dt_2*
	(idef->iparams[type].shake.dB-idef->iparams[type].shake.dA);
    }
    *dvdlambda += dvdl;
  }
#ifdef DEBUG
  fprintf(log,"tnit: %5d  omega: %10.5f\n",tnit,omega);
#endif
  if (ir->bShakeSOR) {
    if (tnit > gamma) {
      delta = -0.5*delta;
    }
    omega = omega + delta;
    gamma = tnit;
  }
  inc_nrnb(nrnb,eNR_SHAKE,tnit);
  inc_nrnb(nrnb,eNR_SHAKE_RIJ,trij);
  if (bCalcVir)
    inc_nrnb(nrnb,eNR_CONSTR_VIR,trij);
  
  return TRUE;
}

