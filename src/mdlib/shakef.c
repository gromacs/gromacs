/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_shakef_c = "$Id$";

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
#include "callf77.h"

static void pv(FILE *log,char *s,rvec x)
{
  int m;

  fprintf(log,"%5s:",s);
  for(m=0; (m<DIM); m++)
    fprintf(log,"  %10.3f",x[m]);
  fprintf(log,"\n");
  fflush(log);
}

int shakef(FILE *log,
	   int natoms,real invmass[],int ncon,
	   t_iparams ip[],t_iatom *iatom,
	   real tol,rvec x[],rvec xp[],
	   real omega,bool bFEP,real lambda,real lagr[])
{
  /*
   *
   *
   *     r.c. van schaik and w.f. van gunsteren
   *     eth zuerich
   *     june 1992
   *     Adapted for use with Gromacs by David van der Spoel november 92.
   *
   *     This version of shake is adapted to conformational search in
   *     hyper-space (more than three dimensions). shake will supply
   *     corrections to given coordinates xp, such that for xp a list of
   *     constraints will be satisfied, each geometrically within a specified
   *     tolerance tol. the corrections are made along vectors derived from x.
   *     when xp = xp(t+dt) results from non-constrained dynamics time step
   *     dt, starting at x=x(t), the corrections are dynamically correct up
   *     to order dt**2. so the dynamical accuracy of shake depends on both
   *     tol and dt, whereas the geometrical accuracy only on tol.
   *
   *     Assume that either PBC has been taken into account already,
   *     or no PBC is used, ie.: no PBC calculations in shake
   *
   *     tol            = relative geometrical tolerance
   *     x[natoms][DIM]
   *                    = reference atom cartesian coordinates
   *     xp[natoms][DIM]
   *                    = atom cartesian coordinates that will be shaken
   *
   *     niter          = delivered with number of iterations
   *     DIM            = dimensionality of the problem
   *
   *     if 1000 or more iterations appear necessary, or if xp deviates too
   *     much from x (dt too large?), the subr. is returned with a message
   *     and niter = 0.
   *
   */
  
  static  bool *skip=NULL;
  static  atom_id *index=NULL;
  static  rvec *rij=NULL;
  static  int maxcon=0;
  
  
  int     nit,ll,i,j,k,k0,m,type,nindex;
  t_iatom *ia;
  rvec    xpij;
  real    *xpi,*xpj,*xij;
  real    L1,tol2,toler,diff,rpij2,rrpr,acor,xh;
  bool    ready;
  
  if (!skip) {
    snew(skip,2*natoms);
    snew(index,natoms);
  }
  if (ncon > maxcon) {
    srenew(rij,ncon);
    maxcon=ncon;
#ifdef DEBUG
    fprintf(log,"shake: maxcon = %d\n",maxcon);
#endif
  }

  L1    = 1.0-lambda;
  tol2  = 2.0*tol;
  nit   = 0;
  ready = FALSE;

  nindex=0;
  ia=iatom;
  for(ll=0; (ll<ncon); ll++,ia+=3) {
    for(m=1; (m<=2); m++) {
      i=ia[m];
      for(k=0; (k<nindex); k++)
	if (index[k]==i)
	  break;
      if (k==nindex) {
#ifdef DEBUG
	if (nindex >= natoms)
	  fatal_error(0,"Range check error in shake\n"
		      "nindex=%d, natoms=%d\n",nindex,natoms);
#endif
	index[nindex]=i;
	nindex++;
      }
    }
    i=ia[1];
    j=ia[2];
    rvec_sub(x[i],x[j],rij[ll]);
  }

  for(k=0; (k<nindex); k++) {
    k0=index[k];
#ifdef DEBUG
    if ((k0 < 0) || (k0 >= natoms))
      fatal_error(0,"Range check error in shake\n"
		  "k0=%d, natoms=%d\n",k0,natoms);
#endif
    skip[k0]        = TRUE;
    skip[natoms+k0] = FALSE;
  }
  
  while (! ready) {
    /*       too many iterations? */
    if (nit > 1000) {
      fprintf(log,
	      " coordinate resetting (shake) was not"
	      " accomplished within 1000 iterations\n\n");
      return 0;
    }
    
    ready=TRUE;
    
    /* loop over all the constraints */
    ia=iatom;
    for(ll=0; (ll<ncon); ll++,ia+=3) {
      i    = ia[1];
      j    = ia[2];

      if ((skip[natoms+i]) && (skip[natoms+j]))
	continue;

#ifdef DEBUGSHAKE
      fprintf(log,"ll: %d\n",ll);
      pv(log,"xij    ",xij);
      pv(log,"rij[ll]",rij[ll]);
#endif
      xij   = rij[ll];
#ifdef DEBUGSHAKE
      fprintf(log,"ll: %d\n",ll);
      pv(log,"xij    ",xij);
      pv(log,"rij[ll]",rij[ll]);
#endif
      type  = ia[0];
      xpi   = xp[i];
      xpj   = xp[j];
      if (bFEP) 
	toler = sqr(L1*ip[type].shake.dA + lambda*ip[type].shake.dB);
      else
	toler = sqr(ip[type].shake.dA);
      diff  = toler;
      rpij2 = 0.0;
      
      rvec_sub(xpi,xpj,xpij);
      rpij2=iprod(xpij,xpij);
      
      diff -= rpij2;
      if (fabs(diff) < toler*tol2)
	continue;
      
      rrpr=iprod(xij,xpij);
      
      if (rrpr < toler*1.e-6) {
	fprintf(log," oops, coordinate resetting cannot be"
		" accomplished, deviation is too large\n"
		" nit = %5d  ll    = %5d\n"
		" i   = %5d  j     = %5d\n"
		" rrpr= %10f\n",
		nit,ll,i+1,j+1,rrpr);
	pv(log,"xij",xij);
	pv(log,"xpij",xpij);
	return 0;
      }
      
      acor=(omega*diff)/(rrpr*(invmass[i]+invmass[j])*2.0);
      lagr[ll] += acor;
      for(m=0; (m<DIM); m++) {
	xh      = xij[m]*acor;
	xpi[m] += xh*invmass[i];
	xpj[m] -= xh*invmass[j];
      }
      
      skip[i] = FALSE;
      skip[j] = FALSE;
      ready   = FALSE;
    }
    
    nit=nit+1;
    for(k=0; (k<nindex); k++) {
      k0=index[k];
      skip[natoms+k0] = skip[k0];
      skip[k0]        = TRUE;
    }
  }
  
  return nit;
}

int vec_shakef(FILE *log,
	       int natoms,real invmass[],int ncon,
	       t_iparams ip[],t_iatom *iatom,
	       real tol,rvec x[],rvec xp[],
	       bool bFEP,real lambda,real lagr[])
{
  static  rvec *rij=NULL;
  static  real *M2=NULL,*tt=NULL,*dist2=NULL;
  static  int  maxcon=0;
  int     maxnit=1000;
  int     nit,ll,i,j,type;
  t_iatom *ia;
  real    L1,tol2,toler;
  real    mm;
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

  /* We have a FORTRAN shake now! */  
#ifdef USE_FORTRAN
  F77_FUNC(fshake,FSHAKE)(iatom,&ncon,&nit,&maxnit,dist2,xp[0],
			  rij[0],M2,invmass,tt,lagr,&error);
#else
  /* And a c shake also ! */
  cshake(iatom,ncon,&nit,maxnit,dist2,xp[0],rij[0],M2,invmass,tt,lagr,&error);
#endif
  if (nit >= maxnit) {
    fprintf(log,"Shake did not converge in %d steps\n",maxnit);
    nit=0;
  }
  else if (error != 0) {
    fprintf(log,"Inproduct between old and new vector = 0.0!\n"
	    "constraint #%d atoms %u and %u\n",
	    error-1,iatom[3*(error-1)+1],iatom[3*(error-1)+2]);
    nit=0;
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

int bshakef(FILE *log,int natoms,real invmass[],int nblocks,int sblock[],
	    t_idef *idef,t_inputrec *ir,matrix box,rvec x_s[],rvec xp[],
	    t_nrnb *nrnb,real lambda,real *dvdlambda)
{
  static  bool bFirst=TRUE;
  static  bool bSafe;

  /* SOR Stuff from Barth et al. JCC 16 (1996) 1192-1209 */
  static  bool bSOR=FALSE;
  static  real delta=0.1;
  static  real omega=1.0;
  static  int  gamma=1000000;
  static  real *lagr;
  
  t_iatom *iatoms;
  real    *lam,dt_2,dvdl;
  int     i,n0,ncons,blen,type;
  int     tnit=0,trij=0;
  
#ifdef DEBUG
  fprintf(log,"nblocks=%d, sblock[0]=%d\n",nblocks,sblock[0]);
#endif
  ncons=idef->il[F_SHAKE].nr/3;
  if (bFirst) {
    bSafe=(getenv("NOVECSHAKE") != NULL);
    bSOR=(getenv("NOSOR") == NULL);
    if (bSOR) 
      please_cite(log,"Barth95a");
    snew(lagr,ncons);
    bFirst=FALSE;
  }
  for(i=0; i<ncons; i++)
    lagr[i] =0;
  
  iatoms = &(idef->il[F_SHAKE].iatoms[sblock[0]]);
  lam    = lagr;
  for(i=0; (i<nblocks); ) {
    blen=(sblock[i+1]-sblock[i]);
    blen/=3;
    
    if (bSafe)
      n0=shakef(log,natoms,invmass,blen,idef->iparams,
		iatoms,ir->shake_tol,x_s,xp,omega,
		ir->efep!=efepNO,lambda,lam);
    else
      n0=vec_shakef(log,natoms,invmass,blen,idef->iparams,
		    iatoms,ir->shake_tol,x_s,xp,
		    ir->efep!=efepNO,lambda,lam);
    
#ifdef DEBUGSHAKE
	check_cons(log,blen,x_s,xp,idef->iparams,iatoms,invmass);
#endif
    
    if (n0==0) {
      check_cons(log,blen,x_s,xp,idef->iparams,iatoms,invmass);
      return -1;
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
  if (bSOR) {
    if (tnit > gamma) {
      delta = -0.5*delta;
    }
    omega = omega + delta;
    gamma = tnit;
  }
  inc_nrnb(nrnb,eNR_SHAKE,tnit);
  inc_nrnb(nrnb,eNR_SHAKE_RIJ,trij);
  
  return 0;
}

