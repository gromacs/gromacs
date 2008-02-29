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
#include "smalloc.h"
#include "mdrun.h"
#include "nrnb.h"

static void do_lincsp(rvec *x,rvec *f,rvec *fp,t_pbc *pbc,
		      t_lincsdata *lincsd,real *invmass,
		      int nrec)
{
  int     b,i,j,k,n,it,rec;
  real    tmp0,tmp1,tmp2,im1,im2,mvb,rlen,len,wfac,lam;  
  real    u0,u1,u2,v0,v1,v2;
  rvec    dx;
  int     ncons,*bla,*blnr,*blbnb;
  rvec    *r;
  real    *blc,*blcc,*blm,*rhs1,*rhs2,*sol,*swap;

  ncons  = lincsd->nc;
  bla    = lincsd->bla;
  r      = lincsd->tmpv;
  blnr   = lincsd->blnr;
  blbnb  = lincsd->blbnb;
  blc    = lincsd->blc;
  blcc   = lincsd->blcc;
  blm    = lincsd->tmpncc;
  rhs1   = lincsd->tmp1;
  rhs2   = lincsd->tmp2;
  sol    = lincsd->tmp3;

  /* Compute normalized i-j vectors */
  if (pbc) {
    for(b=0; b<ncons; b++) {
      pbc_dx(pbc,x[bla[2*b]],x[bla[2*b+1]],dx);
      unitv(dx,r[b]);
    }
  } else {
    for(b=0; b<ncons; b++) {
      rvec_sub(x[bla[2*b]],x[bla[2*b+1]],dx);
      unitv(dx,r[b]);
    } /* 16 ncons flops */
  }
  
  for(b=0; b<ncons; b++) {
    tmp0 = r[b][0];
    tmp1 = r[b][1];
    tmp2 = r[b][2];
    i = bla[2*b];
    j = bla[2*b+1];
    for(n=blnr[b]; n<blnr[b+1]; n++) {
      k = blbnb[n];
      blm[n] = blcc[n]*(tmp0*r[k][0] + tmp1*r[k][1] + tmp2*r[k][2]); 
    } /* 6 nr flops */
    mvb = blc[b]*(tmp0*(f[i][0] - f[j][0]) +
		  tmp1*(f[i][1] - f[j][1]) +    
		  tmp2*(f[i][2] - f[j][2]));
    rhs1[b] = mvb;
    sol[b]  = mvb;
    /* 7 flops */
  }
  /* Together: 23*ncons + 6*nrtot flops */
    
  for(rec=0; rec<nrec; rec++) {
    for(b=0; b<ncons; b++) {
      mvb = 0;
      for(n=blnr[b]; n<blnr[b+1]; n++) {
	j = blbnb[n];
	mvb = mvb + blm[n]*rhs1[j];
      }
      rhs2[b] = mvb;
      sol[b] = sol[b] + mvb;
    }
    swap = rhs1;
    rhs1 = rhs2;
    rhs2 = swap;
  } /* nrec*(ncons+2*nrtot) flops */
  
  for(b=0;b<ncons;b++) {
    i = bla[2*b];
    j = bla[2*b+1];
    mvb = blc[b]*sol[b];
    im1 = invmass[i];
    im2 = invmass[j];
    tmp0 = r[b][0]*mvb;
    tmp1 = r[b][1]*mvb;
    tmp2 = r[b][2]*mvb;
    u0 = fp[i][0] - tmp0*im1;
    u1 = fp[i][1] - tmp1*im1;
    u2 = fp[i][2] - tmp2*im1;
    v0 = fp[j][0] + tmp0*im2;
    v1 = fp[j][1] + tmp1*im2;
    v2 = fp[j][2] + tmp2*im2;
    fp[i][0] = u0;
    fp[i][1] = u1;
    fp[i][2] = u2;
    fp[j][0] = v0;
    fp[j][1] = v1;
    fp[j][2] = v2;
  } /* 16 ncons flops */
}

static void do_lincs(rvec *x,rvec *xp,t_pbc *pbc,
		     t_lincsdata *lincsd,real *invmass,
		     int nit,int nrec,
		     real wangle,int *warn,
		     bool bCalcVir,tensor rmdr)
{
  int     b,i,j,k,n,it,rec;
  real    tmp0,tmp1,tmp2,im1,im2,mvb,rlen,len,wfac,lam;  
  real    u0,u1,u2,v0,v1,v2;
  rvec    dx;
  int     ncons,*bla,*blnr,*blbnb;
  rvec    *r;
  real    *blc,*blcc,*bllen,*blm,*rhs1,*rhs2,*sol,*lambda,*swap;

  ncons  = lincsd->nc;
  bla    = lincsd->bla;
  r      = lincsd->tmpv;
  blnr   = lincsd->blnr;
  blbnb  = lincsd->blbnb;
  blc    = lincsd->blc;
  blcc   = lincsd->blcc;
  bllen  = lincsd->bllen;
  blm    = lincsd->tmpncc;
  rhs1   = lincsd->tmp1;
  rhs2   = lincsd->tmp2;
  sol    = lincsd->tmp3;
  lambda = lincsd->lambda;

  *warn = 0;

  if (pbc) {
    /* Compute normalized i-j vectors */
    for(b=0; b<ncons; b++) {
      pbc_dx(pbc,x[bla[2*b]],x[bla[2*b+1]],dx);
      unitv(dx,r[b]);
    }  
    for(b=0; b<ncons; b++) {
      for(n=blnr[b]; n<blnr[b+1]; n++) {
	blm[n] = blcc[n]*iprod(r[b],r[blbnb[n]]);
      }
      pbc_dx(pbc,xp[bla[2*b]],xp[bla[2*b+1]],dx);
      mvb = blc[b]*(iprod(r[b],dx) - bllen[b]);
      rhs1[b] = mvb;
      sol[b]  = mvb;
    }
  } else {
    /* Compute normalized i-j vectors */
    for(b=0;b<ncons;b++) {
      i = bla[2*b];
      j = bla[2*b+1];
      tmp0 = x[i][0] - x[j][0];
      tmp1 = x[i][1] - x[j][1];
      tmp2 = x[i][2] - x[j][2];
      rlen = invsqrt(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2);
      r[b][0] = rlen*tmp0;
      r[b][1] = rlen*tmp1;
      r[b][2] = rlen*tmp2;
    } /* 16 ncons flops */
    
    for(b=0;b<ncons;b++) {
      tmp0 = r[b][0];
      tmp1 = r[b][1];
      tmp2 = r[b][2];
      len = bllen[b];
      i = bla[2*b];
      j = bla[2*b+1];
      for(n=blnr[b];n<blnr[b+1];n++) {
	k = blbnb[n];
	blm[n] = blcc[n]*(tmp0*r[k][0] + tmp1*r[k][1] + tmp2*r[k][2]); 
      } /* 6 nr flops */
      mvb = blc[b]*(tmp0*(xp[i][0] - xp[j][0]) +
		    tmp1*(xp[i][1] - xp[j][1]) +    
		    tmp2*(xp[i][2] - xp[j][2]) - len);
      rhs1[b]=mvb;
      sol[b]=mvb;
      /* 8 flops */
    }
    /* Together: 24*ncons + 6*nrtot flops */
  }
    
  for(rec=0; rec<nrec; rec++) {
    for(b=0; b<ncons; b++) {
      mvb = 0;
      for(n=blnr[b]; n<blnr[b+1]; n++) {
	j = blbnb[n];
	mvb = mvb + blm[n]*rhs1[j];
      }
      rhs2[b] = mvb;
      sol[b] = sol[b] + mvb;
    }
    swap = rhs1;
    rhs1 = rhs2;
    rhs2 = swap;
  } /* nrec*(ncons+2*nrtot) flops */
  
  for(b=0; b<ncons; b++) {
    i = bla[2*b];
    j = bla[2*b+1];
    mvb = blc[b]*sol[b];
    lambda[b] = -mvb;
    im1 = invmass[i];
    im2 = invmass[j];
    tmp0 = r[b][0]*mvb;
    tmp1 = r[b][1]*mvb;
    tmp2 = r[b][2]*mvb;
    u0 = xp[i][0] - tmp0*im1;
    u1 = xp[i][1] - tmp1*im1;
    u2 = xp[i][2] - tmp2*im1;
    v0 = xp[j][0] + tmp0*im2;
    v1 = xp[j][1] + tmp1*im2;
    v2 = xp[j][2] + tmp2*im2;
    xp[i][0] = u0;
    xp[i][1] = u1;
    xp[i][2] = u2;
    xp[j][0] = v0;
    xp[j][1] = v1;
    xp[j][2] = v2;
  } /* 16 ncons flops */
  
  
  
  /*     
  ********  Correction for centripetal effects  ********  
  */
  
  wfac = cos(DEG2RAD*wangle);
  wfac = wfac*wfac;

  for(it=0; it<nit; it++) {
  
    for(b=0;b<ncons;b++) {
      len = bllen[b];
      if (pbc) {
	pbc_dx(pbc,xp[bla[2*b]],xp[bla[2*b+1]],dx);
      } else {
	rvec_sub(xp[bla[2*b]],xp[bla[2*b+1]],dx);
      }
      u1 = len*len;
      u0 = 2*u1 - norm2(dx);
      if (u0 < wfac*u1) *warn = b;	
      if (u0 < 0) u0 = 0;
      mvb = blc[b]*(len - u0*invsqrt(u0));
      rhs1[b] = mvb;
      sol[b]  = mvb;
    } /* 18*ncons flops */
    
    for(rec=0; rec<nrec; rec++) {
      for(b=0; b<ncons; b++) {
	mvb = 0;
	for(n=blnr[b]; n<blnr[b+1]; n++) {
	  j = blbnb[n];
	  mvb = mvb + blm[n]*rhs1[j];
	}
	rhs2[b] = mvb;
	sol[b] = sol[b] + mvb;
      }
      swap = rhs1;
      rhs1 = rhs2;
      rhs2 = swap;
    } /* nrec*(ncons+2*nrtot) flops */ 
    
    for(b=0; b<ncons; b++) {
      i = bla[2*b];
      j = bla[2*b+1];
      lam = lambda[b];
      mvb = blc[b]*sol[b];
      lambda[b] = lam - mvb;
      im1 = invmass[i];
      im2 = invmass[j];
      tmp0 = r[b][0]*mvb;
      tmp1 = r[b][1]*mvb;
      tmp2 = r[b][2]*mvb;
      u0 = xp[i][0] - tmp0*im1;
      u1 = xp[i][1] - tmp1*im1;
      u2 = xp[i][2] - tmp2*im1;
      v0 = xp[j][0] + tmp0*im2;
      v1 = xp[j][1] + tmp1*im2;
      v2 = xp[j][2] + tmp2*im2;
      xp[i][0] = u0;
      xp[i][1] = u1;
      xp[i][2] = u2;
      xp[j][0] = v0;
      xp[j][1] = v1;
      xp[j][2] = v2;
    } /* 17 ncons flops */
  } /* nit*ncons*(35+9*nrec) flops */

  if (bCalcVir) {
    /* Constraint virial */
    for(b=0; b<ncons; b++) {
      tmp0 = bllen[b]*lambda[b];
      for(i=0; i<DIM; i++) {
	tmp1 = tmp0*r[b][i];
	for(j=0; j<DIM; j++) {
	  rmdr[i][j] -= tmp1*r[b][j];
	}
      }
    } /* 22 ncons flops */
  }

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

void set_lincs_matrix(t_lincsdata *li,real *invmass)
{
  int i,a1,a2,n,k,sign,center;

  for(i=0; (i<li->nc); i++) {
    a1 = li->bla[2*i];
    a2 = li->bla[2*i+1];
    li->blc[i] = invsqrt(invmass[a1] + invmass[a2]);
  }

  /* Construct the coupling coefficient matrix blcc */
  for(i=0; (i<li->nc); i++) {
    a1 = li->bla[2*i];
    a2 = li->bla[2*i+1];
    for(n=li->blnr[i]; (n<li->blnr[i+1]); n++) {
      k = li->blbnb[n];
      if (a1 == li->bla[2*k] || a2 == li->bla[2*k+1])
	sign = -1;
      else
	sign = 1;
      if (a1 == li->bla[2*k] || a1 == li->bla[2*k+1])
	center = a1;
      else
	center = a2;
      li->blcc[n] = sign*invmass[center]*li->blc[i]*li->blc[k];
    }
  }
}

static int int_comp(const void *a,const void *b)
{
  return (*(int *)a) - (*(int *)b);
}

t_lincsdata *init_lincs(FILE *log,t_idef *idef,int start,int homenr,
			bool bDynamics)
{
  t_iatom     *iatom;
  int         i,j,k;
  int         type,a1,a2;
  real        lenA=0,lenB;
  int         **at_c,*at_cn;
  t_lincsdata *li;
  
  snew(li,1);
  li->nc = 0;
  li->nzerolen = 0;
  li->ncc = 0;
  
  if (idef->il[F_SHAKE].nr > 0) {
    /* Make atom-constraint connection list for temporary use */
    snew(at_c,homenr);
    snew(at_cn,homenr);

    iatom=idef->il[F_SHAKE].iatoms;

    /* Make a list of constraint connections */
    for(i=0; i<idef->il[F_SHAKE].nr; i+=3) {
      type = iatom[i];
      a1   = iatom[i+1];
      a2   = iatom[i+2];
      lenA = idef->iparams[type].shake.dA;
      lenB = idef->iparams[type].shake.dB;
      /* Skip the flexible constraints when not doing dynamics */
      if (bDynamics || lenA!=0 || lenB!=0) {
	if (at_cn[a1-start] % 4 == 0)
	  srenew(at_c[a1-start],at_cn[a1-start]+4);
	at_c[a1-start][at_cn[a1-start]] = li->nc;
	at_cn[a1-start]++;
	if (at_cn[a2-start] % 4 == 0)
	  srenew(at_c[a2-start],at_cn[a2-start]+4);
	at_c[a2-start][at_cn[a2-start]] = li->nc;
	at_cn[a2-start]++;
	li->nc++;
      }
    }
    
    snew(li->bllen0,li->nc);
    snew(li->ddist,li->nc);
    snew(li->bla,2*li->nc);

    j = 0;
    for(i=0; i<idef->il[F_SHAKE].nr; i+=3) {
      type = iatom[i];
      a1   = iatom[i+1];
      a2   = iatom[i+2];
      lenA = idef->iparams[type].shake.dA;
      lenB = idef->iparams[type].shake.dB;
      if (lenA == 0 && lenB == 0)
	li->nzerolen++;
      /* Skip the flexible constraints when not doing dynamics */
      if (bDynamics || lenA!=0 || lenB!=0) {
	li->bllen0[j] = lenA;
	li->ddist[j]  = lenB - lenA;
	li->bla[2*j]   = a1;
	li->bla[2*j+1] = a2;
	/* Count the constraint connections */
	li->ncc += at_cn[a1-start] + at_cn[a2-start] - 2;
	j++;
      }
    }      
    
    snew(li->blc,li->nc);
    snew(li->blnr,li->nc+1);
    snew(li->blbnb,li->ncc);
    snew(li->blcc,li->ncc);
    snew(li->bllen,li->nc);
    snew(li->tmpv,li->nc);
    snew(li->tmpncc,li->ncc);
    snew(li->tmp1,li->nc);
    snew(li->tmp2,li->nc);
    snew(li->tmp3,li->nc);
    snew(li->lambda,li->nc);
   
    /* Make constraint-neighbor list */
    li->blnr[0] = 0;
    for(j=0; j<li->nc; j++) {
      a1 = li->bla[2*j];
      a2 = li->bla[2*j+1];
      /* Set the length to the topology A length */
      li->bllen[j] = li->bllen0[j];
      /* Construct the constraint connection matrix blbnb */
      li->blnr[j+1] = li->blnr[j];
      for(k=0; k<at_cn[a1-start]; k++)
	if (at_c[a1-start][k] != j)
	  li->blbnb[(li->blnr[j+1])++] = at_c[a1-start][k];
      for(k=0; k<at_cn[a2-start]; k++)
	if (at_c[a2-start][k] != j)
	  li->blbnb[(li->blnr[j+1])++] = at_c[a2-start][k];
      /* Order the blbnb matrix to optimize memory access */
      qsort(&(li->blbnb[li->blnr[j]]),li->blnr[j+1]-li->blnr[j],
            sizeof(li->blbnb[0]),int_comp);
    }

    sfree(at_cn);
    for(i=0; i<homenr; i++)
      sfree(at_c[i]);
    sfree(at_c);
    
    fprintf(log,"\nInitializing LINear Constraint Solver\n");
    fprintf(log,"  number of constraints is %d\n",li->nc);
    if (li->nc > 0)
      fprintf(log,"  average number of constraints coupled to one constraint is %.1f\n",
	    (real)(li->ncc)/li->nc);
    if (li->nzerolen)
      fprintf(log,"  found %d constraints with zero length\n",li->nzerolen);
    fprintf(log,"\n");
    fflush(log);
  }

  return li;
}

/* Yes. I _am_ awfully ashamed of introducing a static variable after preaching
 * that we should get rid of them. Do as I say, not as I do!
 *
 * In defense, this is a hack just before the 3.3 release. The proper solution
 * is to introduce an abstract constraint object where it is stored. /EL 20050830
 */
static int lincs_warncount = 0;

#define LINCS_MAXWARN 10000

static void lincs_warning(rvec *x,rvec *xprime,t_pbc *pbc,
			  int ncons,int *bla,real *bllen,real wangle)
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
    i = bla[2*b];
    j = bla[2*b+1];
    if (pbc) {
      pbc_dx(pbc,x[i],x[j],v0);
      pbc_dx(pbc,xprime[i],xprime[j],v1);
    } else {
      rvec_sub(x[i],x[j],v0);
      rvec_sub(xprime[i],xprime[j],v1);
    }
    d0 = norm(v0);
    d1 = norm(v1);
    cosine = iprod(v0,v1)/(d0*d1);
    if (cosine<wfac) {
      sprintf(buf," %6d %6d  %5.1f  %8.4f %8.4f    %8.4f\n",
	      i+1,j+1,RAD2DEG*acos(cosine),d0,d1,bllen[b]);
      fprintf(stderr,buf);
      fprintf(stdlog,buf);
      lincs_warncount++;
    }
  }
  if(lincs_warncount > LINCS_MAXWARN)
  {
      gmx_fatal(FARGS,
                "Too many LINCS warnings (%d) - aborting to avoid logfile runaway.\n"
                "This normally happens when your system is not sufficiently equilibrated,"
                "or if you are changing lambda too fast in free energy simulations.\n"
                "If you know what you are doing you can adjust the lincs warning threshold\n"
                "in your mdp file, but normally it is better to fix the problem.\n",
                lincs_warncount);
  }
}

static void cconerr(real *max,real *rms,int *imax,rvec *x,t_pbc *pbc,
		    int ncons,int *bla,real *bllen)
{
  real      len,d,ma,ms,r2;
  int       b,im;
  rvec      dx;
  
  ma = 0;
  ms = 0;
  im = 0;
  for(b=0;b<ncons;b++) {
    if (pbc) {
      pbc_dx(pbc,x[bla[2*b]],x[bla[2*b+1]],dx);
    } else {
      rvec_sub(x[bla[2*b]],x[bla[2*b+1]],dx);
    }
    r2 = norm2(dx);
    len = r2*invsqrt(r2);
    d = fabs(len/bllen[b]-1);
    if (d > ma) {
      ma = d;
      im = b;
    }
    ms = ms+d*d;
  }
  *max = ma;
  *rms = sqrt(ms/ncons);
  *imax = im;
}

bool constrain_lincs(FILE *log,t_inputrec *ir,
		     int step,t_lincsdata *lincsd,t_mdatoms *md,
		     rvec *x,rvec *xprime,rvec *min_proj,matrix box,
		     real lambda,real *dvdlambda,
		     bool bCalcVir,tensor rmdr,
		     bool bCoordinates,
		     t_nrnb *nrnb,bool bDumpOnError)
{
  char  buf[STRLEN];
  int   i,warn,p_imax,error;
  real  p_max,p_rms;
  t_pbc pbc,*pbc_null;
  rvec  dx;
  bool  bOK;
  
  bOK = TRUE;

  if (lincsd->nc == 0)
    return bOK;

  if (ir->ePBC == epbcFULL) {
    /* This is wasting some CPU time as we now do this multiple times
     * per MD step.
     */
    set_pbc_ss(&pbc,box);
    pbc_null = &pbc;
  } else {
    pbc_null = NULL;
  }
  if (bCoordinates) {
    if (ir->efep != efepNO) {
      if (md->bMassPerturbed && lincsd->matlam != lambda) {
	set_lincs_matrix(lincsd,md->invmass);
	lincsd->matlam = lambda;
      }
      
      for(i=0; i<lincsd->nc; i++)
	lincsd->bllen[i] = lincsd->bllen0[i] + lambda*lincsd->ddist[i];
    }
    
    if (lincsd->nzerolen) {
      /* Set the zero lengths to the old lengths */
      if (ir->ePBC == epbcFULL) {
	for(i=0; i<lincsd->nc; i++)
	  if (lincsd->bllen[i] == 0) {
	    pbc_dx(pbc_null,x[lincsd->bla[2*i]],x[lincsd->bla[2*i+1]],dx);
	    lincsd->bllen[i] = norm(dx);
	  }
      } else {
	for(i=0; i<lincsd->nc; i++)
	  if (lincsd->bllen[i] == 0)
	    lincsd->bllen[i] =
	      sqrt(distance2(x[lincsd->bla[2*i]],x[lincsd->bla[2*i+1]]));
      }
    }
    
    if (do_per_step(step,ir->nstlog) || (step < 0))
      cconerr(&p_max,&p_rms,&p_imax,xprime,pbc_null,
	      lincsd->nc,lincsd->bla,lincsd->bllen);
    
    do_lincs(x,xprime,pbc_null,lincsd,md->invmass,
	     ir->nLincsIter,ir->nProjOrder,ir->LincsWarnAngle,&warn,
	     bCalcVir,rmdr);
    
    if (ir->efep != efepNO) {
      real dt_2,dvdl=0;

      dt_2 = 1.0/(ir->delta_t*ir->delta_t);
      for(i=0; (i<lincsd->nc); i++)
	dvdl += lincsd->lambda[i]*dt_2*lincsd->ddist[i];
      *dvdlambda += dvdl;
    }
    
    if (do_per_step(step,ir->nstlog) || (step < 0)) {
      fprintf(stdlog,"   Rel. Constraint Deviation:  Max    between atoms     RMS\n");
      fprintf(stdlog,"       Before LINCS         %.6f %6d %6d   %.6f\n",
	      p_max,lincsd->bla[2*p_imax]+1,lincsd->bla[2*p_imax+1]+1,p_rms);
      cconerr(&p_max,&p_rms,&p_imax,xprime,pbc_null,
	      lincsd->nc,lincsd->bla,lincsd->bllen);
      fprintf(stdlog,"        After LINCS         %.6f %6d %6d   %.6f\n\n",
	      p_max,lincsd->bla[2*p_imax]+1,lincsd->bla[2*p_imax+1]+1,p_rms);
    }
    
    if (warn > 0) {
      if (bDumpOnError) {
	cconerr(&p_max,&p_rms,&p_imax,xprime,pbc_null,
		lincsd->nc,lincsd->bla,lincsd->bllen);
	sprintf(buf,"\nStep %d, time %g (ps)  LINCS WARNING\n"
		"relative constraint deviation after LINCS:\n"
		"max %.6f (between atoms %d and %d) rms %.6f\n",
		step,ir->init_t+step*ir->delta_t,
		p_max,lincsd->bla[2*p_imax]+1,lincsd->bla[2*p_imax+1]+1,
		p_rms);
	fprintf(stdlog,"%s",buf);
	fprintf(stderr,"%s",buf);
	lincs_warning(x,xprime,pbc_null,lincsd->nc,lincsd->bla,lincsd->bllen,
		      ir->LincsWarnAngle);
      }
      bOK = (p_max < 0.5);
    }
    
    if (lincsd->nzerolen) {
      for(i=0; (i<lincsd->nc); i++)
	if (lincsd->bllen0[i] == 0 && lincsd->ddist[i] == 0)
	  lincsd->bllen[i] = 0;
    }
  } 
  else {
    do_lincsp(x,xprime,min_proj,pbc_null,lincsd,md->invmass,ir->nProjOrder);
  }
  
  /* count assuming nit=1 */
  inc_nrnb(nrnb,eNR_LINCS,lincsd->nc);
  inc_nrnb(nrnb,eNR_LINCSMAT,(2+ir->nProjOrder)*lincsd->ncc);
  if (bCalcVir)
    inc_nrnb(nrnb,eNR_CONSTR_VIR,lincsd->nc);

  return bOK;
}
