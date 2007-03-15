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
 *
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
#include "smalloc.h"
#include "mdrun.h"
#include "nrnb.h"
#include "domdec.h"

static void do_lincsp(rvec *x,rvec *f,rvec *fp,t_pbc *pbc,
		      t_lincsdata *lincsd,real *invmass,
		      int nrec)
{
  int     b,i,j,k,n,it,rec;
  real    tmp0,tmp1,tmp2,im1,im2,mvb,rlen,len,wfac,lam;  
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
    fp[i][0] -= tmp0*im1;
    fp[i][1] -= tmp1*im1;
    fp[i][2] -= tmp2*im1;
    fp[j][0] += tmp0*im2;
    fp[j][1] += tmp1*im2;
    fp[j][2] += tmp2*im2;
  } /* 16 ncons flops */
}

static void do_lincs(rvec *x,rvec *xp,matrix box,t_pbc *pbc,
		     t_lincsdata *lincsd,real *invmass,
		     gmx_domdec_t *dd,
		     int nit,int nrec,
		     real wangle,int *warn,
		     real invdt,rvec *v,
		     bool bCalcVir,tensor rmdr)
{
  int     b,i,j,k,n,it,rec;
  real    tmp0,tmp1,tmp2,im1,im2,mvb,rlen,len,len2,dlen2,wfac,lam;  
  rvec    dx;
  int     ncons,*bla,*blnr,*blbnb;
  rvec    *r;
  real    *blc,*blcc,*bllen,*blm,*rhs1,*rhs2,*sol,*lambda,*swap;
  int     *nlocat;

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

  if (dd && dd->constraints)
    nlocat = dd->constraints->con_nlocat;
  else
    nlocat = NULL;

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
    xp[i][0] -= tmp0*im1;
    xp[i][1] -= tmp1*im1;
    xp[i][2] -= tmp2*im1;
    xp[j][0] += tmp0*im2;
    xp[j][1] += tmp1*im2;
    xp[j][2] += tmp2*im2;
  } /* 16 ncons flops */
  
  
  
  /*     
  ********  Correction for centripetal effects  ********  
  */
  
  wfac = cos(DEG2RAD*wangle);
  wfac = wfac*wfac;

  for(it=0; it<nit; it++) {
    if (dd && dd->constraints)
      /* Communicate the corrected non-local coordinates */
      dd_move_x_constraints(dd,box,xp,NULL);

    for(b=0;b<ncons;b++) {
      len = bllen[b];
      if (pbc) {
	pbc_dx(pbc,xp[bla[2*b]],xp[bla[2*b+1]],dx);
      } else {
	rvec_sub(xp[bla[2*b]],xp[bla[2*b+1]],dx);
      }
      len2 = len*len;
      dlen2 = 2*len2 - norm2(dx);
      if (dlen2 < wfac*len2 && (nlocat==NULL || nlocat[b]))
	*warn = b;
      if (dlen2 < 0)
	dlen2 = 0;
      mvb = blc[b]*(len - dlen2*invsqrt(dlen2));
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
      xp[i][0] -= tmp0*im1;
      xp[i][1] -= tmp1*im1;
      xp[i][2] -= tmp2*im1;
      xp[j][0] += tmp0*im2;
      xp[j][1] += tmp1*im2;
      xp[j][2] += tmp2*im2;
    } /* 17 ncons flops */
  } /* nit*ncons*(35+9*nrec) flops */

  if (v) {
    /* Correct the velocities */
    for(b=0; b<ncons; b++) {
      i = bla[2*b];
      j = bla[2*b+1];
      im1 = invmass[i]*lambda[b]*invdt;
      im2 = invmass[j]*lambda[b]*invdt;
      v[i][0] += im1*r[b][0];
      v[i][1] += im1*r[b][1];
      v[i][2] += im1*r[b][2];
      v[j][0] -= im2*r[b][0];
      v[j][1] -= im2*r[b][1];
      v[j][2] -= im2*r[b][2];
    } /* 16 ncons flops */
  }

  if (nlocat) {
    /* Only account for local atoms */
    for(b=0; b<ncons; b++)
      lambda[b] *= 0.5*nlocat[b];
  }

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

t_block make_at2con(int start,int natoms,
		    t_idef *idef,bool bDynamics,
		    int *nconstraints,int *nflexiblecons)
{
  int *count,ncon,con,nflexcon,i,a;
  t_iatom *ia;
  t_block at2con;
  bool bFlexCon;
  
  snew(count,natoms);
  ncon = idef->il[F_CONSTR].nr/3;
  ia   = idef->il[F_CONSTR].iatoms;
  nflexcon = 0;
  for(con=0; con<ncon; con++) {
    bFlexCon = (idef->iparams[ia[0]].constr.dA == 0 &&
		idef->iparams[ia[0]].constr.dB == 0);
    if (bFlexCon)
      nflexcon++;
    if (bDynamics || !bFlexCon) {
      for(i=1; i<3; i++) {
	a = ia[i] - start;
	count[a]++;
      }
    }
    ia += 3;
  }
  *nconstraints = ncon;
  if (!bDynamics)
    *nconstraints -= nflexcon;
  *nflexiblecons = nflexcon;

  snew(at2con.index,natoms+1);
  at2con.index[0] = 0;
  for(a=0; a<natoms; a++) {
    at2con.index[a+1] = at2con.index[a] + count[a];
    count[a] = 0;
  }
  snew(at2con.a,at2con.index[natoms]);

  ia   = idef->il[F_CONSTR].iatoms;
  for(con=0; con<ncon; con++) {
    bFlexCon = (idef->iparams[ia[0]].constr.dA == 0 &&
		idef->iparams[ia[0]].constr.dB == 0);
    if (bDynamics || !bFlexCon) {
      for(i=1; i<3; i++) {
	a = ia[i] - start;
	at2con.a[at2con.index[a]+count[a]++] = con;
      }
    }
    ia += 3;
  }
  
  sfree(count);

  return at2con;
}

static int int_comp(const void *a,const void *b)
{
  return (*(int *)a) - (*(int *)b);
}

#define CONCON_ALLOC_SIZE 1000

void init_lincs(FILE *log,t_idef *idef,int start,int homenr,
		bool bDynamics,gmx_domdec_t *dd,
		t_lincsdata *li)
{
  gmx_domdec_constraints_t *dc;
  t_block     at2con;
  t_iatom     *iatom;
  int         i,k,ncc_alloc,ni,con,cong,nconnect,concon;
  int         type,ag1,ag2,a1,a2;
  real        lenA=0,lenB;
  bool        bLocal;

  li->nc = 0;
  li->ncc = 0;
  
  if (idef->il[F_CONSTR].nr > 0 || (dd && dd->constraints)) {
    if (debug)
      fprintf(debug,"Initializing LINCS\n");
    if (dd == NULL) {
      dc = NULL;
      /* Make atom-constraint connection list for temporary use */
      at2con = make_at2con(start,homenr,idef,bDynamics,&li->nc,&li->nflexcon);
    } else {
      dc = dd->constraints;
      at2con = dc->at2con;
      li->nc = dc->ncon;
      li->nflexcon = dc->nflexcon_global;
    }
    
    if (li->nc > li->nc_alloc || li->nc_alloc == 0) {
      li->nc_alloc = over_alloc(li->nc);
      srenew(li->bllen0,li->nc_alloc);
      srenew(li->ddist,li->nc_alloc);
      srenew(li->bla,2*li->nc_alloc);
      srenew(li->blc,li->nc_alloc);
      srenew(li->blnr,li->nc_alloc+1);
      srenew(li->bllen,li->nc_alloc);
      srenew(li->tmpv,li->nc_alloc);
      srenew(li->tmp1,li->nc_alloc);
      srenew(li->tmp2,li->nc_alloc);
      srenew(li->tmp3,li->nc_alloc);
      srenew(li->lambda,li->nc_alloc);
    }
    
    if (dc == NULL)
      iatom = idef->il[F_CONSTR].iatoms;
    else
      iatom = dc->iatoms;

    ncc_alloc = li->ncc_alloc;
    li->blnr[0] = 0;

    if (dc == NULL) {
      ni = idef->il[F_CONSTR].nr/3;
    } else {
      ni = dc->ncon;
    }
    con = 0;
    nconnect = 0;
    li->blnr[con] = nconnect;
    for(i=0; i<ni; i++) {
      bLocal = TRUE;
      if (dc == NULL) {
	cong = i;
      } else {
	cong = dc->con[i];
      }
      type = iatom[3*cong];
      ag1  = iatom[3*cong+1];
      ag2  = iatom[3*cong+2];
      if (dc == NULL) {
	a1 = ag1;
	a2 = ag2;
      } else {
	/* Find the global atom indices */
	if (dd->ga2la[ag1].cell == 0) {
	  a1 = dd->ga2la[ag1].a;
	} else {
	  a1 = dc->ga2la[ag1];
	  bLocal = FALSE;
	}
	if (dd->ga2la[ag2].cell == 0) {
	  a2 = dd->ga2la[ag2].a;
	} else {
	  a2 = dc->ga2la[ag2];
	  bLocal = FALSE;
	}
      }
      lenA = idef->iparams[type].constr.dA;
      lenB = idef->iparams[type].constr.dB;
      /* Skip the flexible constraints when not doing dynamics */
      if (bDynamics || lenA!=0 || lenB!=0) {
	li->bllen0[con]  = lenA;
	li->ddist[con]   = lenB - lenA;
	/* Set the length to the topology A length */
	li->bllen[con]   = li->bllen0[con];
	li->bla[2*con]   = a1;
	li->bla[2*con+1] = a2;
	/* Construct the constraint connection matrix blbnb */
	for(k=at2con.index[ag1-start]; k<at2con.index[ag1-start+1]; k++) {
	  concon = at2con.a[k];
	  if (concon != cong &&
	      (bLocal || dc->gc2lc[concon] != -1)) {
	    if (nconnect >= ncc_alloc) {
	      ncc_alloc += CONCON_ALLOC_SIZE;
	      srenew(li->blbnb,ncc_alloc);
	    }
	    if (dc == NULL)
	      li->blbnb[nconnect++] = concon;
	    else
	      li->blbnb[nconnect++] = dc->gc2lc[concon];
	  }
	}
	for(k=at2con.index[ag2-start]; k<at2con.index[ag2-start+1]; k++) {
	  concon = at2con.a[k];
	  if (concon != cong &&
	      (bLocal || dc->gc2lc[concon] != -1)) {
	    if (nconnect >= ncc_alloc) {
	      ncc_alloc += CONCON_ALLOC_SIZE;
	      srenew(li->blbnb,ncc_alloc);
	    }
	    if (dc == NULL)
	      li->blbnb[nconnect++] = concon;
	    else
	      li->blbnb[nconnect++] = dc->gc2lc[concon];
	  }
	}
	li->blnr[con+1] = nconnect;

	if (dc == NULL) {
	  /* Order the blbnb matrix to optimize memory access */
	  qsort(&(li->blbnb[li->blnr[con]]),li->blnr[con+1]-li->blnr[con],
		sizeof(li->blbnb[0]),int_comp);
	}
	/* Increase the constraint count */
	con++;
      }
    }
    if (con != li->nc)
      gmx_fatal(FARGS,"Internal error in init_lincs, old (%d) and new (%d) constraint count don't match\n",li->nc,con);
    li->ncc = li->blnr[con];

    if (ncc_alloc > li->ncc_alloc) {
      li->ncc_alloc = ncc_alloc;
      srenew(li->blcc,li->ncc_alloc);
      srenew(li->tmpncc,li->ncc_alloc);
    }

    if (dc == NULL) {
      done_block(&at2con);
      fprintf(log,"\nInitializing LINear Constraint Solver\n");
      fprintf(log,"  number of constraints is %d\n",li->nc);
      if (li->nc > 0)
	fprintf(log,"  average number of constraints coupled to one constraint is %.1f\n",
		(real)(li->ncc)/li->nc);
      if (li->nflexcon)
	fprintf(log,"  found %d flexible constraints\n",li->nflexcon);
      fprintf(log,"\n");
      fflush(log);
    }
    if (debug)
      fprintf(debug,"Number of constraints is %d\n",li->nc);
  }
}

/* Yes. I _am_ awfully ashamed of introducing a static variable after preaching
 * that we should get rid of them. Do as I say, not as I do!
 *
 * In defense, this is a hack just before the 3.3 release. The proper solution
 * is to introduce an abstract constraint object where it is stored. /EL 20050830
 */
static int lincs_warncount = 0;

#define LINCS_MAXWARN 10000

static void lincs_warning(gmx_domdec_t *dd,rvec *x,rvec *xprime,t_pbc *pbc,
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
	      glatnr(dd,i),glatnr(dd,j),RAD2DEG*acos(cosine),d0,d1,bllen[b]);
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

static void cconerr(gmx_domdec_t *dd,
		    real *max,real *rms,int *imax,rvec *x,t_pbc *pbc,
		    int ncons,int *bla,real *bllen)
{
  real      len,d,ma,ms,r2;
  int       *nlocat,count,b,im;
  rvec      dx;
  
  if (dd && dd->constraints)
    nlocat = dd->constraints->con_nlocat;
  else
    nlocat = 0;

  ma = 0;
  ms = 0;
  im = 0;
  count = 0;
  for(b=0;b<ncons;b++) {
    if (pbc) {
      pbc_dx(pbc,x[bla[2*b]],x[bla[2*b+1]],dx);
    } else {
      rvec_sub(x[bla[2*b]],x[bla[2*b+1]],dx);
    }
    r2 = norm2(dx);
    len = r2*invsqrt(r2);
    d = fabs(len/bllen[b]-1);
    if (d > ma && (nlocat==NULL || nlocat[b])) {
      ma = d;
      im = b;
    }
    if (nlocat == NULL) {
      ms = ms + d*d;
      count++;
    } else {
      ms = ms + nlocat[b]*d*d;
      count += nlocat[b];
    }
  }
  *max = ma;
  *rms = sqrt(ms/count);
  *imax = im;
}

static void dump_conf(gmx_domdec_t *dd,t_lincsdata *li,
		      char *name,bool bAll,rvec *x,matrix box)
{
  char str[STRLEN];
  FILE *fp;
  int i,j;
  gmx_domdec_constraints_t *dc;
  bool bPrint;
  
  dc = dd->constraints;

  sprintf(str,"%s%d.pdb",name,dd->sim_nodeid);
  fp = ffopen(str,"w");
  fprintf(fp,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
	  10*norm(box[XX]),10*norm(box[YY]),10*norm(box[ZZ]),
	  90.0,90.0,90.0);
  for(i=0; i<dd->nat_tot_con; i++) {
    if (i<dd->nat_home ||
	(bAll && i>=dd->nat_tot && dc->ga2la[dd->gatindex[i]]>=0)) {
      bPrint = (i < dd->nat_home);
      if (i>=dd->nat_tot) {
	for(j=dc->at2con.index[dd->gatindex[i]+1]; j<dc->at2con.index[dd->gatindex[i]+1]; j++)
	  bPrint = bPrint || (dc->gc2lc[dc->at2con.a[j]]>=0);
      }
      if (bPrint)
	fprintf(fp,"%-6s%5u  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
		"ATOM",glatnr(dd,i),"C","ALA",' ',i+1,
		10*x[i][XX],10*x[i][YY],10*x[i][ZZ],
		1.0,i<dd->nat_tot ? 0.0 : 1.0);
    }
  }
  if (bAll) {
    for(i=0; i<li->nc; i++)
      fprintf(fp,"CONECT%5d%5d\n",
	      glatnr(dd,li->bla[2*i]),
	      glatnr(dd,li->bla[2*i+1]));
  }
  fclose(fp);
}

bool constrain_lincs(FILE *log,bool bLog,
		     t_inputrec *ir,
		     int step,t_lincsdata *lincsd,t_mdatoms *md,
		     gmx_domdec_t *dd,
		     rvec *x,rvec *xprime,rvec *min_proj,matrix box,
		     real lambda,real *dvdlambda,
		     real invdt,rvec *v,
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

  if (lincsd->nc == 0 && dd==NULL)
    return bOK;

  /* We do not need full pbc when constraints do not cross charge groups,
   * i.e. when dd->constraint_comm==NULL
   */
  if ((dd || ir->bPeriodicMols) && !(dd && dd->constraint_comm==NULL)) {
    /* This is wasting some CPU time as we now do this multiple times
     * per MD step.
     */
    pbc_null = set_pbc_ss(&pbc,ir->ePBC,box,dd,FALSE);
  } else {
    pbc_null = NULL;
  }
  if (dd) {
    /* Communicate the coordinates required for the non-local constraints */
    dd_move_x_constraints(dd,box,x,xprime);
    /*dump_conf(dd,lincsd,"con",TRUE,xprime,box);*/
  }
  if (bCoordinates) {
    if (ir->efep != efepNO) {
      if (md->nMassPerturbed && lincsd->matlam != lambda) {
	set_lincs_matrix(lincsd,md->invmass);
	lincsd->matlam = lambda;
      }
      
      for(i=0; i<lincsd->nc; i++)
	lincsd->bllen[i] = lincsd->bllen0[i] + lambda*lincsd->ddist[i];
    }
    
    if (lincsd->nflexcon) {
      /* Set the flexible constraint lengths to the old lengths */
      if (pbc_null) {
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
    
    if (bLog)
      cconerr(dd,&p_max,&p_rms,&p_imax,xprime,pbc_null,
	      lincsd->nc,lincsd->bla,lincsd->bllen);
    
    do_lincs(x,xprime,box,pbc_null,lincsd,md->invmass,dd,
	     ir->nLincsIter,ir->nProjOrder,ir->LincsWarnAngle,&warn,
	     invdt,v,bCalcVir,rmdr);

    if (ir->efep != efepNO) {
      real dt_2,dvdl=0;

      dt_2 = 1.0/(ir->delta_t*ir->delta_t);
      for(i=0; (i<lincsd->nc); i++)
	dvdl += lincsd->lambda[i]*dt_2*lincsd->ddist[i];
      *dvdlambda += dvdl;
    }
    
    if (bLog && lincsd->nc > 0) {
      fprintf(stdlog,"   Rel. Constraint Deviation:  Max    between atoms     RMS\n");
      fprintf(stdlog,"       Before LINCS         %.6f %6d %6d   %.6f\n",
	      p_max,glatnr(dd,lincsd->bla[2*p_imax]),
	      glatnr(dd,lincsd->bla[2*p_imax+1]),p_rms);
      cconerr(dd,&p_max,&p_rms,&p_imax,xprime,pbc_null,
	      lincsd->nc,lincsd->bla,lincsd->bllen);
      fprintf(stdlog,"        After LINCS         %.6f %6d %6d   %.6f\n\n",
	      p_max,glatnr(dd,lincsd->bla[2*p_imax]),
	      glatnr(dd,lincsd->bla[2*p_imax+1]),p_rms);
    }
    
    if (warn > 0) {
      if (bDumpOnError) {
	cconerr(dd,&p_max,&p_rms,&p_imax,xprime,pbc_null,
		lincsd->nc,lincsd->bla,lincsd->bllen);
	sprintf(buf,"\nStep %d, time %g (ps)  LINCS WARNING\n"
		"relative constraint deviation after LINCS:\n"
		"max %.6f (between atoms %d and %d) rms %.6f\n",
		step,ir->init_t+step*ir->delta_t,
		p_max,glatnr(dd,lincsd->bla[2*p_imax]),
		glatnr(dd,lincsd->bla[2*p_imax+1]),p_rms);
	fprintf(stdlog,"%s",buf);
	fprintf(stderr,"%s",buf);
	lincs_warning(dd,x,xprime,pbc_null,
		      lincsd->nc,lincsd->bla,lincsd->bllen,
		      ir->LincsWarnAngle);
      }
      bOK = (p_max < 0.5);
    }
    
    if (lincsd->nflexcon) {
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
  if (v)
    inc_nrnb(nrnb,eNR_CONSTR_V,lincsd->nc*2);
  if (bCalcVir)
    inc_nrnb(nrnb,eNR_CONSTR_VIR,lincsd->nc);

  return bOK;
}
