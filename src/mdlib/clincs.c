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
#include "copyrite.h"
#include "physics.h"
#include "vec.h"
#include "pbc.h"
#include "smalloc.h"
#include "mdrun.h"
#include "nrnb.h"
#include "domdec.h"

typedef struct gmx_lincsdata {
  int  ncg;         /* the global number of constraints */
  int  ncg_flex;    /* the global number of flexible constraints */
  int  ncg_triangle;/* the global number of constraints in triangles */
  int  nIter;       /* the number of iterations */
  int  nOrder;      /* the order of the matrix expansion */
  int  nc;          /* the number of constraints */
  int  nc_alloc;    /* the number we allocated memory for */
  int  ncc;         /* the number of constraint connections */
  int  ncc_alloc;   /* the number we allocated memory for */
  real matlam;      /* the FE lambda value used for filling blc and blmf */
  real *bllen0;     /* the reference distance in topology A */
  real *ddist;      /* the reference distance in top B - the r.d. in top A */
  int  *bla;        /* the atom pairs involved in the constraints */
  real *blc;        /* 1/sqrt(invmass1 + invmass2) */
  int  *blnr;       /* index into blbnb and blmf */
  int  *blbnb;      /* list of constraint connections */
  int  ntriangle;   /* the local number of constraints in triangles */
  int  *triangle;   /* the list of triangle constraints */
  int  *tri_bits;   /* the bits tell if the matrix element should be used */
  int  ncc_triangle;/* the number of constraint connections in triangles */
  real *blmf;       /* matrix of mass factors for constraint connections */
  real *bllen;      /* the reference bond length */
  /* arrays for temporary storage in the LINCS algorithm */
  rvec *tmpv;
  real *tmpncc;
  real *tmp1;
  real *tmp2;
  real *tmp3;
  real *lambda;  /* the Lagrange multipliers */
  /* storage for the constraint RMS relative deviation output */
  real rmsd_data[3];
} t_gmx_lincsdata;

real *lincs_rmsd_data(struct gmx_lincsdata *lincsd)
{
  return lincsd->rmsd_data;
}

real lincs_rmsd(struct gmx_lincsdata *lincsd,bool bSD2)
{
  if (lincsd->rmsd_data[0] > 0)
    return sqrt(lincsd->rmsd_data[bSD2 ? 2 : 1]/lincsd->rmsd_data[0]);
  else
    return 0;
}

static void lincs_matrix_expand(const struct gmx_lincsdata *lincsd,
				const real *blcc,
				real *rhs1,real *rhs2,real *sol)
{
  int  nrec,rec,ncons,b,j,n,nr0,nr1;
  real mvb,*swap;
  int  ntriangle,tb,bits;
  const int *blnr=lincsd->blnr,*blbnb=lincsd->blbnb;
  const int *triangle=lincsd->triangle,*tri_bits=lincsd->tri_bits;

  ncons     = lincsd->nc;
  ntriangle = lincsd->ntriangle;
  nrec      = lincsd->nOrder;

  for(rec=0; rec<nrec; rec++) {
    for(b=0; b<ncons; b++) {
      mvb = 0;
      for(n=blnr[b]; n<blnr[b+1]; n++) {
	j = blbnb[n];
	mvb = mvb + blcc[n]*rhs1[j];
      }
      rhs2[b] = mvb;
      sol[b]  = sol[b] + mvb;
    }
    swap = rhs1;
    rhs1 = rhs2;
    rhs2 = swap;
  } /* nrec*(ncons+2*nrtot) flops */

  if (ntriangle > 0) {
    /* Perform an extra nrec recursions for only the constraints
     * involved in rigid triangles.
     * In this way their accuracy should come close to those of the other
     * constraints, since traingles of constraints can produce eigenvalues
     * around 0.7, while the effective eigenvalue for bond constraints
     * is around 0.4 (and 0.7*0.7=0.5).
     */
    /* We need to copy the temporary array, since only the elements
     * for constraints involved in triangles are updated
     * and then the pointers are swapped.
     */
    for(b=0; b<ncons; b++)
      rhs2[b] = rhs1[b];
    for(rec=0; rec<nrec; rec++) {
      for(tb=0; tb<ntriangle; tb++) {
	b    = triangle[tb];
	bits = tri_bits[tb];
	mvb = 0;
	nr0 = blnr[b];
	nr1 = blnr[b+1];
	for(n=nr0; n<nr1; n++) {
	  if (bits & (1<<(n-nr0))) {
	    j = blbnb[n];
	    mvb = mvb + blcc[n]*rhs1[j];
	  }
	}
	rhs2[b] = mvb;
	sol[b]  = sol[b] + mvb;
      }
      swap = rhs1;
      rhs1 = rhs2;
      rhs2 = swap;
    } /* flops count is missing here */
  }
}

static void do_lincsp(rvec *x,rvec *f,rvec *fp,t_pbc *pbc,
		      struct gmx_lincsdata *lincsd,real *invmass,
		      int econq,real *dvdlambda)
{
  int     b,i,j,k,n;
  real    tmp0,tmp1,tmp2,im1,im2,mvb,rlen,len,wfac,lam;  
  rvec    dx;
  int     ncons,*bla,*blnr,*blbnb;
  rvec    *r;
  real    *blc,*blmf,*blcc,*rhs1,*rhs2,*sol;

  ncons  = lincsd->nc;
  bla    = lincsd->bla;
  r      = lincsd->tmpv;
  blnr   = lincsd->blnr;
  blbnb  = lincsd->blbnb;
  blc    = lincsd->blc;
  blmf   = lincsd->blmf;
  blcc   = lincsd->tmpncc;
  rhs1   = lincsd->tmp1;
  rhs2   = lincsd->tmp2;
  sol    = lincsd->tmp3;

  if (econq != econqForce)
    dvdlambda = NULL;

  /* Compute normalized i-j vectors */
  if (pbc) {
    for(b=0; b<ncons; b++) {
      pbc_dx_aiuc(pbc,x[bla[2*b]],x[bla[2*b+1]],dx);
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
      blcc[n] = blmf[n]*(tmp0*r[k][0] + tmp1*r[k][1] + tmp2*r[k][2]); 
    } /* 6 nr flops */
    mvb = blc[b]*(tmp0*(f[i][0] - f[j][0]) +
		  tmp1*(f[i][1] - f[j][1]) +    
		  tmp2*(f[i][2] - f[j][2]));
    rhs1[b] = mvb;
    sol[b]  = mvb;
    /* 7 flops */
  }
  /* Together: 23*ncons + 6*nrtot flops */
    
  lincs_matrix_expand(lincsd,blcc,rhs1,rhs2,sol);
  /* nrec*(ncons+2*nrtot) flops */
  
  for(b=0;b<ncons;b++) {
    /* With econqDeriv_FlexCon only use the flexible constraints */
    if (econq != econqDeriv_FlexCon ||
	(lincsd->bllen0[b] == 0 && lincsd->ddist[b] == 0)) {
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
      if (dvdlambda) {
	/* This is only correct with forces and invmass=1 */
	*dvdlambda -= mvb*lincsd->ddist[b];
      }
    }
  } /* 16 ncons flops */
}

static void do_lincs(rvec *x,rvec *xp,matrix box,t_pbc *pbc,
		     struct gmx_lincsdata *lincsd,real *invmass,
		     gmx_domdec_t *dd,
		     real wangle,int *warn,
		     real invdt,rvec *v,
		     bool bCalcVir,tensor rmdr)
{
  int     b,i,j,k,n,iter;
  real    tmp0,tmp1,tmp2,im1,im2,mvb,rlen,len,len2,dlen2,wfac,lam;  
  rvec    dx;
  int     ncons,*bla,*blnr,*blbnb;
  rvec    *r;
  real    *blc,*blmf,*bllen,*blcc,*rhs1,*rhs2,*sol,*lambda;
  int     *nlocat;

  ncons  = lincsd->nc;
  bla    = lincsd->bla;
  r      = lincsd->tmpv;
  blnr   = lincsd->blnr;
  blbnb  = lincsd->blbnb;
  blc    = lincsd->blc;
  blmf   = lincsd->blmf;
  bllen  = lincsd->bllen;
  blcc   = lincsd->tmpncc;
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
      pbc_dx_aiuc(pbc,x[bla[2*b]],x[bla[2*b+1]],dx);
      unitv(dx,r[b]);
    }  
    for(b=0; b<ncons; b++) {
      for(n=blnr[b]; n<blnr[b+1]; n++) {
	blcc[n] = blmf[n]*iprod(r[b],r[blbnb[n]]);
      }
      pbc_dx_aiuc(pbc,xp[bla[2*b]],xp[bla[2*b+1]],dx);
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
      for(n=blnr[b]; n<blnr[b+1]; n++) {
	k = blbnb[n];
	blcc[n] = blmf[n]*(tmp0*r[k][0] + tmp1*r[k][1] + tmp2*r[k][2]); 
      } /* 6 nr flops */
      mvb = blc[b]*(tmp0*(xp[i][0] - xp[j][0]) +
		    tmp1*(xp[i][1] - xp[j][1]) +    
		    tmp2*(xp[i][2] - xp[j][2]) - len);
      rhs1[b] = mvb;
      sol[b]  = mvb;
      /* 10 flops */
    }
    /* Together: 26*ncons + 6*nrtot flops */
  }

  lincs_matrix_expand(lincsd,blcc,rhs1,rhs2,sol);
  /* nrec*(ncons+2*nrtot) flops */

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

  for(iter=0; iter<lincsd->nIter; iter++) {
    if (dd && dd->constraints)
      /* Communicate the corrected non-local coordinates */
      dd_move_x_constraints(dd,box,xp,NULL);

    for(b=0;b<ncons;b++) {
      len = bllen[b];
      if (pbc) {
	pbc_dx_aiuc(pbc,xp[bla[2*b]],xp[bla[2*b+1]],dx);
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
    } /* 20*ncons flops */
    
    lincs_matrix_expand(lincsd,blcc,rhs1,rhs2,sol);
    /* nrec*(ncons+2*nrtot) flops */

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
  } /* nit*ncons*(37+9*nrec) flops */

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
   * 26*ncons + 6*nrtot + nrec*(ncons+2*nrtot)
   * + nit * (20*ncons + nrec*(ncons+2*nrtot) + 17 ncons)
   *
   * (26+nrec)*ncons + (6+2*nrec)*nrtot
   * + nit * ((37+nrec)*ncons + 2*nrec*nrtot)
   * if nit=1
   * (63+nrec)*ncons + (6+4*nrec)*nrtot
   */
}

void set_lincs_matrix(struct gmx_lincsdata *li,real *invmass,real lambda)
{
  int i,a1,a2,n,k,sign,center;
  int end,nk,kk;

  for(i=0; (i<li->nc); i++) {
    a1 = li->bla[2*i];
    a2 = li->bla[2*i+1];
    li->blc[i] = invsqrt(invmass[a1] + invmass[a2]);
  }

  /* Construct the coupling coefficient matrix blmf */
  li->ntriangle = 0;
  li->ncc_triangle = 0;
  for(i=0; (i<li->nc); i++) {
    a1 = li->bla[2*i];
    a2 = li->bla[2*i+1];
    for(n=li->blnr[i]; (n<li->blnr[i+1]); n++) {
      k = li->blbnb[n];
      if (a1 == li->bla[2*k] || a2 == li->bla[2*k+1])
	sign = -1;
      else
	sign = 1;
      if (a1 == li->bla[2*k] || a1 == li->bla[2*k+1]) {
	center = a1;
	end    = a2;
      } else {
	center = a2;
	end    = a1;
      }
      li->blmf[n] = sign*invmass[center]*li->blc[i]*li->blc[k];
      if (li->ncg_triangle > 0) {
	/* Look for constraint triangles */
	for(nk=li->blnr[k]; (nk<li->blnr[k+1]); nk++) {
	  kk = li->blbnb[nk];
	  if (kk != i && kk != k &&
	      (li->bla[2*kk] == end || li->bla[2*kk+1] == end)) {
	    if (li->ntriangle == 0 || 
		li->triangle[li->ntriangle-1] < i) {
	      /* Add this constraint to the triangle list */
	      li->triangle[li->ntriangle] = i;
	      li->tri_bits[li->ntriangle] = 0;
	      li->ntriangle++;
	      if (li->blnr[i+1] - li->blnr[i] > sizeof(li->tri_bits[0])*8 - 1)
		gmx_fatal(FARGS,"A constraint is connected to %d constraints, this is more than the %d allowed for constraints participating in triangles",
			  li->blnr[i+1] - li->blnr[i],
			  sizeof(li->tri_bits[0])*8-1);
	    }
	    li->tri_bits[li->ntriangle-1] |= (1<<(n-li->blnr[i]));
	    li->ncc_triangle++;
	  }
	}
      }
    }
  }

  if (debug) {
    fprintf(debug,"Of the %d constraints %d participate in triangles\n",
	    li->nc,li->ntriangle);
    fprintf(debug,"There are %d couplings of which %d in triangles\n",
	    li->ncc,li->ncc_triangle);
  }

  /* Set matlam, so we know with which lambda value the masses have been set */
  li->matlam = lambda;
}

static int count_triangle_constraints(int ncon,t_iatom *ia,t_blocka *at2con)
{
  int  c0,a00,a01,n1,c1,a10,a11,ac1,n2,c2,a20,a21;
  int  ncon_triangle;
  bool bTriangle;

  ncon_triangle = 0;
  for(c0=0; c0<ncon; c0++) {
    bTriangle = FALSE;
    a00 = ia[c0*3+1];
    a01 = ia[c0*3+2];
    for(n1=at2con->index[a01]; n1<at2con->index[a01+1]; n1++) {
      c1 = at2con->a[n1];
      if (c1 != c0) {
	a10 = ia[c1*3+1];
	a11 = ia[c1*3+2];
	if (a10 == a01)
	  ac1 = a11;
	else
	  ac1 = a10;
	for(n2=at2con->index[ac1]; n2<at2con->index[ac1+1]; n2++) {
	  c2 = at2con->a[n2];
	  if (c2 != c0 && c2 != c1) {
	    a20 = ia[c2*3+1];
	    a21 = ia[c2*3+2];
	    if (a20 == a00 || a21 == a00)
	      bTriangle = TRUE;
	  }
	}
      }
    }
    if (bTriangle)
      ncon_triangle++;
  }

  return ncon_triangle;
}

static int int_comp(const void *a,const void *b)
{
  return (*(int *)a) - (*(int *)b);
}

gmx_lincsdata_t init_lincs(FILE *fplog,t_idef *idef,
			   int nflexcon_global,t_blocka *at2con,
			   int bPLINCS,int nIter,int nProjOrder)
{
  struct gmx_lincsdata *li;

  if (fplog)
    fprintf(fplog,"\nInitializing%s LINear Constraint Solver\n",
	    bPLINCS ? " Parallel" : "");

  snew(li,1);

  li->ncg      = idef->il[F_CONSTR].nr/3;
  li->ncg_flex = nflexcon_global;
  li->ncg_triangle  = 
    count_triangle_constraints(li->ncg,idef->il[F_CONSTR].iatoms,at2con);
  
  li->nIter  = nIter;
  li->nOrder = nProjOrder;

  if (bPLINCS || li->ncg_triangle > 0) {
    please_cite(fplog,"Hess2008a");
  } else {
    please_cite(fplog,"Hess97a");
  }

  if (fplog) {
    fprintf(fplog,"The number of constraints is %d\n",li->ncg);
    if (bPLINCS)
      fprintf(fplog,"There are inter charge-group constraints,\n"
	      "will communicate selected coordinates each lincs iteration\n");
    if (li->ncg_triangle > 0)
      fprintf(fplog,
	      "%d constraints are involved in constraint triangles,\n"
	      "will apply an additional matrix expansion of order %d for couplings\n"
	      "between constraints inside triangles\n",
	      li->ncg_triangle,li->nOrder);
  }
  
  return li;
}

void set_lincs(t_idef *idef,int start,int homenr,
	       t_blocka *at2con,
	       bool bDynamics,gmx_domdec_t *dd,
	       struct gmx_lincsdata *li)
{
  gmx_domdec_constraints_t *dc;
  t_iatom     *iatom;
  int         i,k,ncc_alloc,ni,con,cong,nconnect,concon;
  int         type,ag1,ag2,a1,a2;
  real        lenA=0,lenB;
  bool        bLocal;

  li->nc = 0;
  li->ncc = 0;
  
  if (idef->il[F_CONSTR].nr > 0 || (dd && dd->constraints)) {
    if (debug)
      fprintf(debug,"Building the LINCS connectivity\n");
    if (dd == NULL) {
      li->nc = li->ncg;
      if (!bDynamics)
	li->nc -= li->ncg_flex;
      dc = NULL;
    } else {
      dc = dd->constraints;
      li->nc = dc->ncon;
    }
    
    if (li->nc > li->nc_alloc || li->nc_alloc == 0) {
      li->nc_alloc = over_alloc_dd(li->nc);
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
      if (li->ncg_triangle > 0) {
	/* This is allocating too much, but it is difficult too improve */
	srenew(li->triangle,li->nc_alloc);
	srenew(li->tri_bits,li->nc_alloc);
      }
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
	for(k=at2con->index[ag1-start]; k<at2con->index[ag1-start+1]; k++) {
	  concon = at2con->a[k];
	  if (concon != cong &&
	      (bLocal || dc->gc2lc[concon] != -1)) {
	    if (nconnect >= ncc_alloc) {
	      ncc_alloc = over_alloc_small(nconnect+1);
	      srenew(li->blbnb,ncc_alloc);
	    }
	    if (dc == NULL)
	      li->blbnb[nconnect++] = concon;
	    else
	      li->blbnb[nconnect++] = dc->gc2lc[concon];
	  }
	}
	for(k=at2con->index[ag2-start]; k<at2con->index[ag2-start+1]; k++) {
	  concon = at2con->a[k];
	  if (concon != cong &&
	      (bLocal || dc->gc2lc[concon] != -1)) {
	    if (nconnect+1 > ncc_alloc) {
	      ncc_alloc = over_alloc_small(nconnect+1);
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
    if (dd == NULL) {
      /* Since the matrix is static, we can free some memory */
      ncc_alloc = li->ncc;
      srenew(li->blbnb,ncc_alloc);
    }

    if (ncc_alloc > li->ncc_alloc) {
      li->ncc_alloc = ncc_alloc;
      srenew(li->blmf,li->ncc_alloc);
      srenew(li->tmpncc,li->ncc_alloc);
    }

    if (debug)
      fprintf(debug,"Number of constraints is %d, couplings %d\n",
	      li->nc,li->ncc);
  }
}

static void lincs_warning(FILE *fplog,
			  gmx_domdec_t *dd,rvec *x,rvec *xprime,t_pbc *pbc,
			  int ncons,int *bla,real *bllen,real wangle,
			  int maxwarn,int *warncount)
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
  if (fplog)
    fprintf(fplog,buf); 

  for(b=0;b<ncons;b++) {
    i = bla[2*b];
    j = bla[2*b+1];
    if (pbc) {
      pbc_dx_aiuc(pbc,x[i],x[j],v0);
      pbc_dx_aiuc(pbc,xprime[i],xprime[j],v1);
    } else {
      rvec_sub(x[i],x[j],v0);
      rvec_sub(xprime[i],xprime[j],v1);
    }
    d0 = norm(v0);
    d1 = norm(v1);
    cosine = iprod(v0,v1)/(d0*d1);
    if (cosine < wfac) {
      sprintf(buf," %6d %6d  %5.1f  %8.4f %8.4f    %8.4f\n",
	      glatnr(dd,i),glatnr(dd,j),RAD2DEG*acos(cosine),d0,d1,bllen[b]);
      fprintf(stderr,buf);
      if (fplog)
	fprintf(fplog,buf);
      (*warncount)++;
    }
  }
  if (*warncount > maxwarn)
    too_many_constraint_warnings(econtLINCS,*warncount);
}

static void cconerr(gmx_domdec_t *dd,
		    int ncons,int *bla,real *bllen,rvec *x,t_pbc *pbc,
		    real *ncons_loc,real *ssd,real *max,int *imax)
{
  real      len,d,ma,ssd2,r2;
  int       *nlocat,count,b,im;
  rvec      dx;
  
  if (dd && dd->constraints)
    nlocat = dd->constraints->con_nlocat;
  else
    nlocat = 0;

  ma = 0;
  ssd2 = 0;
  im = 0;
  count = 0;
  for(b=0;b<ncons;b++) {
    if (pbc) {
      pbc_dx_aiuc(pbc,x[bla[2*b]],x[bla[2*b+1]],dx);
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
      ssd2 += d*d;
      count++;
    } else {
      ssd2 += nlocat[b]*d*d;
      count += nlocat[b];
    }
  }

  *ncons_loc = (nlocat ? 0.5 : 1)*count;
  *ssd       = (nlocat ? 0.5 : 1)*ssd2;
  *max = ma;
  *imax = im;
}

static void dump_conf(gmx_domdec_t *dd,struct gmx_lincsdata *li,
		      t_blocka *at2con,
		      char *name,bool bAll,rvec *x,matrix box)
{
  char str[STRLEN];
  FILE *fp;
  int  ac0,ac1,i,j;
  gmx_domdec_constraints_t *dc;
  bool bPrint=TRUE;
  
  dc = dd->constraints;

  dd_get_constraint_range(dd,&ac0,&ac1);

  sprintf(str,"%s_%d_%d_%d.pdb",name,dd->ci[XX],dd->ci[YY],dd->ci[ZZ]);
  fp = ffopen(str,"w");
  fprintf(fp,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
	  10*norm(box[XX]),10*norm(box[YY]),10*norm(box[ZZ]),
	  90.0,90.0,90.0);
  for(i=0; i<ac1; i++) {
    if (i<dd->nat_home ||
	(bAll && i>=dd->nat_tot && dc->ga2la[dd->gatindex[i]]>=0)) {
      if (at2con) {
	bPrint = (i < dd->nat_home);
	if (i>=dd->nat_tot) {
	  for(j=at2con->index[dd->gatindex[i]+1]; j<at2con->index[dd->gatindex[i]+1]; j++)
	    bPrint = bPrint || (dc->gc2lc[at2con->a[j]]>=0);
	}
      } else {
	bPrint = TRUE;
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

bool constrain_lincs(FILE *fplog,bool bLog,bool bEner,
		     t_inputrec *ir,
		     int step,struct gmx_lincsdata *lincsd,t_mdatoms *md,
		     gmx_domdec_t *dd,
		     rvec *x,rvec *xprime,rvec *min_proj,matrix box,
		     real lambda,real *dvdlambda,
		     real invdt,rvec *v,
		     bool bCalcVir,tensor rmdr,
		     int econq,
		     t_nrnb *nrnb,
		     int maxwarn,int *warncount)
{
  char  buf[STRLEN];
  int   i,warn,p_imax,error;
  real  ncons_loc,p_ssd,p_max;
  t_pbc pbc,*pbc_null;
  rvec  dx;
  bool  bOK;
  
  bOK = TRUE;

  if (lincsd->nc == 0 && dd == NULL) {
    if (bLog || bEner) {
      lincsd->rmsd_data[0] = 0;
      if (ir->eI == eiSD2 && v == NULL)
	i = 2;
      else
	i = 1;
      lincsd->rmsd_data[i] = 0;
    }
    
    return bOK;
  }

  /* We do not need full pbc when constraints do not cross charge groups,
   * i.e. when dd->constraint_comm==NULL
   */
  if ((dd || ir->bPeriodicMols) && !(dd && dd->constraint_comm==NULL)) {
    /* With pbc=screw the screw has been changed to a shift
     * by the constraint coordinate communication routine,
     * so that here we can use normal pbc.
     */
    pbc_null = set_pbc_dd(&pbc,ir->ePBC,dd,FALSE,box);
  } else {
    pbc_null = NULL;
  }
  if (dd) {
    /* Communicate the coordinates required for the non-local constraints */
    dd_move_x_constraints(dd,box,x,xprime);
    /* dump_conf(dd,lincsd,NULL,"con",TRUE,xprime,box); */
  }
  if (econq == econqCoord) {
    if (ir->efep != efepNO) {
      if (md->nMassPerturbed && lincsd->matlam != md->lambda) {
	set_lincs_matrix(lincsd,md->invmass,md->lambda);
      }
      
      for(i=0; i<lincsd->nc; i++)
	lincsd->bllen[i] = lincsd->bllen0[i] + lambda*lincsd->ddist[i];
    }
    
    if (lincsd->ncg_flex) {
      /* Set the flexible constraint lengths to the old lengths */
      if (pbc_null) {
	for(i=0; i<lincsd->nc; i++)
	  if (lincsd->bllen[i] == 0) {
	    pbc_dx_aiuc(pbc_null,x[lincsd->bla[2*i]],x[lincsd->bla[2*i+1]],dx);
	    lincsd->bllen[i] = norm(dx);
	  }
      } else {
	for(i=0; i<lincsd->nc; i++)
	  if (lincsd->bllen[i] == 0)
	    lincsd->bllen[i] =
	      sqrt(distance2(x[lincsd->bla[2*i]],x[lincsd->bla[2*i+1]]));
      }
    }
    
    if (bLog && fplog)
      cconerr(dd,lincsd->nc,lincsd->bla,lincsd->bllen,xprime,pbc_null,
	      &ncons_loc,&p_ssd,&p_max,&p_imax);
    
    do_lincs(x,xprime,box,pbc_null,lincsd,md->invmass,dd,
	     ir->LincsWarnAngle,&warn,
	     invdt,v,bCalcVir,rmdr);

    if (ir->efep != efepNO) {
      real dt_2,dvdl=0;

      dt_2 = 1.0/(ir->delta_t*ir->delta_t);
      for(i=0; (i<lincsd->nc); i++)
	dvdl += lincsd->lambda[i]*dt_2*lincsd->ddist[i];
      *dvdlambda += dvdl;
    }
    
    if (bLog && fplog && lincsd->nc > 0) {
      fprintf(fplog,"   Rel. Constraint Deviation:  RMS         MAX     between atoms\n");
      fprintf(fplog,"       Before LINCS          %.6f    %.6f %6d %6d\n",
	      sqrt(p_ssd/ncons_loc),p_max,
	      glatnr(dd,lincsd->bla[2*p_imax]),
	      glatnr(dd,lincsd->bla[2*p_imax+1]));
    }
    if (bLog || bEner) {
      cconerr(dd,lincsd->nc,lincsd->bla,lincsd->bllen,xprime,pbc_null,
	      &ncons_loc,&p_ssd,&p_max,&p_imax);
      lincsd->rmsd_data[0] = ncons_loc;
      /* Check if we are doing the second part of SD */
      if (ir->eI == eiSD2 && v == NULL)
	i = 2;
      else
	i = 1;
      lincsd->rmsd_data[0] = ncons_loc;
      lincsd->rmsd_data[i] = p_ssd;
    }
    if (bLog && fplog && lincsd->nc > 0) {
      fprintf(fplog,"        After LINCS          %.6f    %.6f %6d %6d\n\n",
	      sqrt(p_ssd/ncons_loc),p_max,
	      glatnr(dd,lincsd->bla[2*p_imax]),
	      glatnr(dd,lincsd->bla[2*p_imax+1]));
    }
    
    if (warn > 0) {
      if (maxwarn >= 0) {
	cconerr(dd,lincsd->nc,lincsd->bla,lincsd->bllen,xprime,pbc_null,
		&ncons_loc,&p_ssd,&p_max,&p_imax);
	sprintf(buf,"\nStep %d, time %g (ps)  LINCS WARNING\n"
		"relative constraint deviation after LINCS:\n"
		"rms %.6f, max %.6f (between atoms %d and %d)\n",
		step,ir->init_t+step*ir->delta_t,
		sqrt(p_ssd/ncons_loc),p_max,
		glatnr(dd,lincsd->bla[2*p_imax]),
		glatnr(dd,lincsd->bla[2*p_imax+1]));
	if (fplog)
	  fprintf(fplog,"%s",buf);
	fprintf(stderr,"%s",buf);
	lincs_warning(fplog,dd,x,xprime,pbc_null,
		      lincsd->nc,lincsd->bla,lincsd->bllen,
		      ir->LincsWarnAngle,maxwarn,warncount);
      }
      bOK = (p_max < 0.5);
    }
    
    if (lincsd->ncg_flex) {
      for(i=0; (i<lincsd->nc); i++)
	if (lincsd->bllen0[i] == 0 && lincsd->ddist[i] == 0)
	  lincsd->bllen[i] = 0;
    }
  } 
  else {
    do_lincsp(x,xprime,min_proj,pbc_null,lincsd,md->invmass,econq,dvdlambda);
  }
  
  /* count assuming nit=1 */
  inc_nrnb(nrnb,eNR_LINCS,lincsd->nc);
  inc_nrnb(nrnb,eNR_LINCSMAT,(2+lincsd->nOrder)*lincsd->ncc);
  if (lincsd->ntriangle > 0)
    inc_nrnb(nrnb,eNR_LINCSMAT,lincsd->nOrder*lincsd->ncc_triangle);
  if (v)
    inc_nrnb(nrnb,eNR_CONSTR_V,lincsd->nc*2);
  if (bCalcVir)
    inc_nrnb(nrnb,eNR_CONSTR_VIR,lincsd->nc);

  return bOK;
}
