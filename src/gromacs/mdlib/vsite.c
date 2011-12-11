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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "vsite.h"
#include "macros.h"
#include "smalloc.h"
#include "nrnb.h"
#include "vec.h"
#include "mvdata.h"
#include "network.h"
#include "mshift.h"
#include "pbc.h"
#include "domdec.h"
#include "partdec.h"
#include "mtop_util.h"

/* Routines to send/recieve coordinates and force
 * of constructing atoms. 
 */ 

static void move_construct_x(t_comm_vsites *vsitecomm, rvec x[], t_commrec *cr)
{
	rvec *sendbuf;
	rvec *recvbuf;
	int i,ia;
	
	sendbuf = vsitecomm->send_buf;
	recvbuf = vsitecomm->recv_buf;
	

		/* Prepare pulse left by copying to send buffer */
		for(i=0;i<vsitecomm->left_export_nconstruct;i++)
		{
			ia = vsitecomm->left_export_construct[i];
			copy_rvec(x[ia],sendbuf[i]);
		}
	
		/* Pulse coordinates left */
		gmx_tx_rx_real(cr,GMX_LEFT,(real *)sendbuf,3*vsitecomm->left_export_nconstruct,GMX_RIGHT,(real *)recvbuf,3*vsitecomm->right_import_nconstruct);
		
		/* Copy from receive buffer to coordinate array */
		for(i=0;i<vsitecomm->right_import_nconstruct;i++)
		{
			ia = vsitecomm->right_import_construct[i];
			copy_rvec(recvbuf[i],x[ia]);
		}

		/* Prepare pulse right by copying to send buffer */
		for(i=0;i<vsitecomm->right_export_nconstruct;i++)
		{
			ia = vsitecomm->right_export_construct[i];
			copy_rvec(x[ia],sendbuf[i]);
		}
		
		/* Pulse coordinates right */
		gmx_tx_rx_real(cr,GMX_RIGHT,(real *)sendbuf,3*vsitecomm->right_export_nconstruct,GMX_LEFT,(real *)recvbuf,3*vsitecomm->left_import_nconstruct);
		
		/* Copy from receive buffer to coordinate array */
		for(i=0;i<vsitecomm->left_import_nconstruct;i++)
		{
			ia = vsitecomm->left_import_construct[i];
			copy_rvec(recvbuf[i],x[ia]);
		}
}


static void move_construct_f(t_comm_vsites *vsitecomm, rvec f[], t_commrec *cr)
{
	rvec *sendbuf;
	rvec *recvbuf;
	int i,ia;

	sendbuf = vsitecomm->send_buf;
	recvbuf = vsitecomm->recv_buf;	

		/* Prepare pulse right by copying to send buffer */
		for(i=0;i<vsitecomm->right_import_nconstruct;i++)
		{
			ia = vsitecomm->right_import_construct[i];
			copy_rvec(f[ia],sendbuf[i]);
			clear_rvec(f[ia]); /* Zero it here after moving, just to simplify debug book-keeping... */
		}
		
		/* Pulse forces right */
		gmx_tx_rx_real(cr,GMX_RIGHT,(real *)sendbuf,3*vsitecomm->right_import_nconstruct,GMX_LEFT,(real *)recvbuf,3*vsitecomm->left_export_nconstruct);
		
		/* Copy from receive buffer to coordinate array */
		for(i=0;i<vsitecomm->left_export_nconstruct;i++)
		{
			ia = vsitecomm->left_export_construct[i];
			rvec_inc(f[ia],recvbuf[i]);
		}

		/* Prepare pulse left by copying to send buffer */
		for(i=0;i<vsitecomm->left_import_nconstruct;i++)
		{
			ia = vsitecomm->left_import_construct[i];
			copy_rvec(f[ia],sendbuf[i]);
			clear_rvec(f[ia]); /* Zero it here after moving, just to simplify debug book-keeping... */
		}
		
		/* Pulse coordinates left */
		gmx_tx_rx_real(cr,GMX_LEFT,(real *)sendbuf,3*vsitecomm->left_import_nconstruct,GMX_RIGHT,(real *)recvbuf,3*vsitecomm->right_export_nconstruct);
		
		/* Copy from receive buffer to coordinate array */
		for(i=0;i<vsitecomm->right_export_nconstruct;i++)
		{
			ia = vsitecomm->right_export_construct[i];
			rvec_inc(f[ia],recvbuf[i]);
		}
		
	/* All forces are now on the home processors */
}

	
static void
pd_clear_nonlocal_constructs(t_comm_vsites *vsitecomm, rvec f[])
{
	int i,ia;
	
	for(i=0;i<vsitecomm->left_import_nconstruct;i++)
	{
		ia = vsitecomm->left_import_construct[i];
		clear_rvec(f[ia]); 
	}
	for(i=0;i<vsitecomm->right_import_nconstruct;i++)
	{
		ia = vsitecomm->right_import_construct[i];
		clear_rvec(f[ia]); 
	}
}



static int pbc_rvec_sub(const t_pbc *pbc,const rvec xi,const rvec xj,rvec dx)
{
  if (pbc) {
    return pbc_dx_aiuc(pbc,xi,xj,dx);
  }
  else {
    rvec_sub(xi,xj,dx);
    return CENTRAL;
  }
}

/* Vsite construction routines */

static void constr_vsite2(rvec xi,rvec xj,rvec x,real a,t_pbc *pbc)
{
  real b;
  rvec dx;

  b=1.0-a;
  /* 1 flop */
  
  if (pbc) {
    pbc_dx_aiuc(pbc,xj,xi,dx);
    x[XX] = xi[XX] + a*dx[XX];
    x[YY] = xi[YY] + a*dx[YY];
    x[ZZ] = xi[ZZ] + a*dx[ZZ];
  } else {
    x[XX] = b*xi[XX] + a*xj[XX];
    x[YY] = b*xi[YY] + a*xj[YY];
    x[ZZ] = b*xi[ZZ] + a*xj[ZZ];
    /* 9 Flops */
  }
  
  /* TOTAL: 10 flops */
}

static void constr_vsite3(rvec xi,rvec xj,rvec xk,rvec x,real a,real b,
			  t_pbc *pbc)
{
  real c;
  rvec dxj,dxk;

  c=1.0-a-b;
  /* 2 flops */
  
  if (pbc) {
    pbc_dx_aiuc(pbc,xj,xi,dxj);
    pbc_dx_aiuc(pbc,xk,xi,dxk);
    x[XX] = xi[XX] + a*dxj[XX] + b*dxk[XX];
    x[YY] = xi[YY] + a*dxj[YY] + b*dxk[YY];
    x[ZZ] = xi[ZZ] + a*dxj[ZZ] + b*dxk[ZZ];
  } else {
    x[XX] = c*xi[XX] + a*xj[XX] + b*xk[XX];
    x[YY] = c*xi[YY] + a*xj[YY] + b*xk[YY];
    x[ZZ] = c*xi[ZZ] + a*xj[ZZ] + b*xk[ZZ];
  /* 15 Flops */
  }
  
  /* TOTAL: 17 flops */
}

static void constr_vsite3FD(rvec xi,rvec xj,rvec xk,rvec x,real a,real b,
			    t_pbc *pbc)
{
  rvec xij,xjk,temp;
  real c;
  
  pbc_rvec_sub(pbc,xj,xi,xij);
  pbc_rvec_sub(pbc,xk,xj,xjk);
  /* 6 flops */

  /* temp goes from i to a point on the line jk */  
  temp[XX] = xij[XX] + a*xjk[XX];
  temp[YY] = xij[YY] + a*xjk[YY];
  temp[ZZ] = xij[ZZ] + a*xjk[ZZ];
  /* 6 flops */
  
  c=b*gmx_invsqrt(iprod(temp,temp));
  /* 6 + 10 flops */
  
  x[XX] = xi[XX] + c*temp[XX];
  x[YY] = xi[YY] + c*temp[YY];
  x[ZZ] = xi[ZZ] + c*temp[ZZ];
  /* 6 Flops */
  
  /* TOTAL: 34 flops */
}

static void constr_vsite3FAD(rvec xi,rvec xj,rvec xk,rvec x,real a,real b, t_pbc *pbc)
{
  rvec xij,xjk,xp;
  real a1,b1,c1,invdij;
  
  pbc_rvec_sub(pbc,xj,xi,xij);
  pbc_rvec_sub(pbc,xk,xj,xjk);
  /* 6 flops */

  invdij = gmx_invsqrt(iprod(xij,xij));
  c1 = invdij * invdij * iprod(xij,xjk);
  xp[XX] = xjk[XX] - c1*xij[XX];
  xp[YY] = xjk[YY] - c1*xij[YY];
  xp[ZZ] = xjk[ZZ] - c1*xij[ZZ];
  a1 = a*invdij;
  b1 = b*gmx_invsqrt(iprod(xp,xp));
  /* 45 */
  
  x[XX] = xi[XX] + a1*xij[XX] + b1*xp[XX];
  x[YY] = xi[YY] + a1*xij[YY] + b1*xp[YY];
  x[ZZ] = xi[ZZ] + a1*xij[ZZ] + b1*xp[ZZ];
  /* 12 Flops */
  
  /* TOTAL: 63 flops */
}

static void constr_vsite3OUT(rvec xi,rvec xj,rvec xk,rvec x,
			     real a,real b,real c,t_pbc *pbc)
{
  rvec xij,xik,temp;
  
  pbc_rvec_sub(pbc,xj,xi,xij);
  pbc_rvec_sub(pbc,xk,xi,xik);
  cprod(xij,xik,temp);
  /* 15 Flops */
  
  x[XX] = xi[XX] + a*xij[XX] + b*xik[XX] + c*temp[XX];
  x[YY] = xi[YY] + a*xij[YY] + b*xik[YY] + c*temp[YY];
  x[ZZ] = xi[ZZ] + a*xij[ZZ] + b*xik[ZZ] + c*temp[ZZ];
  /* 18 Flops */
  
  /* TOTAL: 33 flops */
}

static void constr_vsite4FD(rvec xi,rvec xj,rvec xk,rvec xl,rvec x,
                            real a,real b,real c,t_pbc *pbc)
{
  rvec xij,xjk,xjl,temp;
  real d;
  
  pbc_rvec_sub(pbc,xj,xi,xij);
  pbc_rvec_sub(pbc,xk,xj,xjk);
  pbc_rvec_sub(pbc,xl,xj,xjl);
  /* 9 flops */

  /* temp goes from i to a point on the plane jkl */  
  temp[XX] = xij[XX] + a*xjk[XX] + b*xjl[XX];
  temp[YY] = xij[YY] + a*xjk[YY] + b*xjl[YY];
  temp[ZZ] = xij[ZZ] + a*xjk[ZZ] + b*xjl[ZZ];
  /* 12 flops */
  
  d=c*gmx_invsqrt(iprod(temp,temp));
  /* 6 + 10 flops */
  
  x[XX] = xi[XX] + d*temp[XX];
  x[YY] = xi[YY] + d*temp[YY];
  x[ZZ] = xi[ZZ] + d*temp[ZZ];
  /* 6 Flops */
  
  /* TOTAL: 43 flops */
}


static void constr_vsite4FDN(rvec xi,rvec xj,rvec xk,rvec xl,rvec x,
                             real a,real b,real c,t_pbc *pbc)
{
    rvec xij,xik,xil,ra,rb,rja,rjb,rm;
    real d;
    
    pbc_rvec_sub(pbc,xj,xi,xij);
    pbc_rvec_sub(pbc,xk,xi,xik);
    pbc_rvec_sub(pbc,xl,xi,xil);
    /* 9 flops */

    ra[XX] = a*xik[XX];
    ra[YY] = a*xik[YY];
    ra[ZZ] = a*xik[ZZ];
    
    rb[XX] = b*xil[XX];
    rb[YY] = b*xil[YY];
    rb[ZZ] = b*xil[ZZ];

    /* 6 flops */

    rvec_sub(ra,xij,rja);
    rvec_sub(rb,xij,rjb);
    /* 6 flops */
    
    cprod(rja,rjb,rm);
    /* 9 flops */
    
    d=c*gmx_invsqrt(norm2(rm));
    /* 5+5+1 flops */
    
    x[XX] = xi[XX] + d*rm[XX];
    x[YY] = xi[YY] + d*rm[YY];
    x[ZZ] = xi[ZZ] + d*rm[ZZ];
    /* 6 Flops */
    
    /* TOTAL: 47 flops */
}


static int constr_vsiten(t_iatom *ia, t_iparams ip[],
			 rvec *x, t_pbc *pbc)
{
  rvec xs,x1,dx;
  dvec dsum;
  int  n3,av,ai,i;
  real a;

  n3 = 3*ip[ia[0]].vsiten.n;
  av = ia[1];
  ai = ia[2];
  copy_rvec(x[ai],x1);
  clear_dvec(dsum);
  for(i=3; i<n3; i+=3) {
    ai = ia[i+2];
    a = ip[ia[i]].vsiten.a;
    if (pbc) {
      pbc_dx_aiuc(pbc,x[ai],x1,dx);
    } else {
      rvec_sub(x[ai],x1,dx);
    }
    dsum[XX] += a*dx[XX];
    dsum[YY] += a*dx[YY];
    dsum[ZZ] += a*dx[ZZ];
    /* 9 Flops */
  }

  x[av][XX] = x1[XX] + dsum[XX];
  x[av][YY] = x1[YY] + dsum[YY];
  x[av][ZZ] = x1[ZZ] + dsum[ZZ];

  return n3;
}


void construct_vsites(FILE *log,gmx_vsite_t *vsite,
		      rvec x[],t_nrnb *nrnb,
		      real dt,rvec *v,
		      t_iparams ip[],t_ilist ilist[],
		      int ePBC,gmx_bool bMolPBC,t_graph *graph,
		      t_commrec *cr,matrix box)
{
  rvec      xpbc,xv,vv,dx;
  real      a1,b1,c1,inv_dt;
  int       i,inc,ii,nra,nr,tp,ftype;
  t_iatom   avsite,ai,aj,ak,al,pbc_atom;
  t_iatom   *ia;
  t_pbc     pbc,*pbc_null,*pbc_null2;
  gmx_bool      bDomDec;
  int       *vsite_pbc,ishift;
  rvec      reftmp,vtmp,rtmp;
	
  bDomDec = cr && DOMAINDECOMP(cr);
		
  /* We only need to do pbc when we have inter-cg vsites */
  if (ePBC != epbcNONE && (bDomDec || bMolPBC) && vsite->n_intercg_vsite) {
    /* This is wasting some CPU time as we now do this multiple times
     * per MD step. But how often do we have vsites with full pbc?
     */
    pbc_null = set_pbc_dd(&pbc,ePBC,cr!=NULL ? cr->dd : NULL,FALSE,box);
  } else {
    pbc_null = NULL;
  }
		
  if (cr) {
    if (bDomDec) {
      dd_move_x_vsites(cr->dd,box,x);
    } else if (vsite->bPDvsitecomm) {
      /* I'm not sure whether the periodicity and shift are guaranteed
       * to be consistent between different nodes when running e.g. polymers
       * in parallel. In this special case we thus unshift/shift,
       * but only when necessary. This is to make sure the coordinates
       * we move don't end up a box away...
       */
		if (graph)
			unshift_self(graph,box,x);
		
		move_construct_x(vsite->vsitecomm,x,cr);
		
		if (graph)
			shift_self(graph,box,x);
    }
  }

  if (v) {
    inv_dt = 1.0/dt;
  } else {
    inv_dt = 1.0;
  }

  pbc_null2 = NULL;
  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (interaction_function[ftype].flags & IF_VSITE) {
      nra    = interaction_function[ftype].nratoms;
      nr     = ilist[ftype].nr;
      ia     = ilist[ftype].iatoms;

      if (pbc_null) {
	vsite_pbc = vsite->vsite_pbc_loc[ftype-F_VSITE2];
      } else {
	vsite_pbc = NULL;
      }
		
      for(i=0; (i<nr); ) {
	tp   = ia[0];
	/*
	if (ftype != idef->functype[tp]) 
	  gmx_incons("Function types for vsites wrong");
	*/
	
	/* The vsite and constructing atoms */
	avsite = ia[1];
	ai   = ia[2];
	aj   = ia[3];

	/* Constants for constructing vsites */
	a1   = ip[tp].vsite.a;
	/* Check what kind of pbc we need to use */
	if (vsite_pbc) {
	  pbc_atom = vsite_pbc[i/(1+nra)];
	  if (pbc_atom > -2) {
	    if (pbc_atom >= 0) {
	      /* We need to copy the coordinates here,
	       * single for single atom cg's pbc_atom is the vsite itself.
	       */
	      copy_rvec(x[pbc_atom],xpbc);
	    }
	    pbc_null2 = pbc_null;
	  } else {
	    pbc_null2 = NULL;
	  }
	} else {
	  pbc_atom = -2;
	}
	/* Copy the old position */
	copy_rvec(x[avsite],xv);

	/* Construct the vsite depending on type */
	inc = nra+1;
	switch (ftype) {
	case F_VSITE2:
	  constr_vsite2(x[ai],x[aj],x[avsite],a1,pbc_null2);
	  break;
	case F_VSITE3:
	  ak = ia[4];
	  b1 = ip[tp].vsite.b;
	  constr_vsite3(x[ai],x[aj],x[ak],x[avsite],a1,b1,pbc_null2);
	  break;
	case F_VSITE3FD:
	  ak = ia[4];
	  b1 = ip[tp].vsite.b;
	  constr_vsite3FD(x[ai],x[aj],x[ak],x[avsite],a1,b1,pbc_null2);
	  break;
	case F_VSITE3FAD:
	  ak = ia[4];
	  b1 = ip[tp].vsite.b;
	  constr_vsite3FAD(x[ai],x[aj],x[ak],x[avsite],a1,b1,pbc_null2);
	  break;
	case F_VSITE3OUT:
	  ak = ia[4];
	  b1 = ip[tp].vsite.b;
	  c1 = ip[tp].vsite.c;
	  constr_vsite3OUT(x[ai],x[aj],x[ak],x[avsite],a1,b1,c1,pbc_null2);
	  break;
	case F_VSITE4FD:
	  ak = ia[4];
	  al = ia[5];
	  b1 = ip[tp].vsite.b;
	  c1 = ip[tp].vsite.c;
	  constr_vsite4FD(x[ai],x[aj],x[ak],x[al],x[avsite],a1,b1,c1,
			  pbc_null2);
	  break;
	case F_VSITE4FDN:
	  ak = ia[4];
	  al = ia[5];
	  b1 = ip[tp].vsite.b;
	  c1 = ip[tp].vsite.c;
	  constr_vsite4FDN(x[ai],x[aj],x[ak],x[al],x[avsite],a1,b1,c1,
                        pbc_null2);
	  break;
	case F_VSITEN:
	  inc = constr_vsiten(ia,ip,x,pbc_null2);
	  break;
	default:
	  gmx_fatal(FARGS,"No such vsite type %d in %s, line %d",
		      ftype,__FILE__,__LINE__);
	}

	if (pbc_atom >= 0) {
	  /* Match the pbc of this vsite to the rest of its charge group */
	  ishift = pbc_dx_aiuc(pbc_null,x[avsite],xpbc,dx);
	  if (ishift != CENTRAL)
	    rvec_add(xpbc,dx,x[avsite]);
	}
	if (v) {
	  /* Calculate velocity of vsite... */
	  rvec_sub(x[avsite],xv,vv);
	  svmul(inv_dt,vv,v[avsite]);
	}

	/* Increment loop variables */
	i  += inc;
	ia += inc;
      }
    }
  }
}

static void spread_vsite2(t_iatom ia[],real a,
			    rvec x[],rvec f[],rvec fshift[],
			    t_pbc *pbc,t_graph *g)
{
  rvec    fi,fj,dx;
  t_iatom av,ai,aj;
  ivec    di;
  real    b;
  int     siv,sij;
  
  av = ia[1];
  ai = ia[2];
  aj = ia[3];
  
  svmul(1-a,f[av],fi);
  svmul(  a,f[av],fj);
  /* 7 flop */
  
  rvec_inc(f[ai],fi);
  rvec_inc(f[aj],fj);
  /* 6 Flops */

  if (g) {
    ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,av),di);
    siv = IVEC2IS(di);
    ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),di);
    sij = IVEC2IS(di);
  } else if (pbc) {
    siv = pbc_dx_aiuc(pbc,x[ai],x[av],dx);
    sij = pbc_dx_aiuc(pbc,x[ai],x[aj],dx);
  } else {
    siv = CENTRAL;
    sij = CENTRAL;
  }

  if (fshift && (siv != CENTRAL || sij != CENTRAL)) {
    rvec_inc(fshift[siv],f[av]);
    rvec_dec(fshift[CENTRAL],fi);
    rvec_dec(fshift[sij],fj);
  }

  /* TOTAL: 13 flops */
}

void construct_vsites_mtop(FILE *log,gmx_vsite_t *vsite,
			   gmx_mtop_t *mtop,rvec x[])
{
  int as,mb,mol;
  gmx_molblock_t *molb;
  gmx_moltype_t  *molt;

  as = 0;
  for(mb=0; mb<mtop->nmolblock; mb++) {
    molb = &mtop->molblock[mb];
    molt = &mtop->moltype[molb->type]; 	  
    for(mol=0; mol<molb->nmol; mol++) {
      construct_vsites(log,vsite,x+as,NULL,0.0,NULL,
		       mtop->ffparams.iparams,molt->ilist,
		       epbcNONE,TRUE,NULL,NULL,NULL);
      as += molt->atoms.nr;
    }
  }
}

static void spread_vsite3(t_iatom ia[],real a,real b,
			    rvec x[],rvec f[],rvec fshift[],
			    t_pbc *pbc,t_graph *g)
{
  rvec    fi,fj,fk,dx;
  atom_id av,ai,aj,ak;
  ivec    di;
  int     siv,sij,sik;

  av = ia[1];
  ai = ia[2];
  aj = ia[3];
  ak = ia[4];
  
  svmul(1-a-b,f[av],fi);
  svmul(    a,f[av],fj);
  svmul(    b,f[av],fk);
  /* 11 flops */

  rvec_inc(f[ai],fi);
  rvec_inc(f[aj],fj);
  rvec_inc(f[ak],fk);
  /* 9 Flops */
  
  if (g) {
    ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,ia[1]),di);
    siv = IVEC2IS(di);
    ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),di);
    sij = IVEC2IS(di);
    ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,ak),di);
    sik = IVEC2IS(di);
  } else if (pbc) {
    siv = pbc_dx_aiuc(pbc,x[ai],x[av],dx);
    sij = pbc_dx_aiuc(pbc,x[ai],x[aj],dx);
    sik = pbc_dx_aiuc(pbc,x[ai],x[ak],dx);
  } else {
    siv = CENTRAL;
    sij = CENTRAL;
    sik = CENTRAL;
  }

  if (fshift && (siv!=CENTRAL || sij!=CENTRAL || sik!=CENTRAL)) {
    rvec_inc(fshift[siv],f[av]);
    rvec_dec(fshift[CENTRAL],fi);
    rvec_dec(fshift[sij],fj);
    rvec_dec(fshift[sik],fk);
  }

  /* TOTAL: 20 flops */
}

static void spread_vsite3FD(t_iatom ia[],real a,real b,
			    rvec x[],rvec f[],rvec fshift[],
			    t_pbc *pbc,t_graph *g)
{
  real fx,fy,fz,c,invl,fproj,a1;
  rvec xvi,xij,xjk,xix,fv,temp;
  t_iatom av,ai,aj,ak;
  int     svi,sji,skj,d;
  ivec    di;

  av = ia[1];
  ai = ia[2];
  aj = ia[3];
  ak = ia[4];
  copy_rvec(f[av],fv);
  
  sji = pbc_rvec_sub(pbc,x[aj],x[ai],xij);
  skj = pbc_rvec_sub(pbc,x[ak],x[aj],xjk);
  /* 6 flops */

  /* xix goes from i to point x on the line jk */  
  xix[XX]=xij[XX]+a*xjk[XX];
  xix[YY]=xij[YY]+a*xjk[YY];
  xix[ZZ]=xij[ZZ]+a*xjk[ZZ];
  /* 6 flops */
  
  invl=gmx_invsqrt(iprod(xix,xix));
  c=b*invl;
  /* 4 + ?10? flops */
  
  fproj=iprod(xix,fv)*invl*invl; /* = (xix . f)/(xix . xix) */
  
  temp[XX]=c*(fv[XX]-fproj*xix[XX]);
  temp[YY]=c*(fv[YY]-fproj*xix[YY]);
  temp[ZZ]=c*(fv[ZZ]-fproj*xix[ZZ]);
  /* 16 */
  
  /* c is already calculated in constr_vsite3FD
     storing c somewhere will save 26 flops!     */
  
  a1=1-a;
  f[ai][XX] += fv[XX] - temp[XX];
  f[ai][YY] += fv[YY] - temp[YY];
  f[ai][ZZ] += fv[ZZ] - temp[ZZ];
  f[aj][XX] += a1*temp[XX];
  f[aj][YY] += a1*temp[YY];
  f[aj][ZZ] += a1*temp[ZZ];
  f[ak][XX] += a*temp[XX];
  f[ak][YY] += a*temp[YY];
  f[ak][ZZ] += a*temp[ZZ];
  /* 19 Flops */

  if (g) {
    ivec_sub(SHIFT_IVEC(g,ia[1]),SHIFT_IVEC(g,ai),di);
    svi = IVEC2IS(di);
    ivec_sub(SHIFT_IVEC(g,aj),SHIFT_IVEC(g,ai),di);
    sji = IVEC2IS(di);
    ivec_sub(SHIFT_IVEC(g,ak),SHIFT_IVEC(g,aj),di);
    skj = IVEC2IS(di);
  } else if (pbc) {
    svi = pbc_rvec_sub(pbc,x[av],x[ai],xvi);
  } else {
    svi = CENTRAL;
  }

  if (fshift && (svi!=CENTRAL || sji!=CENTRAL || skj!=CENTRAL)) {
    rvec_dec(fshift[svi],fv);
    fshift[CENTRAL][XX] += fv[XX] - (1 + a)*temp[XX];
    fshift[CENTRAL][YY] += fv[YY] - (1 + a)*temp[YY];
    fshift[CENTRAL][ZZ] += fv[ZZ] - (1 + a)*temp[ZZ];
    fshift[    sji][XX] += temp[XX];
    fshift[    sji][YY] += temp[YY];
    fshift[    sji][ZZ] += temp[ZZ];
    fshift[    skj][XX] += a*temp[XX];
    fshift[    skj][YY] += a*temp[YY];
    fshift[    skj][ZZ] += a*temp[ZZ];
  }

  /* TOTAL: 61 flops */
}

static void spread_vsite3FAD(t_iatom ia[],real a,real b,
			     rvec x[],rvec f[],rvec fshift[],
			     t_pbc *pbc,t_graph *g)
{
  rvec    xvi,xij,xjk,xperp,Fpij,Fppp,fv,f1,f2,f3;
  real    a1,b1,c1,c2,invdij,invdij2,invdp,fproj;
  t_iatom av,ai,aj,ak;
  int     svi,sji,skj,d;
  ivec    di;
  
  av = ia[1];
  ai = ia[2];
  aj = ia[3];
  ak = ia[4];
  copy_rvec(f[ia[1]],fv);

  sji = pbc_rvec_sub(pbc,x[aj],x[ai],xij);
  skj = pbc_rvec_sub(pbc,x[ak],x[aj],xjk);
  /* 6 flops */
  
  invdij = gmx_invsqrt(iprod(xij,xij));
  invdij2 = invdij * invdij;
  c1 = iprod(xij,xjk) * invdij2;
  xperp[XX] = xjk[XX] - c1*xij[XX];
  xperp[YY] = xjk[YY] - c1*xij[YY];
  xperp[ZZ] = xjk[ZZ] - c1*xij[ZZ];
  /* xperp in plane ijk, perp. to ij */
  invdp = gmx_invsqrt(iprod(xperp,xperp));
  a1 = a*invdij;
  b1 = b*invdp;
  /* 45 flops */
  
  /* a1, b1 and c1 are already calculated in constr_vsite3FAD
     storing them somewhere will save 45 flops!     */
  
  fproj=iprod(xij  ,fv)*invdij2;
  svmul(fproj,                      xij,  Fpij); /* proj. f on xij */
  svmul(iprod(xperp,fv)*invdp*invdp,xperp,Fppp); /* proj. f on xperp */
  svmul(b1*fproj,                   xperp,f3);
  /* 23 flops */
  
  rvec_sub(fv,Fpij,f1); /* f1 = f - Fpij */
  rvec_sub(f1,Fppp,f2); /* f2 = f - Fpij - Fppp */
  for (d=0; (d<DIM); d++) {
    f1[d]*=a1;
    f2[d]*=b1;
  }
  /* 12 flops */
  
  c2=1+c1;
  f[ai][XX] += fv[XX] - f1[XX] + c1*f2[XX] + f3[XX];
  f[ai][YY] += fv[YY] - f1[YY] + c1*f2[YY] + f3[YY];
  f[ai][ZZ] += fv[ZZ] - f1[ZZ] + c1*f2[ZZ] + f3[ZZ];
  f[aj][XX] +=          f1[XX] - c2*f2[XX] - f3[XX];
  f[aj][YY] +=          f1[YY] - c2*f2[YY] - f3[YY];
  f[aj][ZZ] +=          f1[ZZ] - c2*f2[ZZ] - f3[ZZ];
  f[ak][XX] +=                      f2[XX];
  f[ak][YY] +=                      f2[YY];
  f[ak][ZZ] +=                      f2[ZZ];
  /* 30 Flops */

  if (g) {
    ivec_sub(SHIFT_IVEC(g,ia[1]),SHIFT_IVEC(g,ai),di);
    svi = IVEC2IS(di);
    ivec_sub(SHIFT_IVEC(g,aj),SHIFT_IVEC(g,ai),di);
    sji = IVEC2IS(di);
    ivec_sub(SHIFT_IVEC(g,ak),SHIFT_IVEC(g,aj),di);
    skj = IVEC2IS(di);
  } else if (pbc) {
    svi = pbc_rvec_sub(pbc,x[av],x[ai],xvi);
  } else {
    svi = CENTRAL;
  }

  if (fshift && (svi!=CENTRAL || sji!=CENTRAL || skj!=CENTRAL)) {
    rvec_dec(fshift[svi],fv);
    fshift[CENTRAL][XX] += fv[XX] - f1[XX] - (1-c1)*f2[XX] + f3[XX];
    fshift[CENTRAL][YY] += fv[YY] - f1[YY] - (1-c1)*f2[YY] + f3[YY];
    fshift[CENTRAL][ZZ] += fv[ZZ] - f1[ZZ] - (1-c1)*f2[ZZ] + f3[ZZ];
    fshift[    sji][XX] +=          f1[XX] -    c1 *f2[XX] - f3[XX];
    fshift[    sji][YY] +=          f1[YY] -    c1 *f2[YY] - f3[YY];
    fshift[    sji][ZZ] +=          f1[ZZ] -    c1 *f2[ZZ] - f3[ZZ];
    fshift[    skj][XX] +=                          f2[XX];
    fshift[    skj][YY] +=                          f2[YY];
    fshift[    skj][ZZ] +=                          f2[ZZ];
  }
  
  /* TOTAL: 113 flops */
}

static void spread_vsite3OUT(t_iatom ia[],real a,real b,real c,
			     rvec x[],rvec f[],rvec fshift[],
			     t_pbc *pbc,t_graph *g)
{
  rvec    xvi,xij,xik,fv,fj,fk;
  real    cfx,cfy,cfz;
  atom_id av,ai,aj,ak;
  ivec    di;
  int     svi,sji,ski;
  
  av = ia[1];
  ai = ia[2];
  aj = ia[3];
  ak = ia[4];

  sji = pbc_rvec_sub(pbc,x[aj],x[ai],xij);
  ski = pbc_rvec_sub(pbc,x[ak],x[ai],xik);
  /* 6 Flops */
  
  copy_rvec(f[av],fv);

  cfx = c*fv[XX];
  cfy = c*fv[YY];
  cfz = c*fv[ZZ];
  /* 3 Flops */
  
  fj[XX] = a*fv[XX]     -  xik[ZZ]*cfy +  xik[YY]*cfz;
  fj[YY] =  xik[ZZ]*cfx + a*fv[YY]     -  xik[XX]*cfz;
  fj[ZZ] = -xik[YY]*cfx +  xik[XX]*cfy + a*fv[ZZ];
  
  fk[XX] = b*fv[XX]     +  xij[ZZ]*cfy -  xij[YY]*cfz;
  fk[YY] = -xij[ZZ]*cfx + b*fv[YY]     +  xij[XX]*cfz;
  fk[ZZ] =  xij[YY]*cfx -  xij[XX]*cfy + b*fv[ZZ];
  /* 30 Flops */
    
  f[ai][XX] += fv[XX] - fj[XX] - fk[XX];
  f[ai][YY] += fv[YY] - fj[YY] - fk[YY];
  f[ai][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ];
  rvec_inc(f[aj],fj);
  rvec_inc(f[ak],fk);
  /* 15 Flops */

  if (g) {
    ivec_sub(SHIFT_IVEC(g,ia[1]),SHIFT_IVEC(g,ai),di);
    svi = IVEC2IS(di);
    ivec_sub(SHIFT_IVEC(g,aj),SHIFT_IVEC(g,ai),di);
    sji = IVEC2IS(di);
    ivec_sub(SHIFT_IVEC(g,ak),SHIFT_IVEC(g,ai),di);
    ski = IVEC2IS(di);
  } else if (pbc) {
    svi = pbc_rvec_sub(pbc,x[av],x[ai],xvi);
  } else {
    svi = CENTRAL;
  }

  if (fshift && (svi!=CENTRAL || sji!=CENTRAL || ski!=CENTRAL)) {
    rvec_dec(fshift[svi],fv);
    fshift[CENTRAL][XX] += fv[XX] - fj[XX] - fk[XX];
    fshift[CENTRAL][YY] += fv[YY] - fj[YY] - fk[YY];
    fshift[CENTRAL][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ];
    rvec_inc(fshift[sji],fj);
    rvec_inc(fshift[ski],fk);
  }
  
  /* TOTAL: 54 flops */
}

static void spread_vsite4FD(t_iatom ia[],real a,real b,real c,
                            rvec x[],rvec f[],rvec fshift[],
                            t_pbc *pbc,t_graph *g)
{
  real    d,invl,fproj,a1;
  rvec    xvi,xij,xjk,xjl,xix,fv,temp;
  atom_id av,ai,aj,ak,al;
  ivec    di;
  int     svi,sji,skj,slj,m;

  av = ia[1];
  ai = ia[2];
  aj = ia[3];
  ak = ia[4];
  al = ia[5];
 
  sji = pbc_rvec_sub(pbc,x[aj],x[ai],xij);
  skj = pbc_rvec_sub(pbc,x[ak],x[aj],xjk);
  slj = pbc_rvec_sub(pbc,x[al],x[aj],xjl);
  /* 9 flops */
  
  /* xix goes from i to point x on the plane jkl */  
  for(m=0; m<DIM; m++)
    xix[m] = xij[m] + a*xjk[m] + b*xjl[m];
  /* 12 flops */
  
  invl=gmx_invsqrt(iprod(xix,xix));
  d=c*invl;
  /* 4 + ?10? flops */

  copy_rvec(f[av],fv);

  fproj=iprod(xix,fv)*invl*invl; /* = (xix . f)/(xix . xix) */

  for(m=0; m<DIM; m++)
    temp[m] = d*(fv[m] - fproj*xix[m]);
  /* 16 */
  
  /* c is already calculated in constr_vsite3FD
     storing c somewhere will save 35 flops!     */
  
  a1 = 1 - a - b;
  for(m=0; m<DIM; m++) {
    f[ai][m] += fv[m] - temp[m];
    f[aj][m] += a1*temp[m];
    f[ak][m] += a*temp[m];
    f[al][m] += b*temp[m];
  }
  /* 26 Flops */
  
  if (g) {
    ivec_sub(SHIFT_IVEC(g,ia[1]),SHIFT_IVEC(g,ai),di);
    svi = IVEC2IS(di);
    ivec_sub(SHIFT_IVEC(g,aj),SHIFT_IVEC(g,ai),di);
    sji = IVEC2IS(di);
    ivec_sub(SHIFT_IVEC(g,ak),SHIFT_IVEC(g,aj),di);
    skj = IVEC2IS(di);
    ivec_sub(SHIFT_IVEC(g,al),SHIFT_IVEC(g,aj),di);
    slj = IVEC2IS(di);
  } else if (pbc) {
    svi = pbc_rvec_sub(pbc,x[av],x[ai],xvi);
  } else {
    svi = CENTRAL;
  }

  if (fshift &&
      (svi!=CENTRAL || sji!=CENTRAL || skj!=CENTRAL || slj!=CENTRAL)) {
    rvec_dec(fshift[svi],fv);
    for(m=0; m<DIM; m++) {
      fshift[CENTRAL][m] += fv[m] - (1 + a + b)*temp[m];
      fshift[    sji][m] += temp[m];
      fshift[    skj][m] += a*temp[m];
      fshift[    slj][m] += b*temp[m];
    }
  }

  /* TOTAL: 77 flops */
}


static void spread_vsite4FDN(t_iatom ia[],real a,real b,real c,
                             rvec x[],rvec f[],rvec fshift[],
                             t_pbc *pbc,t_graph *g)
{
    rvec xvi,xij,xik,xil,ra,rb,rja,rjb,rab,rm,rt;
    rvec fv,fj,fk,fl;
    real invrm,denom;
    real cfx,cfy,cfz;
    ivec di;
    int  av,ai,aj,ak,al;
    int  svi,sij,sik,sil;

    /* DEBUG: check atom indices */
    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    al = ia[5];

    copy_rvec(f[av],fv);
    
    sij = pbc_rvec_sub(pbc,x[aj],x[ai],xij);
    sik = pbc_rvec_sub(pbc,x[ak],x[ai],xik);
    sil = pbc_rvec_sub(pbc,x[al],x[ai],xil);
    /* 9 flops */
    
    ra[XX] = a*xik[XX];
    ra[YY] = a*xik[YY];
    ra[ZZ] = a*xik[ZZ];
    
    rb[XX] = b*xil[XX];
    rb[YY] = b*xil[YY];
    rb[ZZ] = b*xil[ZZ];
    
    /* 6 flops */
    
    rvec_sub(ra,xij,rja);
    rvec_sub(rb,xij,rjb);
    rvec_sub(rb,ra,rab);
    /* 9 flops */
    
    cprod(rja,rjb,rm);
    /* 9 flops */

    invrm=gmx_invsqrt(norm2(rm));
    denom=invrm*invrm;
    /* 5+5+2 flops */
    
    cfx = c*invrm*fv[XX];
    cfy = c*invrm*fv[YY];
    cfz = c*invrm*fv[ZZ];
    /* 6 Flops */
    
    cprod(rm,rab,rt);
    /* 9 flops */

    rt[XX] *= denom;
    rt[YY] *= denom;
    rt[ZZ] *= denom;
    /* 3flops */
    
    fj[XX] = (        -rm[XX]*rt[XX]) * cfx + ( rab[ZZ]-rm[YY]*rt[XX]) * cfy + (-rab[YY]-rm[ZZ]*rt[XX]) * cfz;
    fj[YY] = (-rab[ZZ]-rm[XX]*rt[YY]) * cfx + (        -rm[YY]*rt[YY]) * cfy + ( rab[XX]-rm[ZZ]*rt[YY]) * cfz;
    fj[ZZ] = ( rab[YY]-rm[XX]*rt[ZZ]) * cfx + (-rab[XX]-rm[YY]*rt[ZZ]) * cfy + (        -rm[ZZ]*rt[ZZ]) * cfz;
    /* 30 flops */
        
    cprod(rjb,rm,rt);
    /* 9 flops */

    rt[XX] *= denom*a;
    rt[YY] *= denom*a;
    rt[ZZ] *= denom*a;
    /* 3flops */
    
    fk[XX] = (          -rm[XX]*rt[XX]) * cfx + (-a*rjb[ZZ]-rm[YY]*rt[XX]) * cfy + ( a*rjb[YY]-rm[ZZ]*rt[XX]) * cfz;
    fk[YY] = ( a*rjb[ZZ]-rm[XX]*rt[YY]) * cfx + (          -rm[YY]*rt[YY]) * cfy + (-a*rjb[XX]-rm[ZZ]*rt[YY]) * cfz;
    fk[ZZ] = (-a*rjb[YY]-rm[XX]*rt[ZZ]) * cfx + ( a*rjb[XX]-rm[YY]*rt[ZZ]) * cfy + (          -rm[ZZ]*rt[ZZ]) * cfz;
    /* 36 flops */
    
    cprod(rm,rja,rt);
    /* 9 flops */
    
    rt[XX] *= denom*b;
    rt[YY] *= denom*b;
    rt[ZZ] *= denom*b;
    /* 3flops */
    
    fl[XX] = (          -rm[XX]*rt[XX]) * cfx + ( b*rja[ZZ]-rm[YY]*rt[XX]) * cfy + (-b*rja[YY]-rm[ZZ]*rt[XX]) * cfz;
    fl[YY] = (-b*rja[ZZ]-rm[XX]*rt[YY]) * cfx + (          -rm[YY]*rt[YY]) * cfy + ( b*rja[XX]-rm[ZZ]*rt[YY]) * cfz;
    fl[ZZ] = ( b*rja[YY]-rm[XX]*rt[ZZ]) * cfx + (-b*rja[XX]-rm[YY]*rt[ZZ]) * cfy + (          -rm[ZZ]*rt[ZZ]) * cfz;
    /* 36 flops */

    f[ai][XX] += fv[XX] - fj[XX] - fk[XX] - fl[XX];
    f[ai][YY] += fv[YY] - fj[YY] - fk[YY] - fl[YY];
    f[ai][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ] - fl[ZZ];
    rvec_inc(f[aj],fj);
    rvec_inc(f[ak],fk);
    rvec_inc(f[al],fl);
    /* 21 flops */

    if (g) {
        ivec_sub(SHIFT_IVEC(g,av),SHIFT_IVEC(g,ai),di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g,aj),SHIFT_IVEC(g,ai),di);
        sij = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g,ak),SHIFT_IVEC(g,ai),di);
        sik = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g,al),SHIFT_IVEC(g,ai),di);
        sil = IVEC2IS(di);
    } else if (pbc) {
        svi = pbc_rvec_sub(pbc,x[av],x[ai],xvi);
    } else {
        svi = CENTRAL;
    }
    
    if (fshift && (svi!=CENTRAL || sij!=CENTRAL || sik!=CENTRAL || sil!=CENTRAL)) {
        rvec_dec(fshift[svi],fv);
        fshift[CENTRAL][XX] += fv[XX] - fj[XX] - fk[XX] - fl[XX];
        fshift[CENTRAL][YY] += fv[YY] - fj[YY] - fk[YY] - fl[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ] - fl[ZZ];
        rvec_inc(fshift[sij],fj);
        rvec_inc(fshift[sik],fk);
        rvec_inc(fshift[sil],fl);
    }
    
    /* Total: 207 flops (Yuck!) */
}


static int spread_vsiten(t_iatom ia[],t_iparams ip[],
			 rvec x[],rvec f[],rvec fshift[],
			 t_pbc *pbc,t_graph *g)
{
  rvec xv,dx,fi;
  int  n3,av,i,ai;
  real a;
  ivec di;
  int  siv;

  n3 = 3*ip[ia[0]].vsiten.n;
  av = ia[1];
  copy_rvec(x[av],xv);
  
  for(i=0; i<n3; i+=3) {
    ai = ia[i+2];
    if (g) {
      ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,av),di);
      siv = IVEC2IS(di);
    } else if (pbc) {
      siv = pbc_dx_aiuc(pbc,x[ai],xv,dx);
    } else {
      siv = CENTRAL;
    }
    a = ip[ia[i]].vsiten.a;
    svmul(a,f[av],fi);
    rvec_inc(f[ai],fi);
    if (fshift && siv != CENTRAL) {
      rvec_inc(fshift[siv],fi);
      rvec_dec(fshift[CENTRAL],fi);
    }
    /* 6 Flops */
  }

  return n3;
}


void spread_vsite_f(FILE *log,gmx_vsite_t *vsite,
		    rvec x[],rvec f[],rvec *fshift,
		    t_nrnb *nrnb,t_idef *idef,
		    int ePBC,gmx_bool bMolPBC,t_graph *g,matrix box,
		    t_commrec *cr)
{
  real      a1,b1,c1;
  int       i,inc,m,nra,nr,tp,ftype;
  int       nd2,nd3,nd3FD,nd3FAD,nd3OUT,nd4FD,nd4FDN,ndN;
  t_iatom   *ia;
  t_iparams *ip;
  t_pbc     pbc,*pbc_null,*pbc_null2;
  int       *vsite_pbc;

  /* We only need to do pbc when we have inter-cg vsites */
  if ((DOMAINDECOMP(cr) || bMolPBC) && vsite->n_intercg_vsite) {
    /* This is wasting some CPU time as we now do this multiple times
     * per MD step. But how often do we have vsites with full pbc?
     */
    pbc_null = set_pbc_dd(&pbc,ePBC,cr->dd,FALSE,box);
  } else {
    pbc_null = NULL;
  }
  
  if (DOMAINDECOMP(cr)) 
  {
    dd_clear_f_vsites(cr->dd,f);
  } 
  else if (PARTDECOMP(cr) && vsite->vsitecomm != NULL)
  {
    pd_clear_nonlocal_constructs(vsite->vsitecomm,f);
  }
	
  ip     = idef->iparams;

  nd2        = 0;
  nd3        = 0;
  nd3FD      = 0;
  nd3FAD     = 0;
  nd3OUT     = 0;
  nd4FD      = 0;
  nd4FDN     = 0;
  ndN        = 0;
   
  /* this loop goes backwards to be able to build *
   * higher type vsites from lower types         */
  pbc_null2 = NULL;
  for(ftype=F_NRE-1; (ftype>=0); ftype--) {
    if (interaction_function[ftype].flags & IF_VSITE) {
      nra    = interaction_function[ftype].nratoms;
      nr     = idef->il[ftype].nr;
      ia     = idef->il[ftype].iatoms;

      if (pbc_null) {
	vsite_pbc = vsite->vsite_pbc_loc[ftype-F_VSITE2];
      } else {
	vsite_pbc = NULL;
      }
     
      for(i=0; (i<nr); ) {
	/* Check if we need to apply pbc for this vsite */
	if (vsite_pbc) {
	  if (vsite_pbc[i/(1+nra)] > -2)
	    pbc_null2 = pbc_null;
	  else
	    pbc_null2 = NULL;
	}

	tp   = ia[0];
	if (ftype != idef->functype[tp])
	  gmx_incons("Functiontypes for vsites wrong");

	/* Constants for constructing */
	a1   = ip[tp].vsite.a; 
	/* Construct the vsite depending on type */
	inc = nra+1;
	switch (ftype) {
	case F_VSITE2:
	  spread_vsite2(ia,a1,x,f,fshift,pbc_null2,g);
	  nd2++;
	  break;
	case F_VSITE3:
	  b1 = ip[tp].vsite.b;
	  spread_vsite3(ia,a1,b1,x,f,fshift,pbc_null2,g);
	  nd3++;
	  break;
	case F_VSITE3FD:
	  b1 = ip[tp].vsite.b;
	  spread_vsite3FD(ia,a1,b1,x,f,fshift,pbc_null2,g);
	  nd3FD++;
	  break;
	case F_VSITE3FAD:
	  b1 = ip[tp].vsite.b;
	  spread_vsite3FAD(ia,a1,b1,x,f,fshift,pbc_null2,g);
	  nd3FAD++;
	  break;
	case F_VSITE3OUT:
	  b1 = ip[tp].vsite.b;
	  c1 = ip[tp].vsite.c;
	  spread_vsite3OUT(ia,a1,b1,c1,x,f,fshift,pbc_null2,g);
	  nd3OUT++;
	  break;
	case F_VSITE4FD:
	  b1 = ip[tp].vsite.b;
	  c1 = ip[tp].vsite.c;
	  spread_vsite4FD(ia,a1,b1,c1,x,f,fshift,pbc_null2,g);
	  nd4FD++;
	  break;
	case F_VSITE4FDN:
	  b1 = ip[tp].vsite.b;
	  c1 = ip[tp].vsite.c;
	  spread_vsite4FDN(ia,a1,b1,c1,x,f,fshift,pbc_null2,g);
	  nd4FDN++;
	  break;
	case F_VSITEN:
	  inc = spread_vsiten(ia,ip,x,f,fshift,pbc_null2,g);
	  ndN += inc;
	  break;
	default:
	  gmx_fatal(FARGS,"No such vsite type %d in %s, line %d",
		      ftype,__FILE__,__LINE__);
	}
	clear_rvec(f[ia[1]]);
	
	/* Increment loop variables */
	i  += inc;
	ia += inc;
      }
    }
  }
	
  inc_nrnb(nrnb,eNR_VSITE2,   nd2     );
  inc_nrnb(nrnb,eNR_VSITE3,   nd3     );
  inc_nrnb(nrnb,eNR_VSITE3FD, nd3FD   );
  inc_nrnb(nrnb,eNR_VSITE3FAD,nd3FAD  );
  inc_nrnb(nrnb,eNR_VSITE3OUT,nd3OUT  );
  inc_nrnb(nrnb,eNR_VSITE4FD, nd4FD   );
  inc_nrnb(nrnb,eNR_VSITE4FDN,nd4FDN  );
  inc_nrnb(nrnb,eNR_VSITEN,   ndN     );

  if (DOMAINDECOMP(cr)) {
    dd_move_f_vsites(cr->dd,f,fshift);
  } else if (vsite->bPDvsitecomm) {
    /* We only move forces here, and they are independent of shifts */
    move_construct_f(vsite->vsitecomm,f,cr);
  }
}

static int *atom2cg(t_block *cgs)
{
  int *a2cg,cg,i;
  
  snew(a2cg,cgs->index[cgs->nr]);
  for(cg=0; cg<cgs->nr; cg++) {
    for(i=cgs->index[cg]; i<cgs->index[cg+1]; i++)
      a2cg[i] = cg;
  }
  
  return a2cg;
}

static int count_intercg_vsite(gmx_mtop_t *mtop)
{
  int  mb,mt,ftype,nral,i,cg,a;
  gmx_molblock_t *molb;
  gmx_moltype_t *molt;
  int  *a2cg;
  t_ilist *il;
  t_iatom *ia;
  int  n_intercg_vsite;

  n_intercg_vsite = 0;
  for(mb=0; mb<mtop->nmolblock; mb++) {
    molb = &mtop->molblock[mb];
    molt = &mtop->moltype[molb->type];
    a2cg = atom2cg(&molt->cgs);
    for(ftype=0; ftype<F_NRE; ftype++) {
      if (interaction_function[ftype].flags & IF_VSITE) {
	nral = NRAL(ftype);
	il = &molt->ilist[ftype];
	ia  = il->iatoms;
	for(i=0; i<il->nr; i+=1+nral) {
	  cg = a2cg[ia[1+i]];
	  for(a=1; a<nral; a++) {
	    if (a2cg[ia[1+a]] != cg) {
	      n_intercg_vsite += molb->nmol;
	      break;
	    }
	  }
	}
      }
    }
    sfree(a2cg);
  }

  return n_intercg_vsite;
}

static int **get_vsite_pbc(t_iparams *iparams,t_ilist *ilist,
			   t_atom *atom,t_mdatoms *md,
			   t_block *cgs,int *a2cg)
{
  int  ftype,nral,i,j,vsi,vsite,cg_v,cg_c,a,nc3=0;
  t_ilist *il;
  t_iatom *ia;
  int  **vsite_pbc,*vsite_pbc_f;
  char *pbc_set;
  gmx_bool bViteOnlyCG_and_FirstAtom;

  /* Make an array that tells if the pbc of an atom is set */
  snew(pbc_set,cgs->index[cgs->nr]);
  /* PBC is set for all non vsites */
  for(a=0; a<cgs->index[cgs->nr]; a++) {
    if ((atom && atom[a].ptype != eptVSite) ||
	(md   && md->ptype[a]  != eptVSite)) {
      pbc_set[a] = 1;
    }
  }

  snew(vsite_pbc,F_VSITEN-F_VSITE2+1);
  
  for(ftype=0; ftype<F_NRE; ftype++) {
    if (interaction_function[ftype].flags & IF_VSITE) {
      nral = NRAL(ftype);
      il = &ilist[ftype];
      ia  = il->iatoms;

      snew(vsite_pbc[ftype-F_VSITE2],il->nr/(1+nral));
      vsite_pbc_f = vsite_pbc[ftype-F_VSITE2];

      i = 0;
      while (i < il->nr) {
	vsi = i/(1+nral);
	vsite = ia[i+1];
	cg_v = a2cg[vsite];
	/* A value of -2 signals that this vsite and its contructing
	 * atoms are all within the same cg, so no pbc is required.
	 */
	vsite_pbc_f[vsi] = -2;
	/* Check if constructing atoms are outside the vsite's cg */
	nc3 = 0;
	if (ftype == F_VSITEN) {
	  nc3 = 3*iparams[ia[i]].vsiten.n;
	  for(j=0; j<nc3; j+=3) {
	    if (a2cg[ia[i+j+2]] != cg_v)
	      vsite_pbc_f[vsi] = -1;
	  }
	} else {
	  for(a=1; a<nral; a++) {
	    if (a2cg[ia[i+1+a]] != cg_v)
	      vsite_pbc_f[vsi] = -1;
	  }
	}
	if (vsite_pbc_f[vsi] == -1) {
	  /* Check if this is the first processed atom of a vsite only cg */
	  bViteOnlyCG_and_FirstAtom = TRUE;
	  for(a=cgs->index[cg_v]; a<cgs->index[cg_v+1]; a++) {
	    /* Non-vsites already have pbc set, so simply check for pbc_set */
	    if (pbc_set[a]) {
	      bViteOnlyCG_and_FirstAtom = FALSE;
	      break;
	    }
	  }
	  if (bViteOnlyCG_and_FirstAtom) {
	    /* First processed atom of a vsite only charge group.
	     * The pbc of the input coordinates to construct_vsites
	     * should be preserved.
	     */
	    vsite_pbc_f[vsi] = vsite;
	  } else if (cg_v != a2cg[ia[1+i+1]]) {
	    /* This vsite has a different charge group index
	     * than it's first constructing atom
	     * and the charge group has more than one atom,
	     * search for the first normal particle
	     * or vsite that already had its pbc defined.
	     * If nothing is found, use full pbc for this vsite.
	     */
	    for(a=cgs->index[cg_v]; a<cgs->index[cg_v+1]; a++) {
	      if (a != vsite && pbc_set[a]) {
		vsite_pbc_f[vsi] = a;
		if (gmx_debug_at)
		  fprintf(debug,"vsite %d match pbc with atom %d\n",
			  vsite+1,a+1);
		break;
	      }
	    }
	    if (gmx_debug_at)
	      fprintf(debug,"vsite atom %d  cg %d - %d pbc atom %d\n",
		      vsite+1,cgs->index[cg_v]+1,cgs->index[cg_v+1],
		      vsite_pbc_f[vsi]+1);
	  }
	}
	if (ftype == F_VSITEN) {
	  /* The other entries in vsite_pbc_f are not used for center vsites */
	  i += nc3;
	} else {
	  i += 1+nral;
	}
	
	/* This vsite now has its pbc defined */
	pbc_set[vsite] = 1;
      }
    }
  }

  sfree(pbc_set);

  return vsite_pbc;
}

gmx_vsite_t *init_vsite(gmx_mtop_t *mtop,t_commrec *cr)
{
  int nvsite,i;
  int *a2cg,cg;
  gmx_vsite_t *vsite;
  int mt;
  gmx_moltype_t *molt;
  
  /* check if there are vsites */
  nvsite = 0;
  for(i=0; i<F_NRE; i++) {
    if (interaction_function[i].flags & IF_VSITE) {
      nvsite += gmx_mtop_ftype_count(mtop,i);
    }
  }

  if (nvsite == 0) {
    return NULL;
  }

  snew(vsite,1);

  vsite->n_intercg_vsite = count_intercg_vsite(mtop);

  if (vsite->n_intercg_vsite > 0 && DOMAINDECOMP(cr)) {
    vsite->nvsite_pbc_molt = mtop->nmoltype;
    snew(vsite->vsite_pbc_molt,vsite->nvsite_pbc_molt);
    for(mt=0; mt<mtop->nmoltype; mt++) {
      molt = &mtop->moltype[mt];
      /* Make an atom to charge group index */
      a2cg = atom2cg(&molt->cgs);
      vsite->vsite_pbc_molt[mt] = get_vsite_pbc(mtop->ffparams.iparams,
						molt->ilist,
						molt->atoms.atom,NULL,
						&molt->cgs,a2cg);
      sfree(a2cg);
    }

    snew(vsite->vsite_pbc_loc_nalloc,F_VSITEN-F_VSITE2+1);
    snew(vsite->vsite_pbc_loc       ,F_VSITEN-F_VSITE2+1);
  }


  return vsite;
}

void set_vsite_top(gmx_vsite_t *vsite,gmx_localtop_t *top,t_mdatoms *md,
		   t_commrec *cr)
{
  int *a2cg;

  /* Make an atom to charge group index */
  a2cg = atom2cg(&top->cgs);

  if (vsite->n_intercg_vsite > 0) {
    vsite->vsite_pbc_loc = get_vsite_pbc(top->idef.iparams,
					 top->idef.il,NULL,md,
					 &top->cgs,a2cg);

    if (PARTDECOMP(cr)) {
      snew(vsite->vsitecomm,1);
		vsite->bPDvsitecomm =
		setup_parallel_vsites(&(top->idef),cr,vsite->vsitecomm);
    }
  }

  sfree(a2cg);
}
