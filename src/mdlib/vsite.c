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

/* Communication buffers */

static rvec *prevbuf=NULL,*nextbuf=NULL;

/* Routines to send/recieve coordinates and force
 * of constructing atoms and vsites. This is necessary
 * when vsite constructs cross node borders (unavoidable
 * for e.g. polymers using anisotropic united atoms).
 */
/* Communication routines for vsites. The coordinates and
 * forces are only move on a need-to-know basis, usually only
 * 2-3 atoms per processor. To achieve this small amount of
 * communication, and to limit it to nearest neighbour messages,
 * we demand that vsites are not spread over nonadjacent nodes.
 * Thus, keep your vsites close to your constructing atoms.
 * (mdrun & grompp will report an error otherwise)
 */ 


static void move_construct_x(t_comm_vsites *vsitecomm, rvec x[], t_commrec *cr)
{
  static bool bFirst=TRUE;
  int i;
   
  if (bFirst) {
    /* Make the larger than necessary to avoid cache sharing */
    snew(nextbuf,2*(vsitecomm->nnextvsite+vsitecomm->nnextconstr)+100);
    snew(prevbuf,2*(vsitecomm->nprevvsite+vsitecomm->nprevconstr)+100);
    bFirst=FALSE;
  }
   
  /* package coords to send left. Vsite coords are needed to create v */
  for(i=0;i<vsitecomm->nprevconstr;i++)
    copy_rvec(x[vsitecomm->idxprevconstr[i]],prevbuf[i]);
  for(i=0;i<vsitecomm->nprevvsite;i++)
    copy_rvec(x[vsitecomm->idxprevvsite[i]],prevbuf[vsitecomm->nprevconstr+i]);
  
  /* send them off, and recieve from the right */
  if(vsitecomm->nprevconstr>0 || vsitecomm->nprevvsite>0)
    gmx_tx(cr->left,prevbuf,
	   sizeof(rvec)*(vsitecomm->nprevconstr+vsitecomm->nprevvsite));
  
  if(vsitecomm->nnextconstr>0 || vsitecomm->nnextvsite>0)
    gmx_rx(cr->right,nextbuf,
	   sizeof(rvec)*(vsitecomm->nnextconstr+vsitecomm->nnextvsite));
  
  if(vsitecomm->nprevconstr>0 || vsitecomm->nprevvsite>0)
    gmx_tx_wait(cr->left);
  
  if(vsitecomm->nnextconstr>0 || vsitecomm->nnextvsite>0)
    gmx_rx_wait(cr->right);
  
  /* Put them where they belong */
  for(i=0;i<vsitecomm->nnextconstr;i++)
    copy_rvec(nextbuf[i],x[vsitecomm->idxnextconstr[i]]);
  for(i=0;i<vsitecomm->nnextvsite;i++)
    copy_rvec(nextbuf[vsitecomm->nnextconstr+i],
	      x[vsitecomm->idxnextvsite[i]]);

  /* Now we are ready to do the constructing business ! */
}


static void move_vsite_xv(t_comm_vsites *vsitecomm, rvec x[], rvec v[],t_commrec *cr)
{
  int i;
  int sendsize,recvsize;

  sendsize=sizeof(rvec)*vsitecomm->nnextvsite;
  recvsize=sizeof(rvec)*vsitecomm->nprevvsite;

  if(v!=NULL) {
    sendsize=sendsize*2;
    recvsize=recvsize*2;
  }

  /* Package nonlocal constructed vsites */
  for(i=0;i<vsitecomm->nnextvsite;i++)
    copy_rvec(x[vsitecomm->idxnextvsite[i]],nextbuf[i]);

  if(v!=NULL)
    for(i=0;i<vsitecomm->nnextvsite;i++)
      copy_rvec(v[vsitecomm->idxnextvsite[i]],nextbuf[vsitecomm->nnextvsite+i]);
  
  /* send them off, and recieve from the right */
  if(vsitecomm->nnextvsite>0)
    gmx_tx(cr->right,nextbuf,sendsize);
  
  if(vsitecomm->nprevvsite>0)
    gmx_rx(cr->left,prevbuf,recvsize);
  
  if(vsitecomm->nnextvsite>0)
    gmx_tx_wait(cr->right);
  
  if(vsitecomm->nprevvsite>0)
    gmx_rx_wait(cr->left);
  
  /* Put them where they belong */
  for(i=0;i<vsitecomm->nprevvsite;i++)
    copy_rvec(prevbuf[i],x[vsitecomm->idxprevvsite[i]]);

  if(v!=NULL)
    for(i=0;i<vsitecomm->nprevvsite;i++)
      copy_rvec(prevbuf[vsitecomm->nprevvsite+i],v[vsitecomm->idxprevvsite[i]]);

  /* All coordinates are in place on the respective home node now */
}

static void move_vsite_f(t_comm_vsites *vsitecomm, rvec f[], t_commrec *cr)
{
  int i;

  /* package vsite particle forces to send left */
  for(i=0;i<vsitecomm->nprevvsite;i++)
    copy_rvec(f[vsitecomm->idxprevvsite[i]],prevbuf[i]);

  /* off they go! - but only if there is something to send! */
  if(vsitecomm->nprevvsite>0)
    gmx_tx(cr->left,prevbuf,sizeof(rvec)*vsitecomm->nprevvsite);

  /* Get our share from the right, if there is anything to have */
  if(vsitecomm->nnextvsite>0)
    gmx_rx(cr->right,nextbuf,sizeof(rvec)*vsitecomm->nnextvsite);
  
  if(vsitecomm->nprevvsite>0)
    gmx_tx_wait(cr->left);
  
  if(vsitecomm->nnextvsite>0)
    gmx_rx_wait(cr->right);

  /* Put them where they belong */
  for(i=0;i<vsitecomm->nnextvsite;i++)
    copy_rvec(nextbuf[i],f[vsitecomm->idxnextvsite[i]]);
  
  /* Zero forces on nonlocal constructing atoms.
   * This is necessary since vsite force spreading is done
   * after the normal force addition, and we don't want
   * to include them twice.
   * (They have already been added on the home node).
   */
  for(i=0;i<vsitecomm->nnextconstr;i++)
    clear_rvec(f[vsitecomm->idxnextconstr[i]]);
}

static void move_construct_f(t_comm_vsites *vsitecomm, rvec f[], t_commrec *cr)
{
  int i;

  /* Spread forces to nonlocal constructing atoms.
   */
  /* package forces to send right */
  for(i=0;i<vsitecomm->nnextconstr;i++)
    copy_rvec(f[vsitecomm->idxnextconstr[i]],nextbuf[i]);
  
  /* send them off, and recieve from the right */
  if(vsitecomm->nnextconstr>0)
    gmx_tx(cr->right,nextbuf,sizeof(rvec)*vsitecomm->nnextconstr);
  
  if(vsitecomm->nprevconstr>0)
    gmx_rx(cr->left,prevbuf,sizeof(rvec)*vsitecomm->nprevconstr);
  
  if(vsitecomm->nnextconstr>0)
    gmx_tx_wait(cr->right);

  if(vsitecomm->nprevconstr>0)
    gmx_rx_wait(cr->left);
  
  /* Add them where they belong */
  for(i=0;i<vsitecomm->nprevconstr;i++)
    rvec_inc(f[vsitecomm->idxprevconstr[i]],prevbuf[i]);
  
  /* Zero nonlocal vsites */
  for(i=0;i<vsitecomm->nprevvsite;i++)
    clear_rvec(f[vsitecomm->idxprevvsite[i]]);
  
  /* All forces are on the home processor now */  
}


/* Vsite construction routines */

static void constr_vsite2(rvec xi,rvec xj,rvec x,real a)
{
  real b;
  
  b=1.0-a;
  /* 1 flop */
  
  x[XX]=b*xi[XX]+a*xj[XX];
  x[YY]=b*xi[YY]+a*xj[YY];
  x[ZZ]=b*xi[ZZ]+a*xj[ZZ];
  /* 9 Flops */
  
  /* TOTAL: 10 flops */
}

static void constr_vsite3(rvec xi,rvec xj,rvec xk,rvec x,real a,real b)
{
  real c;
  
  c=1.0-a-b;
  /* 2 flops */
  
  x[XX] = c*xi[XX] + a*xj[XX] + b*xk[XX];
  x[YY] = c*xi[YY] + a*xj[YY] + b*xk[YY];
  x[ZZ] = c*xi[ZZ] + a*xj[ZZ] + b*xk[ZZ];
  /* 15 Flops */
  
  /* TOTAL: 17 flops */
}

static void constr_vsite3FD(rvec xi,rvec xj,rvec xk,rvec x,real a,real b)
{
  rvec xij,xjk,temp;
  real c;
  
  rvec_sub(xj,xi,xij);
  rvec_sub(xk,xj,xjk);
  /* 6 flops */

  /* temp goes from i to a point on the line jk */  
  temp[XX] = xij[XX] + a*xjk[XX];
  temp[YY] = xij[YY] + a*xjk[YY];
  temp[ZZ] = xij[ZZ] + a*xjk[ZZ];
  /* 6 flops */
  
  c=b*invsqrt(iprod(temp,temp));
  /* 6 + 10 flops */
  
  x[XX] = xi[XX] + c*temp[XX];
  x[YY] = xi[YY] + c*temp[YY];
  x[ZZ] = xi[ZZ] + c*temp[ZZ];
  /* 6 Flops */
  
  /* TOTAL: 34 flops */
}

static void constr_vsite3FAD(rvec xi,rvec xj,rvec xk,rvec x,real a,real b)
{
  rvec xij,xjk,xp;
  real a1,b1,c1,invdij;
  
  rvec_sub(xj,xi,xij);
  rvec_sub(xk,xj,xjk);
  /* 6 flops */

  invdij = invsqrt(iprod(xij,xij));
  c1 = invdij * invdij * iprod(xij,xjk);
  xp[XX] = xjk[XX] - c1*xij[XX];
  xp[YY] = xjk[YY] - c1*xij[YY];
  xp[ZZ] = xjk[ZZ] - c1*xij[ZZ];
  a1 = a*invdij;
  b1 = b*invsqrt(iprod(xp,xp));
  /* 45 */
  
  x[XX] = xi[XX] + a1*xij[XX] + b1*xp[XX];
  x[YY] = xi[YY] + a1*xij[YY] + b1*xp[YY];
  x[ZZ] = xi[ZZ] + a1*xij[ZZ] + b1*xp[ZZ];
  /* 12 Flops */
  
  /* TOTAL: 63 flops */
}

static void constr_vsite3OUT(rvec xi,rvec xj,rvec xk,rvec x,real a,real b,real c)
{
  rvec xij,xik,temp;
  
  rvec_sub(xj,xi,xij);
  rvec_sub(xk,xi,xik);
  oprod(xij,xik,temp);
  /* 15 Flops */
  
  x[XX] = xi[XX] + a*xij[XX] + b*xik[XX] + c*temp[XX];
  x[YY] = xi[YY] + a*xij[YY] + b*xik[YY] + c*temp[YY];
  x[ZZ] = xi[ZZ] + a*xij[ZZ] + b*xik[ZZ] + c*temp[ZZ];
  /* 18 Flops */
  
  /* TOTAL: 33 flops */
}

static void constr_vsite4FD(rvec xi,rvec xj,rvec xk,rvec xl,rvec x,
			  real a,real b,real c)
{
  rvec xij,xjk,xjl,temp;
  real d;
  
  rvec_sub(xj,xi,xij);
  rvec_sub(xk,xj,xjk);
  rvec_sub(xl,xj,xjl);
  /* 9 flops */

  /* temp goes from i to a point on the plane jkl */  
  temp[XX] = xij[XX] + a*xjk[XX] + b*xjl[XX];
  temp[YY] = xij[YY] + a*xjk[YY] + b*xjl[YY];
  temp[ZZ] = xij[ZZ] + a*xjk[ZZ] + b*xjl[ZZ];
  /* 12 flops */
  
  d=c*invsqrt(iprod(temp,temp));
  /* 6 + 10 flops */
  
  x[XX] = xi[XX] + d*temp[XX];
  x[YY] = xi[YY] + d*temp[YY];
  x[ZZ] = xi[ZZ] + d*temp[ZZ];
  /* 6 Flops */
  
  /* TOTAL: 43 flops */
}


void construct_vsites(FILE *log,rvec x[],t_nrnb *nrnb,real dt, 
		       rvec *v,t_idef *idef,t_graph *graph,t_commrec *cr,
		       matrix box,t_comm_vsites *vsitecomm)
{
  rvec      xd,vv;
  real      a1,b1,c1,inv_dt;
  int       i,ii,nra,nrd,tp,ftype;
  t_iatom   avsite,ai,aj,ak,al;
  t_iatom   *ia;
  t_iparams *ip;

  /* Molecules always whole, but I'm not sure whether
   * the periodicity and shift are guaranteed to be consistent
   * between different nodes when running e.g. polymers in
   * parallel. In this special case we thus unshift/shift, but
   * only when necessary. This is to make sure the coordinates
   * we move don't end up a box away...
   */
  if (vsitecomm) {
    if (graph)
      unshift_self(graph,box,x);
    move_construct_x(vsitecomm,x,cr);
    if (graph)
      shift_self(graph,box,x);
  }

  ip     = idef->iparams;
  if (v)
    inv_dt = 1.0/dt;
  else
    inv_dt = 1.0;

  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (interaction_function[ftype].flags & IF_VSITE) {
      nra    = interaction_function[ftype].nratoms;
      nrd    = idef->il[ftype].nr;
      ia     = idef->il[ftype].iatoms;
      
      for(i=0; (i<nrd); ) {
	tp   = ia[0];
	if (ftype != idef->functype[tp]) 
	  gmx_incons("Function types for vsites wrong");
	
	/* The vsite and constructing atoms */
	avsite = ia[1];
	ai   = ia[2];
	aj   = ia[3];

	/* Constants for constructing vsites */
	a1   = ip[tp].vsite.a;
	
	/* Copy the old position */
	copy_rvec(x[avsite],xd);
	
	/* Construct the vsite depending on type */
	switch (ftype) {
	case F_VSITE2:
	  constr_vsite2(x[ai],x[aj],x[avsite],a1);
	  break;
	case F_VSITE3:
	  ak = ia[4];
	  b1 = ip[tp].vsite.b;
	  constr_vsite3(x[ai],x[aj],x[ak],x[avsite],a1,b1);
	  break;
	case F_VSITE3FD:
	  ak = ia[4];
	  b1 = ip[tp].vsite.b;
	  constr_vsite3FD(x[ai],x[aj],x[ak],x[avsite],a1,b1);
	  break;
	case F_VSITE3FAD:
	  ak = ia[4];
	  b1 = ip[tp].vsite.b;
	  constr_vsite3FAD(x[ai],x[aj],x[ak],x[avsite],a1,b1);
	  break;
	case F_VSITE3OUT:
	  ak = ia[4];
	  b1 = ip[tp].vsite.b;
	  c1 = ip[tp].vsite.c;
	  constr_vsite3OUT(x[ai],x[aj],x[ak],x[avsite],a1,b1,c1);
	  break;
	case F_VSITE4FD:
	  ak = ia[4];
	  al = ia[5];
	  b1 = ip[tp].vsite.b;
	  c1 = ip[tp].vsite.c;
	  constr_vsite4FD(x[ai],x[aj],x[ak],x[al],x[avsite],a1,b1,c1);
	  break;
	default:
	  gmx_fatal(FARGS,"No such vsite type %d in %s, line %d",
		      ftype,__FILE__,__LINE__);
	}
	if (v) {
	  /* Calculate velocity of vsite... */
	  rvec_sub(x[avsite],xd,vv);
	  svmul(inv_dt,vv,v[avsite]);
	}
	/* Increment loop variables */
	i  += nra+1;
	ia += nra+1;
      }
    }
  }
  if (vsitecomm) {
    if (graph)
      unshift_self(graph,box,x);
    move_vsite_xv(vsitecomm,x,NULL,cr);
    if (graph)
      shift_self(graph,box,x); /* maybe not necessary */
  }
}

static void spread_vsite2(rvec fi,rvec fj,rvec f,real a)
{
  real fx,fy,fz,b;
  
  b=1.0-a;
  /* 1 flop */
  
  fx=f[XX];
  fy=f[YY];
  fz=f[ZZ];
  fi[XX]+=b*fx;
  fi[YY]+=b*fy;
  fi[ZZ]+=b*fz;
  fj[XX]+=a*fx;
  fj[YY]+=a*fy;
  fj[ZZ]+=a*fz;
  /* 6 Flops */
  
  /* TOTAL: 7 flops */
}

static void spread_vsite3(rvec fi,rvec fj,rvec fk,rvec f,real a,real b)
{
  real fx,fy,fz,c;
  
  c=1.0-a-b;
  /* 2 flops */
  
  fx=f[XX];
  fy=f[YY];
  fz=f[ZZ];
  fi[XX]+=c*fx;
  fi[YY]+=c*fy;
  fi[ZZ]+=c*fz;
  fj[XX]+=a*fx;
  fj[YY]+=a*fy;
  fj[ZZ]+=a*fz;
  fk[XX]+=b*fx;
  fk[YY]+=b*fy;
  fk[ZZ]+=b*fz;
  /* 9 Flops */
  
  /* TOTAL: 11 flops */
}

static void spread_vsite3FD(rvec xi,rvec xj,rvec xk,
			  rvec fi,rvec fj,rvec fk,rvec f,real a,real b)
{
  real fx,fy,fz,c,invl,fproj,a1;
  rvec xij,xjk,xix,temp;
  
  rvec_sub(xj,xi,xij);
  rvec_sub(xk,xj,xjk);
  /* 6 flops */
  
  /* xix goes from i to point x on the line jk */  
  xix[XX]=xij[XX]+a*xjk[XX];
  xix[YY]=xij[YY]+a*xjk[YY];
  xix[ZZ]=xij[ZZ]+a*xjk[ZZ];
  /* 6 flops */
  
  invl=invsqrt(iprod(xix,xix));
  c=b*invl;
  /* 4 + ?10? flops */
  
  fproj=iprod(xix,f)*invl*invl; /* = (xix . f)/(xix . xix) */
  
  fx=f[XX];
  fy=f[YY];
  fz=f[ZZ];
  
  temp[XX]=c*(fx-fproj*xix[XX]);
  temp[YY]=c*(fy-fproj*xix[YY]);
  temp[ZZ]=c*(fz-fproj*xix[ZZ]);
  /* 16 */
  
  /* c is already calculated in constr_vsite3FD
     storing c somewhere will save 26 flops!     */
  
  a1=1-a;
  fi[XX]+=fx-temp[XX];
  fi[YY]+=fy-temp[YY];
  fi[ZZ]+=fz-temp[ZZ];
  fj[XX]+=a1*temp[XX];
  fj[YY]+=a1*temp[YY];
  fj[ZZ]+=a1*temp[ZZ];
  fk[XX]+= a*temp[XX];
  fk[YY]+= a*temp[YY];
  fk[ZZ]+= a*temp[ZZ];
  /* 19 Flops */
  
  /* TOTAL: 61 flops */
}

static void spread_vsite3FAD(rvec xi,rvec xj,rvec xk,
			   rvec fi,rvec fj,rvec fk,rvec f,real a,real b)
{
  rvec xij,xjk,xperp,Fpij,Fppp,f1,f2,f3;
  real a1,b1,c1,c2,invdij,invdij2,invdp,fproj;
  int d;
  
  rvec_sub(xj,xi,xij);
  rvec_sub(xk,xj,xjk);
  /* 6 flops */
  
  invdij = invsqrt(iprod(xij,xij));
  invdij2 = invdij * invdij;
  c1 = iprod(xij,xjk) * invdij2;
  xperp[XX] = xjk[XX] - c1*xij[XX];
  xperp[YY] = xjk[YY] - c1*xij[YY];
  xperp[ZZ] = xjk[ZZ] - c1*xij[ZZ];
  /* xperp in plane ijk, perp. to ij */
  invdp = invsqrt(iprod(xperp,xperp));
  a1 = a*invdij;
  b1 = b*invdp;
  /* 45 flops */
  
  /* a1, b1 and c1 are already calculated in constr_vsite3FAD
     storing them somewhere will save 45 flops!     */
  
  fproj=iprod(xij  ,f)*invdij2;
  svmul(fproj,                     xij,  Fpij); /* proj. f on xij */
  svmul(iprod(xperp,f)*invdp*invdp,xperp,Fppp); /* proj. f on xperp */
  svmul(b1*fproj,                  xperp,f3);
  /* 23 flops */
  
  rvec_sub(f,Fpij,f1);  /* f1 = f - Fpij */
  rvec_sub(f1,Fppp,f2); /* f2 = f - Fpij - Fppp */
  for (d=0; (d<DIM); d++) {
    f1[d]*=a1;
    f2[d]*=b1;
  }
  /* 12 flops */
  
  c2=1+c1;
  fi[XX] += f[XX] - f1[XX] + c1*f2[XX] + f3[XX];
  fi[YY] += f[YY] - f1[YY] + c1*f2[YY] + f3[YY];
  fi[ZZ] += f[ZZ] - f1[ZZ] + c1*f2[ZZ] + f3[ZZ];
  fj[XX] +=         f1[XX] - c2*f2[XX] - f3[XX];
  fj[YY] +=         f1[YY] - c2*f2[YY] - f3[YY];
  fj[ZZ] +=         f1[ZZ] - c2*f2[ZZ] - f3[ZZ];
  fk[XX] +=                     f2[XX];
  fk[YY] +=                     f2[YY];
  fk[ZZ] +=                     f2[ZZ];
  /* 30 Flops */
  
  /* TOTAL: 113 flops */
}

static void spread_vsite3OUT(rvec xi,rvec xj,rvec xk,
			   rvec fi,rvec fj,rvec fk,rvec f,real a,real b,real c)
{
  rvec xij,xik,ffj,ffk;
  real cfx,cfy,cfz;
  int  m;
  
  rvec_sub(xj,xi,xij);
  rvec_sub(xk,xi,xik);
  /* 6 Flops */
  
  cfx=c*f[XX];
  cfy=c*f[YY];
  cfz=c*f[ZZ];
  /* 3 Flops */
  
  ffj[XX] =  a*f[XX]     - xik[ZZ]*cfy + xik[YY]*cfz;
  ffj[YY] =  xik[ZZ]*cfx + a*f[YY]     - xik[XX]*cfz;
  ffj[ZZ] = -xik[YY]*cfx + xik[XX]*cfy + a*f[ZZ];
  
  ffk[XX] =  b*f[XX]     + xij[ZZ]*cfy - xij[YY]*cfz;
  ffk[YY] = -xij[ZZ]*cfx + b*f[YY]     + xij[XX]*cfz;
  ffk[ZZ] =  xij[YY]*cfx - xij[XX]*cfy + b*f[ZZ];
  /* 30 Flops */
    
  for(m=0; (m<DIM); m++) {
    fi[m]+=f[m]-ffj[m]-ffk[m];
    fj[m]+=ffj[m];
    fk[m]+=ffk[m];
  }
  /* 15 Flops */
  
  /* TOTAL: 54 flops */
}

static void spread_vsite4FD(rvec xi,rvec xj,rvec xk,rvec xl,
			  rvec fi,rvec fj,rvec fk,rvec fl,rvec f,
			  real a,real b,real c)
{
  real fx,fy,fz,d,invl,fproj,a1;
  rvec xij,xjk,xjl,xix,temp;
  
  rvec_sub(xj,xi,xij);
  rvec_sub(xk,xj,xjk);
  rvec_sub(xl,xj,xjl);
  /* 9 flops */
  
  /* xix goes from i to point x on the plane jkl */  
  xix[XX]=xij[XX]+a*xjk[XX]+b*xjl[XX];
  xix[YY]=xij[YY]+a*xjk[YY]+b*xjl[YY];
  xix[ZZ]=xij[ZZ]+a*xjk[ZZ]+b*xjl[ZZ];
  /* 12 flops */
  
  invl=invsqrt(iprod(xix,xix));
  d=c*invl;
  /* 4 + ?10? flops */
  
  fproj=iprod(xix,f)*invl*invl; /* = (xix . f)/(xix . xix) */
  
  fx=f[XX];
  fy=f[YY];
  fz=f[ZZ];
  
  temp[XX]=d*(fx-fproj*xix[XX]);
  temp[YY]=d*(fy-fproj*xix[YY]);
  temp[ZZ]=d*(fz-fproj*xix[ZZ]);
  /* 16 */
  
  /* c is already calculated in constr_vsite3FD
     storing c somewhere will save 35 flops!     */
  
  a1=1-a-b;
  fi[XX]+=fx-temp[XX];
  fi[YY]+=fy-temp[YY];
  fi[ZZ]+=fz-temp[ZZ];
  fj[XX]+=a1*temp[XX];
  fj[YY]+=a1*temp[YY];
  fj[ZZ]+=a1*temp[ZZ];
  fk[XX]+= a*temp[XX];
  fk[YY]+= a*temp[YY];
  fk[ZZ]+= a*temp[ZZ];
  fl[XX]+= b*temp[XX];
  fl[YY]+= b*temp[YY];
  fl[ZZ]+= b*temp[ZZ];
  /* 26 Flops */
  
  /* TOTAL: 77 flops */
}

void spread_vsite_f(FILE *log,rvec x[],rvec f[],t_nrnb *nrnb,t_idef *idef,
		    t_comm_vsites *vsitecomm,t_commrec *cr)
{
  real      a1,b1,c1;
  int       i,m,nra,nrd,tp,ftype;
  int       nd2,nd3,nd3FD,nd3FAD,nd3OUT,nd4FD;
  t_iatom   avsite,ai,aj,ak,al;
  t_iatom   *ia;
  t_iparams *ip;
  
  /* We only move forces here, and they are independent of shifts */
  if (vsitecomm)
    move_vsite_f(vsitecomm,f,cr);

  ip     = idef->iparams;

  nd2    = 0;
  nd3    = 0;
  nd3FD  = 0;
  nd3FAD = 0;
  nd3OUT = 0;
  nd4FD  = 0;
   
  /* this loop goes backwards to be able to build *
   * higher type vsites from lower types         */
  for(ftype=F_NRE-1; (ftype>=0); ftype--) {
    if (interaction_function[ftype].flags & IF_VSITE) {
      nra    = interaction_function[ftype].nratoms;
      nrd    = idef->il[ftype].nr;
      ia     = idef->il[ftype].iatoms;
      
      for(i=0; (i<nrd); ) {
	tp   = ia[0];
	if (ftype != idef->functype[tp])
	  gmx_incons("Functiontypes for vsites wrong");
	
	/* The vsite and constructing atoms */
	avsite = ia[1];
	ai   = ia[2];
	aj   = ia[3];
		
	/* Constants for constructing */
	a1   = ip[tp].vsite.a; 
      
	/* Construct the vsite depending on type */
	switch (ftype) {
	case F_VSITE2:
	  spread_vsite2(f[ai],f[aj],f[avsite],a1);
	  nd2++;
	  break;
	case F_VSITE3:
	  ak = ia[4];
	  b1 = ip[tp].vsite.b;
	  spread_vsite3(f[ai],f[aj],f[ak],f[avsite],a1,b1);
	  nd3++;
	  break;
	case F_VSITE3FD:
	  ak = ia[4];
	  b1 = ip[tp].vsite.b;
	  spread_vsite3FD(x[ai],x[aj],x[ak],f[ai],f[aj],f[ak],f[avsite],a1,b1);
	  nd3FD++;
	  break;
	case F_VSITE3FAD:
	  ak = ia[4];
	  b1 = ip[tp].vsite.b;
	  spread_vsite3FAD(x[ai],x[aj],x[ak],f[ai],f[aj],f[ak],f[avsite],a1,b1);
	  nd3FAD++;
	  break;
	case F_VSITE3OUT:
	  ak = ia[4];
	  b1 = ip[tp].vsite.b;
	  c1 = ip[tp].vsite.c;
	  spread_vsite3OUT(x[ai],x[aj],x[ak],f[ai],f[aj],f[ak],f[avsite],a1,b1,c1);
	  nd3OUT++;
	  break;
	case F_VSITE4FD:
	  ak = ia[4];
	  al = ia[5];
	  b1 = ip[tp].vsite.b;
	  c1 = ip[tp].vsite.c;
	  spread_vsite4FD(x[ai],x[aj],x[ak],x[al],
			f[ai],f[aj],f[ak],f[al],f[avsite],a1,b1,c1);
	  nd4FD++;
	  break;
	default:
	  gmx_fatal(FARGS,"No such vsite type %d in %s, line %d",
		      ftype,__FILE__,__LINE__);
	}
	clear_rvec(f[avsite]);
	
	/* Increment loop variables */
	i  += nra+1;
	ia += nra+1;
      }
    }
  }
  
  inc_nrnb(nrnb,eNR_VSITE2,   nd2   );
  inc_nrnb(nrnb,eNR_VSITE3,   nd3   );
  inc_nrnb(nrnb,eNR_VSITE3FD, nd3FD );
  inc_nrnb(nrnb,eNR_VSITE3FAD,nd3FAD);
  inc_nrnb(nrnb,eNR_VSITE3OUT,nd3OUT);
  inc_nrnb(nrnb,eNR_VSITE4FD, nd4FD );

  /* We only move forces here, and they are independent of shifts */
  if(vsitecomm)
    move_construct_f(vsitecomm,f,cr);
}

