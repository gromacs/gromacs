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
#include "dummies.h"
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
 * of constructing atoms and dummies. This is necessary
 * when dummy constructs cross node borders (unavoidable
 * for e.g. polymers using anisotropic united atoms).
 */
/* Communication routines for dummies. The coordinates and
 * forces are only move on a need-to-know basis, usually only
 * 2-3 atoms per processor. To achieve this small amount of
 * communication, and to limit it to nearest neighbour messages,
 * we demand that dummies are not spread over nonadjacent nodes.
 * Thus, keep your dummies close to your constructing atoms.
 * (mdrun & grompp will report an error otherwise)
 */ 


static void move_construct_x(t_comm_dummies *dummycomm, rvec x[], t_commrec *cr)
{
  static bool bFirst=TRUE;
  int i;
   
  if (bFirst) {
    /* Make the larger than necessary to avoid cache sharing */
    snew(nextbuf,2*(dummycomm->nnextdum+dummycomm->nnextconstr)+100);
    snew(prevbuf,2*(dummycomm->nprevdum+dummycomm->nprevconstr)+100);
    bFirst=FALSE;
  }
   
  /* package coords to send left. Dummy coords are needed to create v */
  for(i=0;i<dummycomm->nprevconstr;i++)
    copy_rvec(x[dummycomm->idxprevconstr[i]],prevbuf[i]);
  for(i=0;i<dummycomm->nprevdum;i++)
    copy_rvec(x[dummycomm->idxprevdum[i]],prevbuf[dummycomm->nprevconstr+i]);
  
  /* send them off, and recieve from the right */
  if(dummycomm->nprevconstr>0 || dummycomm->nprevdum>0)
    gmx_tx(cr->left,prevbuf,
	   sizeof(rvec)*(dummycomm->nprevconstr+dummycomm->nprevdum));
  
  if(dummycomm->nnextconstr>0 || dummycomm->nnextdum>0)
    gmx_rx(cr->right,nextbuf,
	   sizeof(rvec)*(dummycomm->nnextconstr+dummycomm->nnextdum));
  
  if(dummycomm->nprevconstr>0 || dummycomm->nprevdum>0)
    gmx_tx_wait(cr->left);
  
  if(dummycomm->nnextconstr>0 || dummycomm->nnextdum>0)
    gmx_rx_wait(cr->right);
  
  /* Put them where they belong */
  for(i=0;i<dummycomm->nnextconstr;i++)
    copy_rvec(nextbuf[i],x[dummycomm->idxnextconstr[i]]);
  for(i=0;i<dummycomm->nnextdum;i++)
    copy_rvec(nextbuf[dummycomm->nnextconstr+i],
	      x[dummycomm->idxnextdum[i]]);

  /* Now we are ready to do the constructing business ! */
}


static void move_dummy_xv(t_comm_dummies *dummycomm, rvec x[], rvec v[],t_commrec *cr)
{
  int i;
  int sendsize,recvsize;

  sendsize=sizeof(rvec)*dummycomm->nnextdum;
  recvsize=sizeof(rvec)*dummycomm->nprevdum;

  if(v!=NULL) {
    sendsize=sendsize*2;
    recvsize=recvsize*2;
  }

  /* Package nonlocal constructed dummies */
  for(i=0;i<dummycomm->nnextdum;i++)
    copy_rvec(x[dummycomm->idxnextdum[i]],nextbuf[i]);

  if(v!=NULL)
    for(i=0;i<dummycomm->nnextdum;i++)
      copy_rvec(v[dummycomm->idxnextdum[i]],nextbuf[dummycomm->nnextdum+i]);
  
  /* send them off, and recieve from the right */
  if(dummycomm->nnextdum>0)
    gmx_tx(cr->right,nextbuf,sendsize);
  
  if(dummycomm->nprevdum>0)
    gmx_rx(cr->left,prevbuf,recvsize);
  
  if(dummycomm->nnextdum>0)
    gmx_tx_wait(cr->right);
  
  if(dummycomm->nprevdum>0)
    gmx_rx_wait(cr->left);
  
  /* Put them where they belong */
  for(i=0;i<dummycomm->nprevdum;i++)
    copy_rvec(prevbuf[i],x[dummycomm->idxprevdum[i]]);

  if(v!=NULL)
    for(i=0;i<dummycomm->nprevdum;i++)
      copy_rvec(prevbuf[dummycomm->nprevdum+i],v[dummycomm->idxprevdum[i]]);

  /* All coordinates are in place on the respective home node now */
}

static void move_dummy_f(t_comm_dummies *dummycomm, rvec f[], t_commrec *cr)
{
  int i;

  /* package dummy particle forces to send left */
  for(i=0;i<dummycomm->nprevdum;i++)
    copy_rvec(f[dummycomm->idxprevdum[i]],prevbuf[i]);

  /* off they go! - but only if there is something to send! */
  if(dummycomm->nprevdum>0)
    gmx_tx(cr->left,prevbuf,sizeof(rvec)*dummycomm->nprevdum);

  /* Get our share from the right, if there is anything to have */
  if(dummycomm->nnextdum>0)
    gmx_rx(cr->right,nextbuf,sizeof(rvec)*dummycomm->nnextdum);
  
  if(dummycomm->nprevdum>0)
    gmx_tx_wait(cr->left);
  
  if(dummycomm->nnextdum>0)
    gmx_rx_wait(cr->right);

  /* Put them where they belong */
  for(i=0;i<dummycomm->nnextdum;i++)
    copy_rvec(nextbuf[i],f[dummycomm->idxnextdum[i]]);
  
  /* Zero forces on nonlocal constructing atoms.
   * This is necessary since dummy force spreading is done
   * after the normal force addition, and we don't want
   * to include them twice.
   * (They have already been added on the home node).
   */
  for(i=0;i<dummycomm->nnextconstr;i++)
    clear_rvec(f[dummycomm->idxnextconstr[i]]);
}

static void move_construct_f(t_comm_dummies *dummycomm, rvec f[], t_commrec *cr)
{
  int i;

  /* Spread forces to nonlocal constructing atoms.
   */
  /* package forces to send right */
  for(i=0;i<dummycomm->nnextconstr;i++)
    copy_rvec(f[dummycomm->idxnextconstr[i]],nextbuf[i]);
  
  /* send them off, and recieve from the right */
  if(dummycomm->nnextconstr>0)
    gmx_tx(cr->right,nextbuf,sizeof(rvec)*dummycomm->nnextconstr);
  
  if(dummycomm->nprevconstr>0)
    gmx_rx(cr->left,prevbuf,sizeof(rvec)*dummycomm->nprevconstr);
  
  if(dummycomm->nnextconstr>0)
    gmx_tx_wait(cr->right);

  if(dummycomm->nprevconstr>0)
    gmx_rx_wait(cr->left);
  
  /* Add them where they belong */
  for(i=0;i<dummycomm->nprevconstr;i++)
    rvec_inc(f[dummycomm->idxprevconstr[i]],prevbuf[i]);
  
  /* Zero nonlocal dummies */
  for(i=0;i<dummycomm->nprevdum;i++)
    clear_rvec(f[dummycomm->idxprevdum[i]]);
  
  /* All forces are on the home processor now */  
}


/* Dummy construction routines */

static void constr_dum2(rvec xi,rvec xj,rvec x,real a)
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

static void constr_dum3(rvec xi,rvec xj,rvec xk,rvec x,real a,real b)
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

static void constr_dum3FD(rvec xi,rvec xj,rvec xk,rvec x,real a,real b)
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

static void constr_dum3FAD(rvec xi,rvec xj,rvec xk,rvec x,real a,real b)
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

static void constr_dum3OUT(rvec xi,rvec xj,rvec xk,rvec x,real a,real b,real c)
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

static void constr_dum4FD(rvec xi,rvec xj,rvec xk,rvec xl,rvec x,
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


void construct_dummies(FILE *log,rvec x[],t_nrnb *nrnb,real dt, 
		       rvec *v,t_idef *idef,t_graph *graph,t_commrec *cr,
		       matrix box,t_comm_dummies *dummycomm)
{
  rvec      xd,vv;
  real      a1,b1,c1,inv_dt;
  int       i,ii,nra,nrd,tp,ftype;
  t_iatom   adum,ai,aj,ak,al;
  t_iatom   *ia;
  t_iparams *ip;

  /* Molecules always whole, but I'm not sure whether
   * the periodicity and shift are guaranteed to be consistent
   * between different nodes when running e.g. polymers in
   * parallel. In this special case we thus unshift/shift, but
   * only when necessary. This is to make sure the coordinates
   * we move don't end up a box away...
   */
  if (dummycomm) {
    if (graph)
      unshift_self(graph,box,x);
    move_construct_x(dummycomm,x,cr);
    if (graph)
      shift_self(graph,box,x);
  }

  ip     = idef->iparams;
  if (v)
    inv_dt = 1.0/dt;
  else
    inv_dt = 1.0;

  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (interaction_function[ftype].flags & IF_DUMMY) {
      nra    = interaction_function[ftype].nratoms;
      nrd    = idef->il[ftype].nr;
      ia     = idef->il[ftype].iatoms;
      
      for(i=0; (i<nrd); ) {
	tp   = ia[0];
	if (ftype != idef->functype[tp]) 
	  gmx_incons("Function types for dummies wrong");
	
	/* The dummy and constructing atoms */
	adum = ia[1];
	ai   = ia[2];
	aj   = ia[3];

	/* Constants for constructing dummies */
	a1   = ip[tp].dummy.a;
	
	/* Copy the old position */
	copy_rvec(x[adum],xd);
	
	/* Construct the dummy depending on type */
	switch (ftype) {
	case F_DUMMY2:
	  constr_dum2(x[ai],x[aj],x[adum],a1);
	  break;
	case F_DUMMY3:
	  ak = ia[4];
	  b1 = ip[tp].dummy.b;
	  constr_dum3(x[ai],x[aj],x[ak],x[adum],a1,b1);
	  break;
	case F_DUMMY3FD:
	  ak = ia[4];
	  b1 = ip[tp].dummy.b;
	  constr_dum3FD(x[ai],x[aj],x[ak],x[adum],a1,b1);
	  break;
	case F_DUMMY3FAD:
	  ak = ia[4];
	  b1 = ip[tp].dummy.b;
	  constr_dum3FAD(x[ai],x[aj],x[ak],x[adum],a1,b1);
	  break;
	case F_DUMMY3OUT:
	  ak = ia[4];
	  b1 = ip[tp].dummy.b;
	  c1 = ip[tp].dummy.c;
	  constr_dum3OUT(x[ai],x[aj],x[ak],x[adum],a1,b1,c1);
	  break;
	case F_DUMMY4FD:
	  ak = ia[4];
	  al = ia[5];
	  b1 = ip[tp].dummy.b;
	  c1 = ip[tp].dummy.c;
	  constr_dum4FD(x[ai],x[aj],x[ak],x[al],x[adum],a1,b1,c1);
	  break;
	default:
	  gmx_fatal(FARGS,"No such dummy type %d in %s, line %d",
		      ftype,__FILE__,__LINE__);
	}
	if (v) {
	  /* Calculate velocity of dummy... */
	  rvec_sub(x[adum],xd,vv);
	  svmul(inv_dt,vv,v[adum]);
	}
	/* Increment loop variables */
	i  += nra+1;
	ia += nra+1;
      }
    }
  }
  if (dummycomm) {
    if (graph)
      unshift_self(graph,box,x);
    move_dummy_xv(dummycomm,x,NULL,cr);
    if (graph)
      shift_self(graph,box,x); /* maybe not necessary */
  }
}

static void spread_dum2(rvec fi,rvec fj,rvec f,real a)
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

static void spread_dum3(rvec fi,rvec fj,rvec fk,rvec f,real a,real b)
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

static void spread_dum3FD(rvec xi,rvec xj,rvec xk,
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
  
  /* c is already calculated in constr_dum3FD
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

static void spread_dum3FAD(rvec xi,rvec xj,rvec xk,
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
  
  /* a1, b1 and c1 are already calculated in constr_dum3FAD
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

static void spread_dum3OUT(rvec xi,rvec xj,rvec xk,
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

static void spread_dum4FD(rvec xi,rvec xj,rvec xk,rvec xl,
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
  
  /* c is already calculated in constr_dum3FD
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

void spread_dummy_f(FILE *log,rvec x[],rvec f[],t_nrnb *nrnb,t_idef *idef,
		    t_comm_dummies *dummycomm,t_commrec *cr)
{
  real      a1,b1,c1;
  int       i,m,nra,nrd,tp,ftype;
  int       nd2,nd3,nd3FD,nd3FAD,nd3OUT,nd4FD;
  t_iatom   adum,ai,aj,ak,al;
  t_iatom   *ia;
  t_iparams *ip;
  
  /* We only move forces here, and they are independent of shifts */
  if (dummycomm)
    move_dummy_f(dummycomm,f,cr);

  ip     = idef->iparams;

  nd2    = 0;
  nd3    = 0;
  nd3FD  = 0;
  nd3FAD = 0;
  nd3OUT = 0;
  nd4FD  = 0;
   
  /* this loop goes backwards to be able to build *
   * higher type dummies from lower types         */
  for(ftype=F_NRE-1; (ftype>=0); ftype--) {
    if (interaction_function[ftype].flags & IF_DUMMY) {
      nra    = interaction_function[ftype].nratoms;
      nrd    = idef->il[ftype].nr;
      ia     = idef->il[ftype].iatoms;
      
      for(i=0; (i<nrd); ) {
	tp   = ia[0];
	if (ftype != idef->functype[tp])
	  gmx_incons("Functiontypes for dummies wrong");
	
	/* The dummy and constructing atoms */
	adum = ia[1];
	ai   = ia[2];
	aj   = ia[3];
		
	/* Constants for constructing */
	a1   = ip[tp].dummy.a; 
      
	/* Construct the dummy depending on type */
	switch (ftype) {
	case F_DUMMY2:
	  spread_dum2(f[ai],f[aj],f[adum],a1);
	  nd2++;
	  break;
	case F_DUMMY3:
	  ak = ia[4];
	  b1 = ip[tp].dummy.b;
	  spread_dum3(f[ai],f[aj],f[ak],f[adum],a1,b1);
	  nd3++;
	  break;
	case F_DUMMY3FD:
	  ak = ia[4];
	  b1 = ip[tp].dummy.b;
	  spread_dum3FD(x[ai],x[aj],x[ak],f[ai],f[aj],f[ak],f[adum],a1,b1);
	  nd3FD++;
	  break;
	case F_DUMMY3FAD:
	  ak = ia[4];
	  b1 = ip[tp].dummy.b;
	  spread_dum3FAD(x[ai],x[aj],x[ak],f[ai],f[aj],f[ak],f[adum],a1,b1);
	  nd3FAD++;
	  break;
	case F_DUMMY3OUT:
	  ak = ia[4];
	  b1 = ip[tp].dummy.b;
	  c1 = ip[tp].dummy.c;
	  spread_dum3OUT(x[ai],x[aj],x[ak],f[ai],f[aj],f[ak],f[adum],a1,b1,c1);
	  nd3OUT++;
	  break;
	case F_DUMMY4FD:
	  ak = ia[4];
	  al = ia[5];
	  b1 = ip[tp].dummy.b;
	  c1 = ip[tp].dummy.c;
	  spread_dum4FD(x[ai],x[aj],x[ak],x[al],
			f[ai],f[aj],f[ak],f[al],f[adum],a1,b1,c1);
	  nd4FD++;
	  break;
	default:
	  gmx_fatal(FARGS,"No such dummy type %d in %s, line %d",
		      ftype,__FILE__,__LINE__);
	}
	clear_rvec(f[adum]);
	
	/* Increment loop variables */
	i  += nra+1;
	ia += nra+1;
      }
    }
  }
  
  inc_nrnb(nrnb,eNR_DUM2,   nd2   );
  inc_nrnb(nrnb,eNR_DUM3,   nd3   );
  inc_nrnb(nrnb,eNR_DUM3FD, nd3FD );
  inc_nrnb(nrnb,eNR_DUM3FAD,nd3FAD);
  inc_nrnb(nrnb,eNR_DUM3OUT,nd3OUT);
  inc_nrnb(nrnb,eNR_DUM4FD, nd4FD );

  /* We only move forces here, and they are independent of shifts */
  if(dummycomm)
    move_construct_f(dummycomm,f,cr);
}

