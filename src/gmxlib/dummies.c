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
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_dummies_c = "$Id$";

#include <stdio.h>
#include "typedefs.h"
#include "assert.h"
#include "dummies.h"
#include "macros.h"
#include "smalloc.h"
#include "nrnb.h"
#include "vec.h"

void set_dummies_ptype(bool bVerbose, t_topology *sys)
{
  int i,ftype;
  int nra,nrd,tp;
  t_iatom *ia,adum;
  
  if (bVerbose)
    fprintf(stderr,"Setting particle type to Dummy for dummy atoms\n");
  if (debug)
    fprintf(stderr,"checking %d functypes\n",F_NRE);
  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (interaction_function[ftype].flags & IF_DUMMY) {
      nra    = interaction_function[ftype].nratoms;
      nrd    = sys->idef.il[ftype].nr;
      ia     = sys->idef.il[ftype].iatoms;
      
      if (debug)
	fprintf(stderr,"doing %d dummies(%d)\n",(nrd / (nra+1)),ftype);
      
      for(i=0; (i<nrd); ) {
	tp   = ia[0];
	assert(ftype == sys->idef.functype[tp]);
	
	/* The dummy atom */
	adum = ia[1];
	sys->atoms.atom[adum].ptype=eptDummy;
	
	i  += nra+1;
	ia += nra+1;
      }
    }
  }
}

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
		       rvec *v,t_idef *idef)
{
  rvec      xd,vv;
  real      a1,b1,c1,inv_dt;
  int       i,nra,nrd,tp,ftype;
  t_iatom   adum,ai,aj,ak,al;
  t_iatom   *ia;
  t_iparams *ip;
  
  ip     = idef->iparams;
  if (v)
    inv_dt = 1.0/dt;
  
  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (interaction_function[ftype].flags & IF_DUMMY) {
      nra    = interaction_function[ftype].nratoms;
      nrd    = idef->il[ftype].nr;
      ia     = idef->il[ftype].iatoms;
      
      for(i=0; (i<nrd); ) {
	tp   = ia[0];
	assert(ftype == idef->functype[tp]);
	
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
	  fatal_error(0,"No such dummy type %d in %s, line %d",
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

void spread_dummy_f(FILE *log,rvec x[],rvec f[],t_nrnb *nrnb,t_idef *idef)
{
  real      a1,b1,c1;
  int       i,nra,nrd,tp,ftype;
  int       nd2,nd3,nd3FD,nd3FAD,nd3OUT,nd4FD;
  t_iatom   adum,ai,aj,ak,al;
  t_iatom   *ia;
  t_iparams *ip;
  
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
	assert(ftype == idef->functype[tp]);
	
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
	  fatal_error(0,"No such dummy type %d in %s, line %d",
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
}

