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

static void constr_dum1(rvec xa,rvec xb,rvec x,real a,real b)
{
  x[XX]=a*xa[XX]+b*xb[XX];
  x[YY]=a*xa[YY]+b*xb[YY];
  x[ZZ]=a*xa[ZZ]+b*xb[ZZ];
  /* 9 Flops */
}

static void constr_dum2(rvec xa,rvec xb,rvec xc,rvec x,real a,real b,real c)
{
  x[XX]=a*xa[XX]+b*xb[XX]+c*xc[XX];
  x[YY]=a*xa[YY]+b*xb[YY]+c*xc[YY];
  x[ZZ]=a*xa[ZZ]+b*xb[ZZ]+c*xc[ZZ];
  /* 15 Flops */
}

static void constr_dum3(rvec x1,rvec x2,rvec x3,rvec x,real a,real b,real c)
{
  rvec x21,x31,temp;
  int  m;
  
  rvec_sub(x2,x1,x21);
  rvec_sub(x3,x1,x31);
  oprod(x21,x31,temp);
  /* 15 Flops */
  
  for(m=0; (m<DIM); m++) 
    x[m]=x1[m]+a*x21[m]+b*x31[m]+c*temp[m];
  /* 18 Flops */
}

void construct_dummies(FILE *log,rvec x[],t_nrnb *nrnb,real dt, 
		       rvec v[],t_idef *idef)
{
  rvec      xd,vv;
  real      a1,b1,c1,inv_dt;
  int       i,nra,nrd,tp,ftype;
  t_iatom   adum,ai,aj,ak;
  t_iatom   *ia;
  t_iparams *ip;
  
  ip     = idef->iparams;
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
	b1   = ip[tp].dummy.b;
	
	/* Copy the old position */
	copy_rvec(x[adum],xd);
	
	/* Construct the dummy depending on type */
	switch (ftype) {
	case F_DUMMY1:
	  constr_dum1(x[ai],x[aj],x[adum],a1,b1);
	  break;
	case F_DUMMY2:
	  ak = ia[4];
	  c1 = ip[tp].dummy.c;
	  constr_dum2(x[ai],x[aj],x[ak],x[adum],a1,b1,c1);
	  break;
	case F_DUMMY3:
	  ak = ia[4];
	  c1 = ip[tp].dummy.c;
	  constr_dum3(x[ai],x[aj],x[ak],x[adum],a1,b1,c1);
	  break;
	default:
	  fatal_error(0,"No such dummy type %d in %s, line %d",
		      ftype,__FILE__,__LINE__);
	}
	/* Calculate velocity of dummy... */
	rvec_sub(x[adum],xd,vv);
	svmul(inv_dt,vv,v[adum]);
	
	/* Increment loop variables */
	i  += nra+1;
	ia += nra+1;
      }
    }
  }
}


static void spread_dum1(rvec fa,rvec fb,rvec f,real a,real b)
{
  real fx,fy,fz;
  
  fx=f[XX];
  fy=f[YY];
  fz=f[ZZ];
  fa[XX]+=a*fx;
  fa[YY]+=a*fy;
  fa[ZZ]+=a*fz;
  fb[XX]+=b*fx;
  fb[YY]+=b*fy;
  fb[ZZ]+=b*fz;
  /* 6 Flops */
}

static void spread_dum2(rvec fa,rvec fb,rvec fc,rvec f,real a,real b,real c)
{
  real fx,fy,fz;
  
  fx=f[XX];
  fy=f[YY];
  fz=f[ZZ];
  fa[XX]+=a*fx;
  fa[YY]+=a*fy;
  fa[ZZ]+=a*fz;
  fb[XX]+=b*fx;
  fb[YY]+=b*fy;
  fb[ZZ]+=b*fz;
  fc[XX]+=c*fx;
  fc[YY]+=c*fy;
  fc[ZZ]+=c*fz;
  /* 9 Flops */
}

static void spread_dum3(rvec x1,rvec x2,rvec x3,
			rvec fa,rvec fb,rvec fc,rvec f,real a,real b,real c)
{
  rvec x21,x31,ffb,ffc;
  real cfx,cfy,cfz;
  int  m;
  
  rvec_sub(x2,x1,x21);
  rvec_sub(x3,x1,x31);
  /* 6 Flops */
    
  cfx=c*f[XX];
  cfy=c*f[YY];
  cfz=c*f[ZZ];
  /* 3 Flops */
  
  ffb[XX] =  a*f[XX]     - x31[ZZ]*cfy + x31[YY]*cfz;
  ffb[YY] =  x31[ZZ]*cfx + a*f[YY]     - x31[XX]*cfz;
  ffb[ZZ] = -x31[YY]*cfx + x31[XX]*cfy + a*f[ZZ];
  
  ffc[XX] =  b*f[XX]     + x21[ZZ]*cfy - x21[YY]*cfz;
  ffc[YY] = -x21[ZZ]*cfx + b*f[YY]     + x21[XX]*cfz;
  ffc[ZZ] =  x21[YY]*cfx - x21[XX]*cfy + b*f[ZZ];
  /* 30 Flops */
    
  for(m=0; (m<DIM); m++) {
    fb[m]+=ffb[m];
    fc[m]+=ffc[m];
    fa[m]+=f[m]-ffb[m]-ffc[m];
  }
  /* 15 Flops */
}

void spread_dummy_f(FILE *log,rvec x[],rvec f[],t_nrnb *nrnb,t_idef *idef)
{
  real      a1,b1,c1;
  int       i,nra,nrd,tp,ftype;
  int       nd1,nd2,nd3;
  t_iatom   adum,ai,aj,ak;
  t_iatom   *ia;
  t_iparams *ip;
  
  ip     = idef->iparams;
  
  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (interaction_function[ftype].flags & IF_DUMMY) {
      nra    = interaction_function[ftype].nratoms;
      nrd    = idef->il[ftype].nr;
      ia     = idef->il[ftype].iatoms;
      nd1    = 0;
      nd2    = 0;
      nd3    = 0;
      
      for(i=0; (i<nrd); ) {
	tp   = ia[0];
	assert(ftype == idef->functype[tp]);
	
	/* The dummy and constructing atoms */
	adum = ia[1];
	ai   = ia[2];
	aj   = ia[3];
	
	/* Constants for constructing */
	a1   = ip[tp].dummy.a;
	b1   = ip[tp].dummy.b;
	
	/* Construct the dummy depending on type */
	switch (ftype) {
	case F_DUMMY1:
	  spread_dum1(f[ai],f[aj],f[adum],a1,b1);
	  nd1++;
	  break;
	case F_DUMMY2:
	  ak = ia[4];
	  c1 = ip[tp].dummy.c;
	  spread_dum2(f[ai],f[aj],f[ak],f[adum],a1,b1,c1);
	  nd2++;
	  break;
	case F_DUMMY3:
	  ak = ia[4];
	  c1 = ip[tp].dummy.c;
	  spread_dum3(x[ai],x[aj],x[ak],f[ai],f[aj],f[ak],f[adum],a1,b1,c1);
	  nd3++;
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
  inc_nrnb(nrnb,eNR_DUM1,nd1);
  inc_nrnb(nrnb,eNR_DUM2,nd2);
  inc_nrnb(nrnb,eNR_DUM3,nd3);
}

