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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "physics.h"
#include "vec.h"
#include "maths.h"
#include "txtdump.h"
#include "bondf.h"
#include "smalloc.h"
#include "pbc.h"
#include "ns.h"
#include "macros.h"
#include "names.h"
#include "fatal.h"
#include "mshift.h"
#include "main.h"
#include "disre.h"
#include "orires.h"
#include "force.h"
#include "mdrun.h"

static int pbc_rvec_sub(const t_pbc *pbc,const rvec xi,const rvec xj,rvec dx)
{
  if (pbc) {
    return pbc_dx(pbc,xi,xj,dx);
  }
  else {
    rvec_sub(xi,xj,dx);
    return CENTRAL;
  }
}

void calc_bonds(FILE *fplog,const t_commrec *cr,const t_commrec *mcr,
		const t_idef *idef,
		rvec x[],rvec f[],
		t_forcerec *fr,const t_pbc *pbc,const t_graph *g,
		real epot[],t_nrnb *nrnb,
		real lambda,
		const t_mdatoms *md,int ngrp,real egnb[],real egcoul[],
		t_fcdata *fcd,
		int step,bool bSepDVDL)
{
  int    ftype,nbonds,ind,nat;
  real   v,dvdl;
  const  t_pbc *pbc_null;

  if (fr->ePBC == epbcFULL)
    pbc_null = pbc;
  else
    pbc_null = NULL;

  if (bSepDVDL)
    fprintf(fplog,"Step %d: bonded V and dVdl for node %d:\n",step,cr->nodeid);

#ifdef DEBUG
  if (g)
    p_graph(debug,"Bondage is fun",g);
#endif
  
  /* Do pre force calculation stuff which might require communication */
  if (idef->il[F_ORIRES].nr)
    epot[F_ORIRESDEV] = calc_orires_dev(mcr,idef->il[F_ORIRES].nr,
					idef->il[F_ORIRES].iatoms,
					idef->iparams,md,(const rvec*)x,
					pbc_null,fcd);
  if (idef->il[F_DISRES].nr)
    calc_disres_R_6(mcr,idef->il[F_DISRES].nr,
		    idef->il[F_DISRES].iatoms,
		    idef->iparams,(const rvec*)x,pbc_null,fcd);
  
  /* Loop over all bonded force types to calculate the bonded forces */
  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (interaction_function[ftype].flags & IF_BOND && ftype!=F_CONNBONDS) {
      nbonds=idef->il[ftype].nr;
      if (nbonds > 0) {
	ind = interaction_function[ftype].nrnb_ind;
	nat = interaction_function[ftype].nratoms+1;
	dvdl = 0;
	if (ftype != F_LJ14) {
	  v = interaction_function[ftype].ifunc
	    (nbonds,idef->il[ftype].iatoms,
	     idef->iparams,
	     (const rvec*)x,f,fr->fshift,
	     (ftype == F_POSRES) ? pbc : pbc_null,
	     g,lambda,&dvdl,md,fcd);
	  if (bSepDVDL) {
	    fprintf(fplog,"  %-23s #%4d  V %12.5e  dVdl %12.5e\n",
		    interaction_function[ftype].longname,nbonds/nat,v,dvdl);
	  }
	} else {
	  v = do_14(nbonds,idef->il[ftype].iatoms,
		    idef->iparams,
		    (const rvec*)x,f,fr->fshift,
		    pbc_null,g,
		    lambda,&dvdl,
		    md,fr,ngrp,egnb,egcoul);
	  if (bSepDVDL) {
	    fprintf(fplog,"  %-5s + %-15s #%4d                  dVdl %12.5e\n",
		    interaction_function[ftype].longname,
		    interaction_function[F_COUL14].longname,nbonds/nat,dvdl);
	  }
	}
	if (ind != -1)
	  inc_nrnb(nrnb,ind,nbonds/nat);
	epot[ftype]  += v;
	epot[F_DVDL] += dvdl;
      }
    }
  }
  /* Copy the sum of violations for the distance restraints from fcd */
  epot[F_DISRESVIOL] = fcd->disres.sumviol;
}

/*
 * Morse potential bond by Frank Everdij
 *
 * Three parameters needed:
 *
 * b0 = equilibrium distance in nm
 * be = beta in nm^-1 (actually, it's nu_e*Sqrt(2*pi*pi*mu/D_e))
 * cb = well depth in kJ/mol
 *
 * Note: the potential is referenced to be +cb at infinite separation
 *       and zero at the equilibrium distance!
 */

real morse_bonds(int nbonds,
		 const t_iatom forceatoms[],const t_iparams forceparams[],
		 const rvec x[],rvec f[],rvec fshift[],
		 const t_pbc *pbc,const t_graph *g,
		 real lambda,real *dvdl,
		 const t_mdatoms *md,t_fcdata *fcd)
{
  const real one=1.0;
  const real two=2.0;
  real  dr,dr2,temp,omtemp,cbomtemp,fbond,vbond,fij,b0,be,cb,vtot;
  rvec  dx;
  int   i,m,ki,type,ai,aj;
  ivec  dt;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    
    b0   = forceparams[type].morse.b0;
    be   = forceparams[type].morse.beta;
    cb   = forceparams[type].morse.cb;

    ki   = pbc_rvec_sub(pbc,x[ai],x[aj],dx);            /*   3          */
    dr2  = iprod(dx,dx);                            /*   5          */
    dr   = dr2*invsqrt(dr2);                        /*  10          */
    temp = exp(-be*(dr-b0));                        /*  12          */
    
    if (temp == one)
      continue;

    omtemp   = one-temp;                               /*   1          */
    cbomtemp = cb*omtemp;                              /*   1          */
    vbond    = cbomtemp*omtemp;                        /*   1          */
    fbond    = -two*be*temp*cbomtemp*invsqrt(dr2);      /*   9          */
    vtot    += vbond;       /* 1 */
    
    if (g) {
      ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
      ki = IVEC2IS(dt);
    }

    for (m=0; (m<DIM); m++) {                          /*  15          */
      fij=fbond*dx[m];
      f[ai][m]+=fij;
      f[aj][m]-=fij;
      fshift[ki][m]+=fij;
      fshift[CENTRAL][m]-=fij;
    }
  }                                           /*  58 TOTAL    */
  return vtot;
}

real cubic_bonds(int nbonds,
		 const t_iatom forceatoms[],const t_iparams forceparams[],
		 const rvec x[],rvec f[],rvec fshift[],
		 const t_pbc *pbc,const t_graph *g,
		 real lambda,real *dvdl,
		 const t_mdatoms *md,t_fcdata *fcd)
{
  const real three = 3.0;
  const real two   = 2.0;
  real  kb,b0,kcub;
  real  dr,dr2,dist,kdist,kdist2,fbond,vbond,fij,vtot;
  rvec  dx;
  int   i,m,ki,type,ai,aj;
  ivec  dt;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    
    b0   = forceparams[type].cubic.b0;
    kb   = forceparams[type].cubic.kb;
    kcub = forceparams[type].cubic.kcub;

    ki   = pbc_rvec_sub(pbc,x[ai],x[aj],dx);                /*   3          */
    dr2  = iprod(dx,dx);                                /*   5          */
    
    if (dr2 == 0.0)
      continue;
      
    dr         = dr2*invsqrt(dr2);                      /*  10          */
    dist       = dr-b0;
    kdist      = kb*dist;
    kdist2     = kdist*dist;
    
    vbond      = kdist2 + kcub*kdist2*dist;
    fbond      = -(two*kdist + three*kdist2*kcub)/dr;

    vtot      += vbond;       /* 21 */
    
    if (g) {
      ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
      ki=IVEC2IS(dt);
    }
    for (m=0; (m<DIM); m++) {                          /*  15          */
      fij=fbond*dx[m];
      f[ai][m]+=fij;
      f[aj][m]-=fij;
      fshift[ki][m]+=fij;
      fshift[CENTRAL][m]-=fij;
    }
  }                                           /*  54 TOTAL    */
  return vtot;
}

real FENE_bonds(int nbonds,
		const t_iatom forceatoms[],const t_iparams forceparams[],
		const rvec x[],rvec f[],rvec fshift[],
		const t_pbc *pbc,const t_graph *g,
		real lambda,real *dvdl,
		const t_mdatoms *md,t_fcdata *fcd)
{
  const real half=0.5;
  const real one=1.0;
  real  bm,kb;
  real  dr,dr2,bm2,omdr2obm2,fbond,vbond,fij,vtot;
  rvec  dx;
  int   i,m,ki,type,ai,aj;
  ivec  dt;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    
    bm   = forceparams[type].fene.bm;
    kb   = forceparams[type].fene.kb;

    ki   = pbc_rvec_sub(pbc,x[ai],x[aj],dx);            /*   3          */
    dr2  = iprod(dx,dx);                                /*   5          */
    
    if (dr2 == 0.0)
      continue;

    bm2 = bm*bm;

    if (dr2 >= bm2)
      gmx_fatal(FARGS,
		"r^2 (%f) >= bm^2 (%f) in FENE bond between atoms %d and %d",
		dr2,bm2,ai+1,aj+1);
      
    omdr2obm2  = one - dr2/bm2;
    
    vbond      = -half*bm2*log(omdr2obm2);
    fbond      = -bm/omdr2obm2;

    vtot      += vbond;       /* 35 */
    
    if (g) {
      ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
      ki=IVEC2IS(dt);
    }
    for (m=0; (m<DIM); m++) {                          /*  15          */
      fij=fbond*dx[m];
      f[ai][m]+=fij;
      f[aj][m]-=fij;
      fshift[ki][m]+=fij;
      fshift[CENTRAL][m]-=fij;
    }
  }                                           /*  58 TOTAL    */
  return vtot;
}

real harmonic(real kA,real kB,real xA,real xB,real x,real lambda,
	      real *V,real *F)
{
  const real half=0.5;
  real  L1,kk,x0,dx,dx2;
  real  v,f,dvdl;
  
  L1    = 1.0-lambda;
  kk    = L1*kA+lambda*kB;
  x0    = L1*xA+lambda*xB;
  
  dx    = x-x0;
  dx2   = dx*dx;
  
  f     = -kk*dx;
  v     = half*kk*dx2;
  dvdl  = half*(kB-kA)*dx2 + (xA-xB)*kk*dx;
  
  *F    = f;
  *V    = v;
  
  return dvdl;
  
  /* That was 19 flops */
}


real bonds(int nbonds,
	   const t_iatom forceatoms[],const t_iparams forceparams[],
	   const rvec x[],rvec f[],rvec fshift[],
	   const t_pbc *pbc,const t_graph *g,
	   real lambda,real *dvdlambda,
	   const t_mdatoms *md,t_fcdata *fcd)
{
  int  i,m,ki,ai,aj,type;
  real dr,dr2,fbond,vbond,fij,vtot;
  rvec dx;
  ivec dt;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
  
    ki   = pbc_rvec_sub(pbc,x[ai],x[aj],dx);	/*   3 		*/
    dr2  = iprod(dx,dx);			/*   5		*/
    dr   = dr2*invsqrt(dr2);		        /*  10		*/

    *dvdlambda += harmonic(forceparams[type].harmonic.krA,
			   forceparams[type].harmonic.krB,
			   forceparams[type].harmonic.rA,
			   forceparams[type].harmonic.rB,
			   dr,lambda,&vbond,&fbond);  /*  19  */

    if (dr2 == 0.0)
      continue;

    
    vtot  += vbond;/* 1*/
    fbond *= invsqrt(dr2);			/*   6		*/
#ifdef DEBUG
    if (debug)
      fprintf(debug,"BONDS: dr = %10g  vbond = %10g  fbond = %10g\n",
	      dr,vbond,fbond);
#endif
    if (g) {
      ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
      ki=IVEC2IS(dt);
    }
    for (m=0; (m<DIM); m++) {			/*  15		*/
      fij=fbond*dx[m];
      f[ai][m]+=fij;
      f[aj][m]-=fij;
      fshift[ki][m]+=fij;
      fshift[CENTRAL][m]-=fij;
    }
  }					/* 59 TOTAL	*/
  return vtot;
}

real polarize(int nbonds,
	      const t_iatom forceatoms[],const t_iparams forceparams[],
	      const rvec x[],rvec f[],rvec fshift[],
	      const t_pbc *pbc,const t_graph *g,
	      real lambda,real *dvdlambda,
	      const t_mdatoms *md,t_fcdata *fcd)
{
  int  i,m,ki,ai,aj,type;
  real dr,dr2,fbond,vbond,fij,vtot,ksh;
  rvec dx;
  ivec dt;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ksh  = sqr(md->chargeA[aj])*ONE_4PI_EPS0/forceparams[type].polarize.alpha;
  
    ki   = pbc_rvec_sub(pbc,x[ai],x[aj],dx);	/*   3 		*/
    dr2  = iprod(dx,dx);			/*   5		*/
    dr   = dr2*invsqrt(dr2);		        /*  10		*/

    *dvdlambda += harmonic(ksh,ksh,0,0,dr,lambda,&vbond,&fbond);  /*  19  */

    if (dr2 == 0.0)
      continue;
    
    vtot  += vbond;/* 1*/
    fbond *= invsqrt(dr2);			/*   6		*/

    if (g) {
      ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
      ki=IVEC2IS(dt);
    }
    for (m=0; (m<DIM); m++) {			/*  15		*/
      fij=fbond*dx[m];
      f[ai][m]+=fij;
      f[aj][m]-=fij;
      fshift[ki][m]+=fij;
      fshift[CENTRAL][m]-=fij;
    }
  }					/* 59 TOTAL	*/
  return vtot;
}

real water_pol(int nbonds,
	       const t_iatom forceatoms[],const t_iparams forceparams[],
	       const rvec x[],rvec f[],rvec fshift[],
	       const t_pbc *pbc,const t_graph *g,
	       real lambda,real *dvdlambda,
	       const t_mdatoms *md,t_fcdata *fcd)
{
  /* This routine implements anisotropic polarizibility for water, through
   * a shell connected to a dummy with spring constant that differ in the
   * three spatial dimensions in the molecular frame.
   */
  int  i,m,aO,aH1,aH2,aD,aS,type,type0;
  rvec dOH1,dOH2,dHH,dOD,dDS,nW,kk,dx,kdx,proj;
#ifdef DEBUG
  rvec df;
#endif
  real vtot,fij,r_HH,r_OD,r_nW,tx,ty,tz,qS;

  vtot = 0.0;
  if (nbonds > 0) {
    type0  = forceatoms[0];
    aS     = forceatoms[5];
    qS     = md->chargeA[aS];
    kk[XX] = sqr(qS)*ONE_4PI_EPS0/forceparams[type0].wpol.al_x;
    kk[YY] = sqr(qS)*ONE_4PI_EPS0/forceparams[type0].wpol.al_y;
    kk[ZZ] = sqr(qS)*ONE_4PI_EPS0/forceparams[type0].wpol.al_z;
    r_HH   = 1.0/forceparams[type0].wpol.rHH;
    r_OD   = 1.0/forceparams[type0].wpol.rOD;
    if (debug) {
      fprintf(debug,"WPOL: qS  = %10.5f aS = %5d\n",qS,aS);
      fprintf(debug,"WPOL: kk  = %10.3f        %10.3f        %10.3f\n",
	      kk[XX],kk[YY],kk[ZZ]);
      fprintf(debug,"WPOL: rOH = %10.3f  rHH = %10.3f  rOD = %10.3f\n",
	      forceparams[type0].wpol.rOH,
	      forceparams[type0].wpol.rHH,
	      forceparams[type0].wpol.rOD);
    }
    for(i=0; (i<nbonds); i+=6) {
      type = forceatoms[i];
      if (type != type0)
	gmx_fatal(FARGS,"Sorry, type = %d, type0 = %d, file = %s, line = %d",
		    type,type0,__FILE__,__LINE__);
      aO   = forceatoms[i+1];
      aH1  = forceatoms[i+2];
      aH2  = forceatoms[i+3];
      aD   = forceatoms[i+4];
      aS   = forceatoms[i+5];
      
      /* Compute vectors describing the water frame */
      rvec_sub(x[aH1],x[aO], dOH1);
      rvec_sub(x[aH2],x[aO], dOH2);
      rvec_sub(x[aH2],x[aH1],dHH);
      rvec_sub(x[aD], x[aO], dOD);
      rvec_sub(x[aS], x[aD], dDS);
      oprod(dOH1,dOH2,nW);
      
      /* Compute inverse length of normal vector 
       * (this one could be precomputed, but I'm too lazy now)
       */
      r_nW = invsqrt(iprod(nW,nW));
      /* This is for precision, but does not make a big difference,
       * it can go later.
       */
      r_OD = invsqrt(iprod(dOD,dOD)); 
      
      /* Normalize the vectors in the water frame */
      svmul(r_nW,nW,nW);
      svmul(r_HH,dHH,dHH);
      svmul(r_OD,dOD,dOD);
      
      /* Compute displacement of shell along components of the vector */
      dx[ZZ] = iprod(dDS,dOD);
      /* Compute projection on the XY plane: dDS - dx[ZZ]*dOD */
      for(m=0; (m<DIM); m++)
	proj[m] = dDS[m]-dx[ZZ]*dOD[m];
      
      /*dx[XX] = iprod(dDS,nW);
	dx[YY] = iprod(dDS,dHH);*/
      dx[XX] = iprod(proj,nW);
      for(m=0; (m<DIM); m++)
	proj[m] -= dx[XX]*nW[m];
      dx[YY] = iprod(proj,dHH);
      /*#define DEBUG*/
#ifdef DEBUG
      if (debug) {
	fprintf(debug,"WPOL: dx2=%10g  dy2=%10g  dz2=%10g  sum=%10g  dDS^2=%10g\n",
		sqr(dx[XX]),sqr(dx[YY]),sqr(dx[ZZ]),iprod(dx,dx),iprod(dDS,dDS));
	fprintf(debug,"WPOL: dHH=(%10g,%10g,%10g)\n",dHH[XX],dHH[YY],dHH[ZZ]);
	fprintf(debug,"WPOL: dOD=(%10g,%10g,%10g), 1/r_OD = %10g\n",
		dOD[XX],dOD[YY],dOD[ZZ],1/r_OD);
	fprintf(debug,"WPOL: nW =(%10g,%10g,%10g), 1/r_nW = %10g\n",
		nW[XX],nW[YY],nW[ZZ],1/r_nW);
	fprintf(debug,"WPOL: dx  =%10g, dy  =%10g, dz  =%10g\n",
		dx[XX],dx[YY],dx[ZZ]);
 	fprintf(debug,"WPOL: dDSx=%10g, dDSy=%10g, dDSz=%10g\n",
		dDS[XX],dDS[YY],dDS[ZZ]);
      }
#endif
      /* Now compute the forces and energy */
      kdx[XX] = kk[XX]*dx[XX];
      kdx[YY] = kk[YY]*dx[YY];
      kdx[ZZ] = kk[ZZ]*dx[ZZ];
      vtot   += iprod(dx,kdx);
      for(m=0; (m<DIM); m++) {
	/* This is a tensor operation but written out for speed */
	tx        =  nW[m]*kdx[XX];
	ty        = dHH[m]*kdx[YY];
	tz        = dOD[m]*kdx[ZZ];
	fij       = -tx-ty-tz;
#ifdef DEBUG
	df[m] = fij;
#endif
	f[aS][m] += fij;
	f[aD][m] -= fij;
      }
#ifdef DEBUG
      if (debug) {
	fprintf(debug,"WPOL: vwpol=%g\n",0.5*iprod(dx,kdx));
	fprintf(debug,"WPOL: df = (%10g, %10g, %10g)\n",df[XX],df[YY],df[ZZ]);
      }
#endif
    }	
  }
  return 0.5*vtot;
}

real bond_angle(const rvec xi,const rvec xj,const rvec xk,const t_pbc *pbc,
		rvec r_ij,rvec r_kj,real *costh,
		int *t1,int *t2)
/* Return value is the angle between the bonds i-j and j-k */
{
  /* 41 FLOPS */
  real th;
  
  *t1 = pbc_rvec_sub(pbc,xi,xj,r_ij);			/*  3		*/
  *t2 = pbc_rvec_sub(pbc,xk,xj,r_kj);			/*  3		*/

  *costh=cos_angle(r_ij,r_kj);		/* 25		*/
  th=acos(*costh);			/* 10		*/
					/* 41 TOTAL	*/
  return th;
}

real angles(int nbonds,
	    const t_iatom forceatoms[],const t_iparams forceparams[],
	    const rvec x[],rvec f[],rvec fshift[],
	    const t_pbc *pbc,const t_graph *g,
	    real lambda,real *dvdlambda,
	    const t_mdatoms *md,t_fcdata *fcd)
{
  int  i,ai,aj,ak,t1,t2,type;
  rvec r_ij,r_kj;
  real cos_theta,theta,dVdt,va,vtot;
  ivec jt,dt_ij,dt_kj;
  
  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    
    theta  = bond_angle(x[ai],x[aj],x[ak],pbc,
			r_ij,r_kj,&cos_theta,&t1,&t2);	/*  41		*/
  
    *dvdlambda += harmonic(forceparams[type].harmonic.krA,
			   forceparams[type].harmonic.krB,
			   forceparams[type].harmonic.rA*DEG2RAD,
			   forceparams[type].harmonic.rB*DEG2RAD,
			   theta,lambda,&va,&dVdt);  /*  21  */
    vtot += va;
    
    {
      int  m;
      real snt,st,sth;
      real cik,cii,ckk;
      real nrkj2,nrij2;
      rvec f_i,f_j,f_k;
      
      snt=sin(theta);				/*  10		*/
      if (fabs(snt) < 1e-12)
	snt=1e-12;
      st  = dVdt/snt;				/*  10		*/
      sth = st*cos_theta;			/*   1		*/
#ifdef DEBUG
      if (debug)
	fprintf(debug,"ANGLES: theta = %10g  vth = %10g  dV/dtheta = %10g\n",
		theta*RAD2DEG,va,dVdt);
#endif
      nrkj2=iprod(r_kj,r_kj);			/*   5		*/
      nrij2=iprod(r_ij,r_ij);
      
      cik=st*invsqrt(nrkj2*nrij2);		/*  12		*/ 
      cii=sth/nrij2;				/*  10		*/
      ckk=sth/nrkj2;				/*  10		*/
      
      for (m=0; (m<DIM); m++) {			/*  39		*/
	f_i[m]=-(cik*r_kj[m]-cii*r_ij[m]);
	f_k[m]=-(cik*r_ij[m]-ckk*r_kj[m]);
	f_j[m]=-f_i[m]-f_k[m];
	f[ai][m]+=f_i[m];
	f[aj][m]+=f_j[m];
	f[ak][m]+=f_k[m];
      }
      if (g) {
	copy_ivec(SHIFT_IVEC(g,aj),jt);
      
	ivec_sub(SHIFT_IVEC(g,ai),jt,dt_ij);
	ivec_sub(SHIFT_IVEC(g,ak),jt,dt_kj);
	t1=IVEC2IS(dt_ij);
	t2=IVEC2IS(dt_kj);
      }
      rvec_inc(fshift[t1],f_i);
      rvec_inc(fshift[CENTRAL],f_j);
      rvec_inc(fshift[t2],f_k);
    }                                           /* 168 TOTAL	*/
  }
  return vtot;
}

real urey_bradley(int nbonds,
		  const t_iatom forceatoms[],const t_iparams forceparams[],
		  const rvec x[],rvec f[],rvec fshift[],
		  const t_pbc *pbc,const t_graph *g,
		  real lambda,real *dvdlambda,
		  const t_mdatoms *md,t_fcdata *fcd)
{
  int  i,m,ai,aj,ak,t1,t2,type,ki;
  rvec r_ij,r_kj,r_ik;
  real cos_theta,theta,dVdt,va,vtot,kth,th0,kUB,r13,dr,dr2,vbond,fbond,fik;
  ivec jt,dt_ij,dt_kj,dt_ik;
  
  gmx_impl("Not implemented yet");
  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    th0  = forceparams[type].u_b.theta*DEG2RAD;
    kth  = forceparams[type].u_b.ktheta;
    r13  = forceparams[type].u_b.r13;
    kUB  = forceparams[type].u_b.kUB;
    
    theta  = bond_angle(x[ai],x[aj],x[ak],pbc,
			r_ij,r_kj,&cos_theta,&t1,&t2);	/*  41		*/
  
    *dvdlambda += harmonic(kth,kth,th0,th0,theta,lambda,&va,&dVdt);  /*  21  */
    vtot += va;
    
    ki   = pbc_rvec_sub(pbc,x[ai],x[ak],r_ik);	/*   3 		*/
    dr2  = iprod(r_ik,r_ik);			/*   5		*/
    dr   = dr2*invsqrt(dr2);		        /*  10		*/

    *dvdlambda += harmonic(kUB,kUB,r13,r13,dr,lambda,&vbond,&fbond); /*  19  */

    {
      real snt,st,sth;
      real cik,cii,ckk;
      real nrkj2,nrij2;
      rvec f_i,f_j,f_k;
      
      snt=sin(theta);				/*  10		*/
      if (fabs(snt) < 1e-12)
	snt=1e-12;
      st  = dVdt/snt;				/*  10		*/
      sth = st*cos_theta;			/*   1		*/
#ifdef DEBUG
      if (debug)
	fprintf(debug,"ANGLES: theta = %10g  vth = %10g  dV/dtheta = %10g\n",
		theta*RAD2DEG,va,dVdt);
#endif
      nrkj2=iprod(r_kj,r_kj);			/*   5		*/
      nrij2=iprod(r_ij,r_ij);
      
      cik=st*invsqrt(nrkj2*nrij2);		/*  12		*/ 
      cii=sth/nrij2;				/*  10		*/
      ckk=sth/nrkj2;				/*  10		*/
      
      for (m=0; (m<DIM); m++) {			/*  39		*/
	f_i[m]=-(cik*r_kj[m]-cii*r_ij[m]);
	f_k[m]=-(cik*r_ij[m]-ckk*r_kj[m]);
	f_j[m]=-f_i[m]-f_k[m];
	f[ai][m]+=f_i[m];
	f[aj][m]+=f_j[m];
	f[ak][m]+=f_k[m];
      }
      if (g) {
	copy_ivec(SHIFT_IVEC(g,aj),jt);
      
	ivec_sub(SHIFT_IVEC(g,ai),jt,dt_ij);
	ivec_sub(SHIFT_IVEC(g,ak),jt,dt_kj);
	t1=IVEC2IS(dt_ij);
	t2=IVEC2IS(dt_kj);
      }
      rvec_inc(fshift[t1],f_i);
      rvec_inc(fshift[CENTRAL],f_j);
      rvec_inc(fshift[t2],f_k);
    }                                           /* 168 TOTAL	*/
    /* Time for the bond calculations */
    if (dr2 == 0.0)
      continue;

    vtot  += vbond;  /* 1*/
    fbond *= invsqrt(dr2);			/*   6		*/
    
    if (g) {
      ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,ak),dt_ik);
      ki=IVEC2IS(dt_ik);
    }
    for (m=0; (m<DIM); m++) {			/*  15		*/
      fik=fbond*r_ik[m];
      f[ai][m]+=fik;
      f[ak][m]-=fik;
      fshift[ki][m]+=fik;
      fshift[CENTRAL][m]-=fik;
    }
  }
  return vtot;
}

real quartic_angles(int nbonds,
		    const t_iatom forceatoms[],const t_iparams forceparams[],
		    const rvec x[],rvec f[],rvec fshift[],
		    const t_pbc *pbc,const t_graph *g,
		    real lambda,real *dvdlambda,
		    const t_mdatoms *md,t_fcdata *fcd)
{
  int  i,j,ai,aj,ak,t1,t2,type;
  rvec r_ij,r_kj;
  real cos_theta,theta,dt,dVdt,va,dtp,c,vtot;
  ivec jt,dt_ij,dt_kj;
  
  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];

    theta  = bond_angle(x[ai],x[aj],x[ak],pbc,
			r_ij,r_kj,&cos_theta,&t1,&t2);	/*  41		*/

    dt = theta - forceparams[type].qangle.theta*DEG2RAD; /* 2          */

    dVdt = 0;
    va = forceparams[type].qangle.c[0];
    dtp = 1.0;
    for(j=1; j<=4; j++) {
      c = forceparams[type].qangle.c[j];
      dVdt -= j*c*dtp;
      dtp *= dt;
      va += c*dtp;
    }
    /* 20 */

    vtot += va;
    
    {
      int  m;
      real snt,st,sth;
      real cik,cii,ckk;
      real nrkj2,nrij2;
      rvec f_i,f_j,f_k;
      
      snt=sin(theta);				/*  10		*/
      if (fabs(snt) < 1e-12)
	snt=1e-12;
      st  = dVdt/snt;				/*  10		*/
      sth = st*cos_theta;			/*   1		*/
#ifdef DEBUG
      if (debug)
	fprintf(debug,"ANGLES: theta = %10g  vth = %10g  dV/dtheta = %10g\n",
		theta*RAD2DEG,va,dVdt);
#endif
      nrkj2=iprod(r_kj,r_kj);			/*   5		*/
      nrij2=iprod(r_ij,r_ij);
      
      cik=st*invsqrt(nrkj2*nrij2);		/*  12		*/ 
      cii=sth/nrij2;				/*  10		*/
      ckk=sth/nrkj2;				/*  10		*/
      
      for (m=0; (m<DIM); m++) {			/*  39		*/
	f_i[m]=-(cik*r_kj[m]-cii*r_ij[m]);
	f_k[m]=-(cik*r_ij[m]-ckk*r_kj[m]);
	f_j[m]=-f_i[m]-f_k[m];
	f[ai][m]+=f_i[m];
	f[aj][m]+=f_j[m];
	f[ak][m]+=f_k[m];
      }
      if (g) {
	copy_ivec(SHIFT_IVEC(g,aj),jt);
      
	ivec_sub(SHIFT_IVEC(g,ai),jt,dt_ij);
	ivec_sub(SHIFT_IVEC(g,ak),jt,dt_kj);
	t1=IVEC2IS(dt_ij);
	t2=IVEC2IS(dt_kj);
      }
      rvec_inc(fshift[t1],f_i);
      rvec_inc(fshift[CENTRAL],f_j);
      rvec_inc(fshift[t2],f_k);
    }                                           /* 160 TOTAL	*/
  }
  return vtot;
}

real dih_angle(const rvec xi,const rvec xj,const rvec xk,const rvec xl,
	       const t_pbc *pbc,
	       rvec r_ij,rvec r_kj,rvec r_kl,rvec m,rvec n,
	       real *cos_phi,real *sign,int *t1,int *t2,int *t3)
{
  real ipr,phi;

  *t1 = pbc_rvec_sub(pbc,xi,xj,r_ij);       		/*  3 		*/
  *t2 = pbc_rvec_sub(pbc,xk,xj,r_kj);			/*  3		*/
  *t3 = pbc_rvec_sub(pbc,xk,xl,r_kl);			/*  3		*/

  oprod(r_ij,r_kj,m); 			/*  9 		*/
  oprod(r_kj,r_kl,n);			/*  9		*/
  *cos_phi=cos_angle(m,n); 		/* 41 		*/
  phi=acos(*cos_phi); 			/* 10 		*/
  ipr=iprod(r_ij,n); 			/*  5 		*/
  (*sign)=(ipr<0.0)?-1.0:1.0;
  phi=(*sign)*phi; 			/*  1		*/
					/* 84 TOTAL	*/
  return phi;
}

void do_dih_fup(int i,int j,int k,int l,real ddphi,
		rvec r_ij,rvec r_kj,rvec r_kl,
		rvec m,rvec n,rvec f[],rvec fshift[],
		const t_pbc *pbc,const t_graph *g,
		const rvec x[],int t1,int t2,int t3)
{
  /* 143 FLOPS */
  rvec f_i,f_j,f_k,f_l;
  rvec uvec,vvec,svec,dx_jl;
  real ipr,nrkj,nrkj2;
  real a,p,q;
  ivec jt,dt_ij,dt_kj,dt_lj;  
  
  ipr   = iprod(m,m);		/*  5 	*/
  nrkj2 = iprod(r_kj,r_kj);	/*  5	*/
  nrkj  = nrkj2*invsqrt(nrkj2);	/* 10	*/
  a     = -ddphi*nrkj/ipr;	/* 11	*/
  svmul(a,m,f_i);		/*  3	*/
  ipr   = iprod(n,n);		/*  5	*/
  a     = ddphi*nrkj/ipr;	/* 11	*/
  svmul(a,n,f_l);		/*  3 	*/
  p     = iprod(r_ij,r_kj);	/*  5	*/
  p    /= nrkj2;		/* 10	*/
  q     = iprod(r_kl,r_kj);	/*  5	*/
  q    /= nrkj2;		/* 10	*/
  svmul(p,f_i,uvec);		/*  3	*/
  svmul(q,f_l,vvec);		/*  3	*/
  rvec_sub(uvec,vvec,svec);	/*  3	*/
  rvec_sub(f_i,svec,f_j);	/*  3	*/
  rvec_add(f_l,svec,f_k);	/*  3	*/
  rvec_inc(f[i],f_i);   	/*  3	*/
  rvec_dec(f[j],f_j);   	/*  3	*/
  rvec_dec(f[k],f_k);   	/*  3	*/
  rvec_inc(f[l],f_l);   	/*  3	*/

  if (g) {
    copy_ivec(SHIFT_IVEC(g,j),jt);
    ivec_sub(SHIFT_IVEC(g,i),jt,dt_ij);
    ivec_sub(SHIFT_IVEC(g,k),jt,dt_kj);
    ivec_sub(SHIFT_IVEC(g,l),jt,dt_lj);
    t1=IVEC2IS(dt_ij);
    t2=IVEC2IS(dt_kj);
    t3=IVEC2IS(dt_lj);
  }    
  else
    t3 = pbc_rvec_sub(pbc,x[l],x[j],dx_jl);
    
  rvec_inc(fshift[t1],f_i);
  rvec_dec(fshift[CENTRAL],f_j);
  rvec_dec(fshift[t2],f_k);
  rvec_inc(fshift[t3],f_l);
  
  /* 112 TOTAL 	*/
}


real dopdihs(real cpA,real cpB,real phiA,real phiB,int mult,
	     real phi,real lambda,real *V,real *F)
{
  real v,dvdl,mdphi,v1,sdphi,ddphi;
  real L1   = 1.0-lambda;
  real ph0  = DEG2RAD*(L1*phiA+lambda*phiB);
  real cp   = L1*cpA + lambda*cpB;
  
  mdphi =  mult*phi-ph0;
  sdphi = sin(mdphi);
  ddphi = -cp*mult*sdphi;
  v1    = 1.0+cos(mdphi);
  v     = cp*v1;
  
  dvdl  = (cpB-cpA)*v1 - cp*(phiA-phiB)*sdphi;
  
  *V = v;
  *F = ddphi;
  
  return dvdl;
  
  /* That was 40 flops */
}

static real dopdihs_min(real cpA,real cpB,real phiA,real phiB,int mult,
		       real phi,real lambda,real *V,real *F)
     /* similar to dopdihs, except for a minus sign  *
      * and a different treatment of mult/phi0       */
{
  real v,dvdl,mdphi,v1,sdphi,ddphi;
  real L1   = 1.0-lambda;
  real ph0  = DEG2RAD*(L1*phiA+lambda*phiB);
  real cp   = L1*cpA + lambda*cpB;
  
  mdphi = mult*(phi-ph0);
  sdphi = sin(mdphi);
  ddphi = cp*mult*sdphi;
  v1    = 1.0-cos(mdphi);
  v     = cp*v1;
  
  dvdl  = (cpB-cpA)*v1 - cp*(phiA-phiB)*sdphi;
  
  *V = v;
  *F = ddphi;
  
  return dvdl;
  
  /* That was 40 flops */
}

real pdihs(int nbonds,
	   const t_iatom forceatoms[],const t_iparams forceparams[],
	   const rvec x[],rvec f[],rvec fshift[],
	   const t_pbc *pbc,const t_graph *g,
	   real lambda,real *dvdlambda,
	   const t_mdatoms *md,t_fcdata *fcd)
{
  int  i,type,ai,aj,ak,al;
  int  t1,t2,t3;
  rvec r_ij,r_kj,r_kl,m,n;
  real phi,cos_phi,sign,ddphi,vpd,vtot;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    al   = forceatoms[i++];
    
    phi=dih_angle(x[ai],x[aj],x[ak],x[al],pbc,r_ij,r_kj,r_kl,m,n,
		  &cos_phi,&sign,&t1,&t2,&t3);			/*  84 		*/
		
    *dvdlambda += dopdihs(forceparams[type].pdihs.cpA,
			  forceparams[type].pdihs.cpB,
			  forceparams[type].pdihs.phiA,
			  forceparams[type].pdihs.phiB,
			  forceparams[type].pdihs.mult,
			  phi,lambda,&vpd,&ddphi);
		       
    vtot += vpd;
    do_dih_fup(ai,aj,ak,al,ddphi,r_ij,r_kj,r_kl,m,n,
	       f,fshift,pbc,g,x,t1,t2,t3);			/* 112		*/

#ifdef DEBUG
    fprintf(debug,"pdih: (%d,%d,%d,%d) cp=%g, phi=%g\n",
	    ai,aj,ak,al,cos_phi,phi);
#endif
  } /* 223 TOTAL 	*/

  return vtot;
}



real idihs(int nbonds,
	   const t_iatom forceatoms[],const t_iparams forceparams[],
	   const rvec x[],rvec f[],rvec fshift[],
	   const t_pbc *pbc,const t_graph *g,
	   real lambda,real *dvdlambda,
	   const t_mdatoms *md,t_fcdata *fcd)
{
  int  i,type,ai,aj,ak,al;
  int  t1,t2,t3;
  real phi,phi0,cos_phi,ddphi,sign,vtot;
  rvec r_ij,r_kj,r_kl,m,n;
  real L1,kk,dp,dp2,kA,kB,pA,pB,dvdl;

  L1 = 1.0-lambda;
  dvdl = 0;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    al   = forceatoms[i++];
    
    phi=dih_angle(x[ai],x[aj],x[ak],x[al],pbc,r_ij,r_kj,r_kl,m,n,
		  &cos_phi,&sign,&t1,&t2,&t3);			/*  84		*/
    
    /* phi can jump if phi0 is close to Pi/-Pi, which will cause huge
     * force changes if we just apply a normal harmonic.
     * Instead, we first calculate phi-phi0 and take it modulo (-Pi,Pi).
     * This means we will never have the periodicity problem, unless
     * the dihedral is Pi away from phiO, which is very unlikely due to
     * the potential.
     */
    kA = forceparams[type].harmonic.krA;
    kB = forceparams[type].harmonic.krB;
    pA = forceparams[type].harmonic.rA;
    pB = forceparams[type].harmonic.rB;

    kk   = L1*kA+lambda*kB;
    phi0 = (L1*pA+lambda*pB)*DEG2RAD;
    
    /* dp = (phi-phi0), modulo (-pi,pi) */
    dp = phi-phi0;  
    /* dp cannot be outside (-2*pi,2*pi) */
    if (dp>=M_PI)
      dp -= 2*M_PI;
    else if(dp<-M_PI)
      dp += 2*M_PI;
    
    dp2 = dp*dp;

    vtot += 0.5*kk*dp2;
    ddphi = -kk*dp;
    
    dvdl += 0.5*(kB-kA)*dp2 + (pA-pB)*kk*dp;

    do_dih_fup(ai,aj,ak,al,(real)(-ddphi),r_ij,r_kj,r_kl,m,n,
	       f,fshift,pbc,g,x,t1,t2,t3);			/* 112		*/
    /* 217 TOTAL	*/
#ifdef DEBUG
    fprintf("idih: (%d,%d,%d,%d) cp=%g, phi=%g\n",
	    ai,aj,ak,al,cos_phi,phi);
#endif
  }
  
  *dvdlambda += dvdl;
  return vtot;
}


real posres(int nbonds,
	    const t_iatom forceatoms[],const t_iparams forceparams[],
	    const rvec x[],rvec f[],rvec fshift[],
	    const t_pbc *pbc,const t_graph *g,
	    real lambda,real *dvdlambda,
	    const t_mdatoms *md,t_fcdata *fcd)
{
  int  i,ai,m,type,ki;
  const t_iparams *pr;
  real v,vtot,fm,*fc;
  rvec dx;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    pr   = &forceparams[type];
    
    /*
    if (pbc == NULL)
      rvec_sub(x[ai],forceparams[type].posres.pos0A,dx);
    else
    pbc_dx(pbc,x[ai],forceparams[type].posres.pos0A,dx);*/
    
    ki = pbc_rvec_sub(pbc,x[ai],forceparams[type].posres.pos0A,dx);
    
    v=0;
    for (m=0; (m<DIM); m++) {
      *dvdlambda += harmonic(pr->posres.fcA[m],pr->posres.fcB[m],
			     0,pr->posres.pos0B[m]-pr->posres.pos0A[m],
			     dx[m],lambda,&v,&fm);
      vtot += v;
      f[ai][m] += fm;
      /* Position restraints should not contribute to the virial.
       * We need to subtract this contribution !!!
       * Needs to be implemented, but requires an extra virial matrix.
       */
    }
  }

  return vtot;
}

static real low_angres(int nbonds,
		       const t_iatom forceatoms[],const t_iparams forceparams[],
		       const rvec x[],rvec f[],rvec fshift[],
		       const t_pbc *pbc,const t_graph *g,
		       real lambda,real *dvdlambda,
		       bool bZAxis)
{
  int  i,m,type,ai,aj,ak,al;
  int  t1,t2;
  real phi,cos_phi,vid,vtot,dVdphi;
  rvec r_ij,r_kl,f_i,f_k;
  real st,sth,sin_phi,nrij2,nrkl2,c,cij,ckl;

  ivec dt;  
  t2 = 0; /* avoid warning with gcc-3.3. It is never used uninitialized */

  vtot = 0.0;
  ak=al=0; /* to avoid warnings */
  for(i=0; i<nbonds; ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    t1   = pbc_rvec_sub(pbc,x[aj],x[ai],r_ij);            	/*  3		*/
    if (!bZAxis) {      
      ak   = forceatoms[i++];
      al   = forceatoms[i++];
      t2   = pbc_rvec_sub(pbc,x[al],x[ak],r_kl);           /*  3		*/
    } else {
      r_kl[XX] = 0;
      r_kl[YY] = 0;
      r_kl[ZZ] = 1;
    }

    cos_phi = cos_angle(r_ij,r_kl);		/* 25		*/
    phi     = acos(cos_phi);                    /* 10           */

    *dvdlambda += dopdihs_min(forceparams[type].pdihs.cpA,
			      forceparams[type].pdihs.cpB,
			      forceparams[type].pdihs.phiA,
			      forceparams[type].pdihs.phiB,
			      forceparams[type].pdihs.mult,
			      phi,lambda,&vid,&dVdphi); /*  40  */
    
    vtot += vid;

    sin_phi = sin(phi);			        /*  10		*/
    if (fabs(sin_phi) < 1e-12)
      sin_phi=1e-12;
    st  = -dVdphi/sin_phi;		        /*  10		*/
    sth = st*cos_phi;				/*   1		*/
    nrij2 = iprod(r_ij,r_ij);			/*   5		*/
    nrkl2 = iprod(r_kl,r_kl);                   /*   5          */
    
    c   = st*invsqrt(nrij2*nrkl2);		/*  11		*/ 
    cij = sth/nrij2;				/*  10		*/
    ckl = sth/nrkl2;				/*  10		*/
    
    for (m=0; m<DIM; m++) {			/*  18+18       */
      f_i[m] = (c*r_kl[m]-cij*r_ij[m]);
      f[ai][m] += f_i[m];
      f[aj][m] -= f_i[m];
      if (!bZAxis) {
	f_k[m] = (c*r_ij[m]-ckl*r_kl[m]);
	f[ak][m] += f_k[m];
	f[al][m] -= f_k[m];
      }
    }
    
    if (g) {
      ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
      t1=IVEC2IS(dt);
    }
    rvec_inc(fshift[t1],f_i);
    rvec_dec(fshift[CENTRAL],f_i);
    if (!bZAxis) {
      if (g) {
	ivec_sub(SHIFT_IVEC(g,ak),SHIFT_IVEC(g,al),dt);
	t2=IVEC2IS(dt);
      }
      rvec_inc(fshift[t2],f_k);
      rvec_dec(fshift[CENTRAL],f_k);
    }
  }

  return vtot;  /*  191 / 164 (bZAxis)  total  */
}

real angres(int nbonds,
	    const t_iatom forceatoms[],const t_iparams forceparams[],
	    const rvec x[],rvec f[],rvec fshift[],
	    const t_pbc *pbc,const t_graph *g,
	    real lambda,real *dvdlambda,
	    const t_mdatoms *md,t_fcdata *fcd)
{
  return low_angres(nbonds,forceatoms,forceparams,x,f,fshift,pbc,g,
		    lambda,dvdlambda,FALSE);
}

real angresz(int nbonds,
	     const t_iatom forceatoms[],const t_iparams forceparams[],
	     const rvec x[],rvec f[],rvec fshift[],
	     const t_pbc *pbc,const t_graph *g,
	     real lambda,real *dvdlambda,
	     const t_mdatoms *md,t_fcdata *fcd)
{
  return low_angres(nbonds,forceatoms,forceparams,x,f,fshift,pbc,g,
		    lambda,dvdlambda,TRUE);
}


real unimplemented(int nbonds,
		   const t_iatom forceatoms[],const t_iparams forceparams[],
		   const rvec x[],rvec f[],rvec fshift[],
		   const t_pbc *pbc,const t_graph *g,
		   real lambda,real *dvdlambda,
		   const t_mdatoms *md,t_fcdata *fcd)
{
  gmx_impl("*** you are using a not implemented function");

  return 0.0; /* To make the compiler happy */
}

real rbdihs(int nbonds,
	    const t_iatom forceatoms[],const t_iparams forceparams[],
	    const rvec x[],rvec f[],rvec fshift[],
	    const t_pbc *pbc,const t_graph *g,
	    real lambda,real *dvdlambda,
	    const t_mdatoms *md,t_fcdata *fcd)
{
  const real c0=0.0,c1=1.0,c2=2.0,c3=3.0,c4=4.0,c5=5.0;
  int  type,ai,aj,ak,al,i,j;
  int  t1,t2,t3;
  rvec r_ij,r_kj,r_kl,m,n;
  real parmA[NR_RBDIHS];
  real parmB[NR_RBDIHS];
  real parm[NR_RBDIHS];
  real phi,cos_phi,rbp,rbpBA;
  real v,sign,ddphi,sin_phi;
  real cosfac,vtot;
  real L1   = 1.0-lambda;
  real dvdl=0;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    al   = forceatoms[i++];

    phi=dih_angle(x[ai],x[aj],x[ak],x[al],pbc,r_ij,r_kj,r_kl,m,n,
		  &cos_phi,&sign,&t1,&t2,&t3);			/*  84		*/

    /* Change to polymer convention */
    if (phi < c0)
      phi += M_PI;
    else
      phi -= M_PI;			/*   1		*/
    cos_phi = -cos_phi;					/*   1		*/
    
    sin_phi=sin(phi);

    for(j=0; (j<NR_RBDIHS); j++) {
      parmA[j] = forceparams[type].rbdihs.rbcA[j];
      parmB[j] = forceparams[type].rbdihs.rbcB[j];
      parm[j]  = L1*parmA[j]+lambda*parmB[j];
    }
    /* Calculate cosine powers */
    /* Calculate the energy */
    /* Calculate the derivative */

    v       = parm[0];
    dvdl   += (parmB[0]-parmA[0]);
    ddphi   = c0;
    cosfac  = c1;
    
    rbp     = parm[1];
    rbpBA   = parmB[1]-parmA[1];
    ddphi  += rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    dvdl   += cosfac*rbpBA;
    rbp     = parm[2];
    rbpBA   = parmB[2]-parmA[2];    
    ddphi  += c2*rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    dvdl   += cosfac*rbpBA;
    rbp     = parm[3];
    rbpBA   = parmB[3]-parmA[3];
    ddphi  += c3*rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    dvdl   += cosfac*rbpBA;
    rbp     = parm[4];
    rbpBA   = parmB[4]-parmA[4];
    ddphi  += c4*rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    dvdl   += cosfac*rbpBA;
    rbp     = parm[5];
    rbpBA   = parmB[5]-parmA[5];
    ddphi  += c5*rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    dvdl   += cosfac*rbpBA;
   
    ddphi = -ddphi*sin_phi;				/*  11		*/
    
    do_dih_fup(ai,aj,ak,al,ddphi,r_ij,r_kj,r_kl,m,n,
	       f,fshift,pbc,g,x,t1,t2,t3);		/* 112		*/
    vtot += v;
  }  
  *dvdlambda += dvdl;

  return vtot;
}

/***********************************************************
 *
 *   G R O M O S  9 6   F U N C T I O N S
 *
 ***********************************************************/
real g96harmonic(real kA,real kB,real xA,real xB,real x,real lambda,
		 real *V,real *F)
{
  const real half=0.5;
  real  L1,kk,x0,dx,dx2;
  real  v,f,dvdl;
  
  L1    = 1.0-lambda;
  kk    = L1*kA+lambda*kB;
  x0    = L1*xA+lambda*xB;
  
  dx    = x-x0;
  dx2   = dx*dx;
  
  f     = -kk*dx;
  v     = half*kk*dx2;
  dvdl  = half*(kB-kA)*dx2 + (xA-xB)*kk*dx;
  
  *F    = f;
  *V    = v;
  
  return dvdl;
  
  /* That was 21 flops */
}

real g96bonds(int nbonds,
	      const t_iatom forceatoms[],const t_iparams forceparams[],
	      const rvec x[],rvec f[],rvec fshift[],
	      const t_pbc *pbc,const t_graph *g,
	      real lambda,real *dvdlambda,
	      const t_mdatoms *md,t_fcdata *fcd)
{
  int  i,m,ki,ai,aj,type;
  real dr2,fbond,vbond,fij,vtot;
  rvec dx;
  ivec dt;
  
  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
  
    ki   = pbc_rvec_sub(pbc,x[ai],x[aj],dx);		/*   3 		*/
    dr2  = iprod(dx,dx);				/*   5		*/
      
    *dvdlambda += g96harmonic(forceparams[type].harmonic.krA,
			      forceparams[type].harmonic.krB,
			      forceparams[type].harmonic.rA,
			      forceparams[type].harmonic.rB,
			      dr2,lambda,&vbond,&fbond);

    vtot  += 0.5*vbond;                             /* 1*/
#ifdef DEBUG
    if (debug)
      fprintf(debug,"G96-BONDS: dr = %10g  vbond = %10g  fbond = %10g\n",
	      sqrt(dr2),vbond,fbond);
#endif
   
    if (g) {
      ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
      ki=IVEC2IS(dt);
    }
    for (m=0; (m<DIM); m++) {			/*  15		*/
      fij=fbond*dx[m];
      f[ai][m]+=fij;
      f[aj][m]-=fij;
      fshift[ki][m]+=fij;
      fshift[CENTRAL][m]-=fij;
    }
  }					/* 44 TOTAL	*/
  return vtot;
}

real g96bond_angle(const rvec xi,const rvec xj,const rvec xk,const t_pbc *pbc,
		   rvec r_ij,rvec r_kj,
		   int *t1,int *t2)
/* Return value is the angle between the bonds i-j and j-k */
{
  real costh;
  
  *t1 = pbc_rvec_sub(pbc,xi,xj,r_ij);			/*  3		*/
  *t2 = pbc_rvec_sub(pbc,xk,xj,r_kj);			/*  3		*/

  costh=cos_angle(r_ij,r_kj);		/* 25		*/
					/* 41 TOTAL	*/
  return costh;
}

real g96angles(int nbonds,
	       const t_iatom forceatoms[],const t_iparams forceparams[],
	       const rvec x[],rvec f[],rvec fshift[],
	       const t_pbc *pbc,const t_graph *g,
	       real lambda,real *dvdlambda,
	       const t_mdatoms *md,t_fcdata *fcd)
{
  int  i,ai,aj,ak,type,m,t1,t2;
  rvec r_ij,r_kj;
  real cos_theta,dVdt,va,vtot;
  real rij_1,rij_2,rkj_1,rkj_2,rijrkj_1;
  rvec f_i,f_j,f_k;
  ivec jt,dt_ij,dt_kj;
  
  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    
    cos_theta  = g96bond_angle(x[ai],x[aj],x[ak],pbc,r_ij,r_kj,&t1,&t2);

    *dvdlambda += g96harmonic(forceparams[type].harmonic.krA,
			      forceparams[type].harmonic.krB,
			      forceparams[type].harmonic.rA,
			      forceparams[type].harmonic.rB,
			      cos_theta,lambda,&va,&dVdt);
    vtot    += va;
    
    rij_1    = invsqrt(iprod(r_ij,r_ij));
    rkj_1    = invsqrt(iprod(r_kj,r_kj));
    rij_2    = rij_1*rij_1;
    rkj_2    = rkj_1*rkj_1;
    rijrkj_1 = rij_1*rkj_1;                     /* 23 */
    
#ifdef DEBUG
    if (debug)
      fprintf(debug,"G96ANGLES: costheta = %10g  vth = %10g  dV/dct = %10g\n",
	      cos_theta,va,dVdt);
#endif
    for (m=0; (m<DIM); m++) {			/*  42	*/
      f_i[m]=dVdt*(r_kj[m]*rijrkj_1 - r_ij[m]*rij_2*cos_theta);
      f_k[m]=dVdt*(r_ij[m]*rijrkj_1 - r_kj[m]*rkj_2*cos_theta);
      f_j[m]=-f_i[m]-f_k[m];
      f[ai][m]+=f_i[m];
      f[aj][m]+=f_j[m];
      f[ak][m]+=f_k[m];
    }
    
    if (g) {
      copy_ivec(SHIFT_IVEC(g,aj),jt);
      
      ivec_sub(SHIFT_IVEC(g,ai),jt,dt_ij);
      ivec_sub(SHIFT_IVEC(g,ak),jt,dt_kj);
      t1=IVEC2IS(dt_ij);
      t2=IVEC2IS(dt_kj);
    }      
    rvec_inc(fshift[t1],f_i);
    rvec_inc(fshift[CENTRAL],f_j);
    rvec_inc(fshift[t2],f_k);               /* 9 */
    /* 163 TOTAL	*/
  }
  return vtot;
}

real cross_bond_bond(int nbonds,
		     const t_iatom forceatoms[],const t_iparams forceparams[],
		     const rvec x[],rvec f[],rvec fshift[],
		     const t_pbc *pbc,const t_graph *g,
		     real lambda,real *dvdlambda,
		     const t_mdatoms *md,t_fcdata *fcd)
{
  /* Potential from Lawrence and Skimmer, Chem. Phys. Lett. 372 (2003)
   * pp. 842-847
   */
  int  i,ai,aj,ak,type,m,t1,t2;
  rvec r_ij,r_kj;
  real vtot,vrr,s1,s2,r1,r2,r1e,r2e,krr;
  rvec f_i,f_j,f_k;
  ivec jt,dt_ij,dt_kj;
  
  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    r1e  = forceparams[type].cross_bb.r1e;
    r2e  = forceparams[type].cross_bb.r2e;
    krr  = forceparams[type].cross_bb.krr;
    
    /* Compute distance vectors ... */
    t1 = pbc_rvec_sub(pbc,x[ai],x[aj],r_ij);
    t2 = pbc_rvec_sub(pbc,x[ak],x[aj],r_kj);
    
    /* ... and their lengths */
    r1 = norm(r_ij);
    r2 = norm(r_kj);
    
    /* Deviations from ideality */
    s1 = r1-r1e;
    s2 = r2-r2e;
    
    /* Energy (can be negative!) */
    vrr   = krr*s1*s2;
    vtot += vrr;
    
    /* Forces */
    svmul(-krr*s2/r1,r_ij,f_i);
    svmul(-krr*s1/r2,r_kj,f_k);
    
    for (m=0; (m<DIM); m++) {			/*  12	*/
      f_j[m]    = -f_i[m] - f_k[m];
      f[ai][m] += f_i[m];
      f[aj][m] += f_j[m];
      f[ak][m] += f_k[m];
    }
    
    /* Virial stuff */
    if (g) {
      copy_ivec(SHIFT_IVEC(g,aj),jt);
      
      ivec_sub(SHIFT_IVEC(g,ai),jt,dt_ij);
      ivec_sub(SHIFT_IVEC(g,ak),jt,dt_kj);
      t1=IVEC2IS(dt_ij);
      t2=IVEC2IS(dt_kj);
    }      
    rvec_inc(fshift[t1],f_i);
    rvec_inc(fshift[CENTRAL],f_j);
    rvec_inc(fshift[t2],f_k);               /* 9 */
    /* 163 TOTAL	*/
  }
  return vtot;
}

real cross_bond_angle(int nbonds,
		      const t_iatom forceatoms[],const t_iparams forceparams[],
		      const rvec x[],rvec f[],rvec fshift[],
		      const t_pbc *pbc,const t_graph *g,
		      real lambda,real *dvdlambda,
		      const t_mdatoms *md,t_fcdata *fcd)
{
  /* Potential from Lawrence and Skimmer, Chem. Phys. Lett. 372 (2003)
   * pp. 842-847
   */
  int  i,ai,aj,ak,type,m,t1,t2,t3;
  rvec r_ij,r_kj,r_ik;
  real vtot,vrt,s1,s2,s3,r1,r2,r3,r1e,r2e,r3e,krt,k1,k2,k3;
  rvec f_i,f_j,f_k;
  ivec jt,dt_ij,dt_kj;
  
  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    r1e  = forceparams[type].cross_ba.r1e;
    r2e  = forceparams[type].cross_ba.r2e;
    r3e  = forceparams[type].cross_ba.r3e;
    krt  = forceparams[type].cross_ba.krt;
    
    /* Compute distance vectors ... */
    t1 = pbc_rvec_sub(pbc,x[ai],x[aj],r_ij);
    t2 = pbc_rvec_sub(pbc,x[ak],x[aj],r_kj);
    t3 = pbc_rvec_sub(pbc,x[ai],x[ak],r_ik);
    
    /* ... and their lengths */
    r1 = norm(r_ij);
    r2 = norm(r_kj);
    r3 = norm(r_ik);
    
    /* Deviations from ideality */
    s1 = r1-r1e;
    s2 = r2-r2e;
    s3 = r3-r3e;
    
    /* Energy (can be negative!) */
    vrt   = krt*s3*(s1+s2);
    vtot += vrt;
    
    /* Forces */
    k1 = -krt*(s3/r1);
    k2 = -krt*(s3/r2);
    k3 = -krt*(s1+s2)/r3;
    for(m=0; (m<DIM); m++) {
      f_i[m] = k1*r_ij[m] + k3*r_ik[m];
      f_k[m] = k2*r_kj[m] - k3*r_ik[m];
      f_j[m] = -f_i[m] - f_k[m];
    }
    
    for (m=0; (m<DIM); m++) {			/*  12	*/
      f[ai][m] += f_i[m];
      f[aj][m] += f_j[m];
      f[ak][m] += f_k[m];
    }
    
    /* Virial stuff */
    if (g) {
      copy_ivec(SHIFT_IVEC(g,aj),jt);
      
      ivec_sub(SHIFT_IVEC(g,ai),jt,dt_ij);
      ivec_sub(SHIFT_IVEC(g,ak),jt,dt_kj);
      t1=IVEC2IS(dt_ij);
      t2=IVEC2IS(dt_kj);
    }      
    rvec_inc(fshift[t1],f_i);
    rvec_inc(fshift[CENTRAL],f_j);
    rvec_inc(fshift[t2],f_k);               /* 9 */
    /* 163 TOTAL	*/
  }
  return vtot;
}

