/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
#include "gmx_fatal.h"
#include "mshift.h"
#include "main.h"
#include "disre.h"
#include "orires.h"
#include "force.h"
#include "nonbonded.h"
#include "mdrun.h"

/* Find a better place for this? */
const int cmap_coeff_matrix[] = {
1, 0, -3,  2, 0, 0,  0,  0, -3,  0,  9, -6,  2,  0, -6,  4 ,
0, 0,  0,  0, 0, 0,  0,  0,  3,  0, -9,  6, -2,  0,  6, -4,
0, 0,  0,  0, 0, 0,  0,  0,  0,  0,  9, -6,  0,  0, -6,  4 ,
0, 0,  3, -2, 0, 0,  0,  0,  0,  0, -9,  6,  0,  0,  6, -4,
0, 0,  0,  0, 1, 0, -3,  2, -2,  0,  6, -4,  1,  0, -3,  2 ,
0, 0,  0,  0, 0, 0,  0,  0, -1,  0,  3, -2,  1,  0, -3,  2 ,
0, 0,  0,  0, 0, 0,  0,  0,  0,  0, -3,  2,  0,  0,  3, -2,
0, 0,  0,  0, 0, 0,  3, -2,  0,  0, -6,  4,  0,  0,  3, -2,
0, 1, -2,  1, 0, 0,  0,  0,  0, -3,  6, -3,  0,  2, -4,  2 ,
0, 0,  0,  0, 0, 0,  0,  0,  0,  3, -6,  3,  0, -2,  4, -2,
0, 0,  0,  0, 0, 0,  0,  0,  0,  0, -3,  3,  0,  0,  2, -2,
0, 0, -1,  1, 0, 0,  0,  0,  0,  0,  3, -3,  0,  0, -2,  2 ,
0, 0,  0,  0, 0, 1, -2,  1,  0, -2,  4, -2,  0,  1, -2,  1,
0, 0,  0,  0, 0, 0,  0,  0,  0, -1,  2, -1,  0,  1, -2,  1,
0, 0,  0,  0, 0, 0,  0,  0,  0,  0,  1, -1,  0,  0, -1,  1,
0, 0,  0,  0, 0, 0, -1,  1,  0,  0,  2, -2,  0,  0, -1,  1
};



int glatnr(int *global_atom_index,int i)
{
    int atnr;

    if (global_atom_index == NULL) {
        atnr = i + 1;
    } else {
        atnr = global_atom_index[i] + 1;
    }

    return atnr;
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
		 const t_mdatoms *md,t_fcdata *fcd,
		 int *global_atom_index)
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
    dr   = dr2*gmx_invsqrt(dr2);                        /*  10          */
    temp = exp(-be*(dr-b0));                        /*  12          */
    
    if (temp == one)
      continue;

    omtemp   = one-temp;                               /*   1          */
    cbomtemp = cb*omtemp;                              /*   1          */
    vbond    = cbomtemp*omtemp;                        /*   1          */
    fbond    = -two*be*temp*cbomtemp*gmx_invsqrt(dr2);      /*   9          */
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
		 const t_mdatoms *md,t_fcdata *fcd,
		 int *global_atom_index)
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
      
    dr         = dr2*gmx_invsqrt(dr2);                      /*  10          */
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
		const t_mdatoms *md,t_fcdata *fcd,
		int *global_atom_index)
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
		dr2,bm2,
		glatnr(global_atom_index,ai),
		glatnr(global_atom_index,aj));
      
    omdr2obm2  = one - dr2/bm2;
    
    vbond      = -half*kb*bm2*log(omdr2obm2);
    fbond      = -kb/omdr2obm2;

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
	   const t_mdatoms *md,t_fcdata *fcd,
	   int *global_atom_index)
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
    dr   = dr2*gmx_invsqrt(dr2);		        /*  10		*/

    *dvdlambda += harmonic(forceparams[type].harmonic.krA,
			   forceparams[type].harmonic.krB,
			   forceparams[type].harmonic.rA,
			   forceparams[type].harmonic.rB,
			   dr,lambda,&vbond,&fbond);  /*  19  */

    if (dr2 == 0.0)
      continue;

    
    vtot  += vbond;/* 1*/
    fbond *= gmx_invsqrt(dr2);			/*   6		*/
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

real restraint_bonds(int nbonds,
                     const t_iatom forceatoms[],const t_iparams forceparams[],
                     const rvec x[],rvec f[],rvec fshift[],
                     const t_pbc *pbc,const t_graph *g,
                     real lambda,real *dvdlambda,
                     const t_mdatoms *md,t_fcdata *fcd,
                     int *global_atom_index)
{
    int  i,m,ki,ai,aj,type;
    real dr,dr2,fbond,vbond,fij,vtot;
    real L1;
    real low,dlow,up1,dup1,up2,dup2,k,dk;
    real drh,drh2;
    rvec dx;
    ivec dt;

    L1   = 1.0 - lambda;

    vtot = 0.0;
    for(i=0; (i<nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        aj   = forceatoms[i++];
        
        ki   = pbc_rvec_sub(pbc,x[ai],x[aj],dx);	/*   3 		*/
        dr2  = iprod(dx,dx);		             	/*   5		*/
        dr   = dr2*gmx_invsqrt(dr2);		        /*  10		*/

        low  = L1*forceparams[type].restraint.lowA + lambda*forceparams[type].restraint.lowB;
        dlow =   -forceparams[type].restraint.lowA +        forceparams[type].restraint.lowB;
        up1  = L1*forceparams[type].restraint.up1A + lambda*forceparams[type].restraint.up1B;
        dup1 =   -forceparams[type].restraint.up1A +        forceparams[type].restraint.up1B;
        up2  = L1*forceparams[type].restraint.up2A + lambda*forceparams[type].restraint.up2B;
        dup2 =   -forceparams[type].restraint.up2A +        forceparams[type].restraint.up2B;
        k    = L1*forceparams[type].restraint.kA   + lambda*forceparams[type].restraint.kB;
        dk   =   -forceparams[type].restraint.kA   +        forceparams[type].restraint.kB;
        /* 24 */

        if (dr < low)
        {
            drh   = dr - low;
            drh2  = drh*drh;
            vbond = 0.5*k*drh2;
            fbond = -k*drh;
            *dvdlambda += 0.5*dk*drh2 - k*dlow*drh;
        } /* 11 */
        else if (dr <= up1)
        {
            vbond = 0;
            fbond = 0;
        }
        else if (dr <= up2)
        {
            drh   = dr - up1;
            drh2  = drh*drh;
            vbond = 0.5*k*drh2;
            fbond = -k*drh;
            *dvdlambda += 0.5*dk*drh2 - k*dup1*drh;
        } /* 11	*/
        else
        {
            drh   = dr - up2;
            vbond = k*(up2 - up1)*(0.5*(up2 - up1) + drh);
            fbond = -k*(up2 - up1);
            *dvdlambda += dk*(up2 - up1)*(0.5*(up2 - up1) + drh)
                + k*(dup2 - dup1)*(up2 - up1 + drh)
                - k*(up2 - up1)*dup2;
        }
   
        if (dr2 == 0.0)
            continue;
        
        vtot  += vbond;/* 1*/
        fbond *= gmx_invsqrt(dr2);			/*   6		*/
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
	      const t_mdatoms *md,t_fcdata *fcd,
	      int *global_atom_index)
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
    if (debug)
      fprintf(debug,"POL: local ai = %d aj = %d ksh = %.3f\n",ai,aj,ksh);
  
    ki   = pbc_rvec_sub(pbc,x[ai],x[aj],dx);	/*   3 		*/
    dr2  = iprod(dx,dx);			/*   5		*/
    dr   = dr2*gmx_invsqrt(dr2);		        /*  10		*/

    *dvdlambda += harmonic(ksh,ksh,0,0,dr,lambda,&vbond,&fbond);  /*  19  */

    if (dr2 == 0.0)
      continue;
    
    vtot  += vbond;/* 1*/
    fbond *= gmx_invsqrt(dr2);			/*   6		*/

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
	       const t_mdatoms *md,t_fcdata *fcd,
	       int *global_atom_index)
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
      cprod(dOH1,dOH2,nW);
      
      /* Compute inverse length of normal vector 
       * (this one could be precomputed, but I'm too lazy now)
       */
      r_nW = gmx_invsqrt(iprod(nW,nW));
      /* This is for precision, but does not make a big difference,
       * it can go later.
       */
      r_OD = gmx_invsqrt(iprod(dOD,dOD)); 
      
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

static real do_1_thole(const rvec xi,const rvec xj,rvec fi,rvec fj,
		       const t_pbc *pbc,real qq,
		       rvec fshift[],real afac)
{
  rvec r12;
  real r12sq,r12_1,r12n,r12bar,v0,v1,fscal,ebar,fff;
  int  m,t;
    
  t      = pbc_rvec_sub(pbc,xi,xj,r12); /*  3 */
  
  r12sq  = iprod(r12,r12);              /*  5 */
  r12_1  = gmx_invsqrt(r12sq);              /*  5 */
  r12bar = afac/r12_1;                  /*  5 */
  v0     = qq*ONE_4PI_EPS0*r12_1;       /*  2 */
  ebar   = exp(-r12bar);                /*  5 */
  v1     = (1-(1+0.5*r12bar)*ebar);     /*  4 */
  fscal  = ((v0*r12_1)*v1 - v0*0.5*afac*ebar*(r12bar+1))*r12_1; /* 9 */
  if (debug)
    fprintf(debug,"THOLE: v0 = %.3f v1 = %.3f r12= % .3f r12bar = %.3f fscal = %.3f  ebar = %.3f\n",v0,v1,1/r12_1,r12bar,fscal,ebar);
  
  for(m=0; (m<DIM); m++) {
    fff    = fscal*r12[m];
    fi[m] += fff;
    fj[m] -= fff;             
    fshift[t][m]       += fff;
    fshift[CENTRAL][m] -= fff;
  } /* 15 */
  
  return v0*v1; /* 1 */
  /* 54 */
}

real thole_pol(int nbonds,
	       const t_iatom forceatoms[],const t_iparams forceparams[],
	       const rvec x[],rvec f[],rvec fshift[],
	       const t_pbc *pbc,const t_graph *g,
	       real lambda,real *dvdlambda,
	       const t_mdatoms *md,t_fcdata *fcd,
	       int *global_atom_index)
{
  /* Interaction between two pairs of particles with opposite charge */
  int i,type,a1,da1,a2,da2;
  real q1,q2,qq,a,al1,al2,afac;
  real V=0;
  
  for(i=0; (i<nbonds); ) {
    type  = forceatoms[i++];
    a1    = forceatoms[i++];
    da1   = forceatoms[i++];
    a2    = forceatoms[i++];
    da2   = forceatoms[i++];
    q1    = md->chargeA[da1];
    q2    = md->chargeA[da2];
    a     = forceparams[type].thole.a;
    al1   = forceparams[type].thole.alpha1;
    al2   = forceparams[type].thole.alpha2;
    qq    = q1*q2;
    afac  = a*pow(al1*al2,-1.0/6.0);
    V += do_1_thole(x[a1], x[a2], f[a1], f[a2], pbc, qq,fshift,afac);
    V += do_1_thole(x[da1],x[a2], f[da1],f[a2], pbc,-qq,fshift,afac);
    V += do_1_thole(x[a1], x[da2],f[a1], f[da2],pbc,-qq,fshift,afac);
    V += do_1_thole(x[da1],x[da2],f[da1],f[da2],pbc, qq,fshift,afac);
  }
  /* 290 flops */
  return V;
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
	    const t_mdatoms *md,t_fcdata *fcd,
	    int *global_atom_index)
{
  int  i,ai,aj,ak,t1,t2,type;
  rvec r_ij,r_kj;
  real cos_theta,cos_theta2,theta,dVdt,va,vtot;
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
    
    cos_theta2 = sqr(cos_theta);
    if (cos_theta2 < 1) {
      int  m;
      real st,sth;
      real cik,cii,ckk;
      real nrkj2,nrij2;
      rvec f_i,f_j,f_k;
      
      st  = dVdt*gmx_invsqrt(1 - cos_theta2);	/*  12		*/
      sth = st*cos_theta;			/*   1		*/
#ifdef DEBUG
      if (debug)
	fprintf(debug,"ANGLES: theta = %10g  vth = %10g  dV/dtheta = %10g\n",
		theta*RAD2DEG,va,dVdt);
#endif
      nrkj2=iprod(r_kj,r_kj);			/*   5		*/
      nrij2=iprod(r_ij,r_ij);
      
      cik=st*gmx_invsqrt(nrkj2*nrij2);		/*  12		*/ 
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
    }                                           /* 161 TOTAL	*/
  }
  return vtot;
}

real urey_bradley(int nbonds,
		  const t_iatom forceatoms[],const t_iparams forceparams[],
		  const rvec x[],rvec f[],rvec fshift[],
		  const t_pbc *pbc,const t_graph *g,
		  real lambda,real *dvdlambda,
		  const t_mdatoms *md,t_fcdata *fcd,
		  int *global_atom_index)
{
  int  i,m,ai,aj,ak,t1,t2,type,ki;
  rvec r_ij,r_kj,r_ik;
  real cos_theta,cos_theta2,theta;
  real dVdt,va,vtot,kth,th0,kUB,r13,dr,dr2,vbond,fbond,fik;
  ivec jt,dt_ij,dt_kj,dt_ik;
  
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
    dr   = dr2*gmx_invsqrt(dr2);		        /*  10		*/

    *dvdlambda += harmonic(kUB,kUB,r13,r13,dr,lambda,&vbond,&fbond); /*  19  */

    cos_theta2 = sqr(cos_theta);                /*   1		*/
    if (cos_theta2 < 1) {
      real st,sth;
      real cik,cii,ckk;
      real nrkj2,nrij2;
      rvec f_i,f_j,f_k;
      
      st  = dVdt*gmx_invsqrt(1 - cos_theta2);	/*  12		*/
      sth = st*cos_theta;			/*   1		*/
#ifdef DEBUG
      if (debug)
	fprintf(debug,"ANGLES: theta = %10g  vth = %10g  dV/dtheta = %10g\n",
		theta*RAD2DEG,va,dVdt);
#endif
      nrkj2=iprod(r_kj,r_kj);			/*   5		*/
      nrij2=iprod(r_ij,r_ij);
      
      cik=st*gmx_invsqrt(nrkj2*nrij2);		/*  12		*/ 
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
    }                                           /* 161 TOTAL	*/
    /* Time for the bond calculations */
    if (dr2 == 0.0)
      continue;

    vtot  += vbond;  /* 1*/
    fbond *= gmx_invsqrt(dr2);			/*   6		*/
    
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
		    const t_mdatoms *md,t_fcdata *fcd,
		    int *global_atom_index)
{
  int  i,j,ai,aj,ak,t1,t2,type;
  rvec r_ij,r_kj;
  real cos_theta,cos_theta2,theta,dt,dVdt,va,dtp,c,vtot;
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
    
    cos_theta2 = sqr(cos_theta);                /*   1		*/
    if (cos_theta2 < 1) {
      int  m;
      real st,sth;
      real cik,cii,ckk;
      real nrkj2,nrij2;
      rvec f_i,f_j,f_k;
      
      st  = dVdt*gmx_invsqrt(1 - cos_theta2);    	/*  12		*/
      sth = st*cos_theta;			/*   1		*/
#ifdef DEBUG
      if (debug)
	fprintf(debug,"ANGLES: theta = %10g  vth = %10g  dV/dtheta = %10g\n",
		theta*RAD2DEG,va,dVdt);
#endif
      nrkj2=iprod(r_kj,r_kj);			/*   5		*/
      nrij2=iprod(r_ij,r_ij);
      
      cik=st*gmx_invsqrt(nrkj2*nrij2);		/*  12		*/ 
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
    }                                           /* 153 TOTAL	*/
  }
  return vtot;
}

real dih_angle(const rvec xi,const rvec xj,const rvec xk,const rvec xl,
               const t_pbc *pbc,
               rvec r_ij,rvec r_kj,rvec r_kl,rvec m,rvec n,
               real *sign,int *t1,int *t2,int *t3)
{
  real ipr,phi;

  *t1 = pbc_rvec_sub(pbc,xi,xj,r_ij);       		/*  3 		*/
  *t2 = pbc_rvec_sub(pbc,xk,xj,r_kj);			/*  3		*/
  *t3 = pbc_rvec_sub(pbc,xk,xl,r_kl);			/*  3		*/

  cprod(r_ij,r_kj,m); 			/*  9 		*/
  cprod(r_kj,r_kl,n);			/*  9		*/
  phi=gmx_angle(m,n); 			/* 49 (assuming 25 for atan2) */
  ipr=iprod(r_ij,n); 			/*  5 		*/
  (*sign)=(ipr<0.0)?-1.0:1.0;
  phi=(*sign)*phi; 			/*  1		*/
					/* 82 TOTAL	*/
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
  real iprm,iprn,nrkj,nrkj2;
  real a,p,q,toler;
  ivec jt,dt_ij,dt_kj,dt_lj;  
  
  iprm  = iprod(m,m);		/*  5 	*/
  iprn  = iprod(n,n);		/*  5	*/
  nrkj2 = iprod(r_kj,r_kj);	/*  5	*/
  toler = nrkj2*GMX_REAL_EPS;
  if ((iprm > toler) && (iprn > toler)) {
    nrkj  = nrkj2*gmx_invsqrt(nrkj2);	/* 10	*/
    a     = -ddphi*nrkj/iprm;	/* 11	*/
    svmul(a,m,f_i);		/*  3	*/
    a     = ddphi*nrkj/iprn;	/* 11	*/
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
    } else if (pbc) {
      t3 = pbc_rvec_sub(pbc,x[l],x[j],dx_jl);
    } else {
      t3 = CENTRAL;
    }
    
    rvec_inc(fshift[t1],f_i);
    rvec_dec(fshift[CENTRAL],f_j);
    rvec_dec(fshift[t2],f_k);
    rvec_inc(fshift[t3],f_l);
  }
  /* 112 TOTAL 	*/
}


real dopdihs(real cpA,real cpB,real phiA,real phiB,int mult,
	     real phi,real lambda,real *V,real *F)
{
  real v,dvdl,mdphi,v1,sdphi,ddphi;
  real L1   = 1.0 - lambda;
  real ph0  = (L1*phiA + lambda*phiB)*DEG2RAD;
  real dph0 = (phiB - phiA)*DEG2RAD;
  real cp   = L1*cpA + lambda*cpB;
  
  mdphi =  mult*phi - ph0;
  sdphi = sin(mdphi);
  ddphi = -cp*mult*sdphi;
  v1    = 1.0 + cos(mdphi);
  v     = cp*v1;
  
  dvdl  = (cpB - cpA)*v1 + cp*dph0*sdphi;
  
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
  real L1   = 1.0 - lambda;
  real ph0  = (L1*phiA + lambda*phiB)*DEG2RAD;
  real dph0 = (phiB - phiA)*DEG2RAD;
  real cp   = L1*cpA + lambda*cpB;
  
  mdphi = mult*(phi-ph0);
  sdphi = sin(mdphi);
  ddphi = cp*mult*sdphi;
  v1    = 1.0-cos(mdphi);
  v     = cp*v1;
  
  dvdl  = (cpB-cpA)*v1 + cp*dph0*sdphi;
  
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
	   const t_mdatoms *md,t_fcdata *fcd,
	   int *global_atom_index)
{
  int  i,type,ai,aj,ak,al;
  int  t1,t2,t3;
  rvec r_ij,r_kj,r_kl,m,n;
  real phi,sign,ddphi,vpd,vtot;

  vtot = 0.0;

  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    al   = forceatoms[i++];
    
    phi=dih_angle(x[ai],x[aj],x[ak],x[al],pbc,r_ij,r_kj,r_kl,m,n,
                  &sign,&t1,&t2,&t3);			/*  84 		*/

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
    fprintf(debug,"pdih: (%d,%d,%d,%d) phi=%g\n",
	    ai,aj,ak,al,phi);
#endif
  } /* 223 TOTAL 	*/

  return vtot;
}



real idihs(int nbonds,
	   const t_iatom forceatoms[],const t_iparams forceparams[],
	   const rvec x[],rvec f[],rvec fshift[],
	   const t_pbc *pbc,const t_graph *g,
	   real lambda,real *dvdlambda,
	   const t_mdatoms *md,t_fcdata *fcd,
	   int *global_atom_index)
{
  int  i,type,ai,aj,ak,al;
  int  t1,t2,t3;
  real phi,phi0,dphi0,ddphi,sign,vtot;
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
                  &sign,&t1,&t2,&t3);			/*  84		*/
    
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

    kk    = L1*kA + lambda*kB;
    phi0  = (L1*pA + lambda*pB)*DEG2RAD;
    dphi0 = (pB - pA)*DEG2RAD;

    /* dp = (phi-phi0), modulo (-pi,pi) */
    dp = phi-phi0;  
    /* dp cannot be outside (-2*pi,2*pi) */
    if (dp >= M_PI)
      dp -= 2*M_PI;
    else if(dp < -M_PI)
      dp += 2*M_PI;
    
    dp2 = dp*dp;

    vtot += 0.5*kk*dp2;
    ddphi = -kk*dp;
    
    dvdl += 0.5*(kB - kA)*dp2 - kk*dphi0*dp;

    do_dih_fup(ai,aj,ak,al,(real)(-ddphi),r_ij,r_kj,r_kl,m,n,
	       f,fshift,pbc,g,x,t1,t2,t3);			/* 112		*/
    /* 217 TOTAL	*/
#ifdef DEBUG
    if (debug)
      fprintf(debug,"idih: (%d,%d,%d,%d) phi=%g\n",
	      ai,aj,ak,al,phi);
#endif
  }
  
  *dvdlambda += dvdl;
  return vtot;
}


real posres(int nbonds,
            const t_iatom forceatoms[],const t_iparams forceparams[],
            const rvec x[],rvec f[],rvec vir_diag,
            t_pbc *pbc,
            real lambda,real *dvdlambda,
            int refcoord_scaling,int ePBC,rvec comA,rvec comB)
{
    int  i,ai,m,d,type,ki,npbcdim=0;
    const t_iparams *pr;
    real L1;
    real vtot,kk,fm;
    real posA,posB,ref=0;
    rvec comA_sc,comB_sc,rdist,dpdl,pos,dx;

    npbcdim = ePBC2npbcdim(ePBC);

    if (refcoord_scaling == erscCOM)
    {
        clear_rvec(comA_sc);
        clear_rvec(comB_sc);
        for(m=0; m<npbcdim; m++)
        {
            for(d=m; d<npbcdim; d++)
            {
                comA_sc[m] += comA[d]*pbc->box[d][m];
                comB_sc[m] += comB[d]*pbc->box[d][m];
            }
        }
    }

    L1 = 1.0 - lambda;

    vtot = 0.0;
    for(i=0; (i<nbonds); )
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        pr   = &forceparams[type];
        
        for(m=0; m<DIM; m++)
        {
            posA = forceparams[type].posres.pos0A[m];
            posB = forceparams[type].posres.pos0B[m];
            if (m < npbcdim)
            {
                switch (refcoord_scaling)
                {
                case erscNO:
                    ref      = 0;
                    rdist[m] = L1*posA + lambda*posB;
                    dpdl[m]  = posB - posA;
                    break;
                case erscALL:
                    /* Box relative coordinates are stored for dimensions with pbc */
                    posA *= pbc->box[m][m];
                    posB *= pbc->box[m][m];
                    for(d=m+1; d<npbcdim; d++)
                    {
                        posA += forceparams[type].posres.pos0A[d]*pbc->box[d][m];
                        posB += forceparams[type].posres.pos0B[d]*pbc->box[d][m];
                    }
                    ref      = L1*posA + lambda*posB;
                    rdist[m] = 0;
                    dpdl[m]  = posB - posA;
                    break;
                case erscCOM:
                    ref      = L1*comA_sc[m] + lambda*comB_sc[m];
                    rdist[m] = L1*posA       + lambda*posB;
                    dpdl[m]  = comB_sc[m] - comA_sc[m] + posB - posA;
                    break;
                }
            }
            else
            {
                ref      = L1*posA + lambda*posB;
                rdist[m] = 0;
                dpdl[m]  = posB - posA;
            }

            /* We do pbc_dx with ref+rdist,
             * since with only ref we can be up to half a box vector wrong.
             */
            pos[m] = ref + rdist[m];
        }

        if (pbc)
        {
            pbc_dx(pbc,x[ai],pos,dx);
        }
        else
        {
            rvec_sub(x[ai],pos,dx);
        }

        for (m=0; (m<DIM); m++)
        {
            kk          = L1*pr->posres.fcA[m] + lambda*pr->posres.fcB[m];
            fm          = -kk*dx[m];
            f[ai][m]   += fm;
            vtot       += 0.5*kk*dx[m]*dx[m];
            *dvdlambda +=
                0.5*(pr->posres.fcB[m] - pr->posres.fcA[m])*dx[m]*dx[m]
                -fm*dpdl[m];

            /* Here we correct for the pbc_dx which included rdist */
            vir_diag[m] -= 0.5*(dx[m] + rdist[m])*fm;
        }
    }

    return vtot;
}

static real low_angres(int nbonds,
		       const t_iatom forceatoms[],const t_iparams forceparams[],
		       const rvec x[],rvec f[],rvec fshift[],
		       const t_pbc *pbc,const t_graph *g,
		       real lambda,real *dvdlambda,
		       gmx_bool bZAxis)
{
  int  i,m,type,ai,aj,ak,al;
  int  t1,t2;
  real phi,cos_phi,cos_phi2,vid,vtot,dVdphi;
  rvec r_ij,r_kl,f_i,f_k={0,0,0};
  real st,sth,nrij2,nrkl2,c,cij,ckl;

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

    cos_phi2 = sqr(cos_phi);                    /*   1		*/
    if (cos_phi2 < 1) {
      st  = -dVdphi*gmx_invsqrt(1 - cos_phi2);      /*  12		*/
      sth = st*cos_phi;				/*   1		*/
      nrij2 = iprod(r_ij,r_ij);			/*   5		*/
      nrkl2 = iprod(r_kl,r_kl);                 /*   5          */
      
      c   = st*gmx_invsqrt(nrij2*nrkl2);		/*  11		*/ 
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
  }

  return vtot;  /*  184 / 157 (bZAxis)  total  */
}

real angres(int nbonds,
	    const t_iatom forceatoms[],const t_iparams forceparams[],
	    const rvec x[],rvec f[],rvec fshift[],
	    const t_pbc *pbc,const t_graph *g,
	    real lambda,real *dvdlambda,
	    const t_mdatoms *md,t_fcdata *fcd,
	    int *global_atom_index)
{
  return low_angres(nbonds,forceatoms,forceparams,x,f,fshift,pbc,g,
		    lambda,dvdlambda,FALSE);
}

real angresz(int nbonds,
	     const t_iatom forceatoms[],const t_iparams forceparams[],
	     const rvec x[],rvec f[],rvec fshift[],
	     const t_pbc *pbc,const t_graph *g,
	     real lambda,real *dvdlambda,
	     const t_mdatoms *md,t_fcdata *fcd,
	     int *global_atom_index)
{
  return low_angres(nbonds,forceatoms,forceparams,x,f,fshift,pbc,g,
		    lambda,dvdlambda,TRUE);
}


real unimplemented(int nbonds,
		   const t_iatom forceatoms[],const t_iparams forceparams[],
		   const rvec x[],rvec f[],rvec fshift[],
		   const t_pbc *pbc,const t_graph *g,
		   real lambda,real *dvdlambda,
		   const t_mdatoms *md,t_fcdata *fcd,
		   int *global_atom_index)
{
  gmx_impl("*** you are using a not implemented function");

  return 0.0; /* To make the compiler happy */
}

real rbdihs(int nbonds,
	    const t_iatom forceatoms[],const t_iparams forceparams[],
	    const rvec x[],rvec f[],rvec fshift[],
	    const t_pbc *pbc,const t_graph *g,
	    real lambda,real *dvdlambda,
	    const t_mdatoms *md,t_fcdata *fcd,
	    int *global_atom_index)
{
  const real c0=0.0,c1=1.0,c2=2.0,c3=3.0,c4=4.0,c5=5.0;
  int  type,ai,aj,ak,al,i,j;
  int  t1,t2,t3;
  rvec r_ij,r_kj,r_kl,m,n;
  real parmA[NR_RBDIHS];
  real parmB[NR_RBDIHS];
  real parm[NR_RBDIHS];
  real cos_phi,phi,rbp,rbpBA;
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
                    &sign,&t1,&t2,&t3);			/*  84		*/

    /* Change to polymer convention */
    if (phi < c0)
      phi += M_PI;
    else
      phi -= M_PI;			/*   1		*/
      
    cos_phi = cos(phi);		
    /* Beware of accuracy loss, cannot use 1-sqrt(cos^2) ! */
    sin_phi = sin(phi);

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

int cmap_setup_grid_index(int ip, int grid_spacing, int *ipm1, int *ipp1, int *ipp2)
{
	int im1, ip1, ip2;
	
	if(ip<0)
	{
		ip = ip + grid_spacing - 1;
	}
	else if(ip > grid_spacing)
	{
		ip = ip - grid_spacing - 1;
	}
	
	im1 = ip - 1;
	ip1 = ip + 1;
	ip2 = ip + 2;
	
	if(ip == 0)
	{
		im1 = grid_spacing - 1;
	}
	else if(ip == grid_spacing-2)
	{
		ip2 = 0;
	}
	else if(ip == grid_spacing-1)
	{
		ip1 = 0;
		ip2 = 1;
	}
	
	*ipm1 = im1;
	*ipp1 = ip1;
	*ipp2 = ip2;
	
	return ip;
	
}

real cmap_dihs(int nbonds,
			   const t_iatom forceatoms[],const t_iparams forceparams[],
               const gmx_cmap_t *cmap_grid,
			   const rvec x[],rvec f[],rvec fshift[],
			   const t_pbc *pbc,const t_graph *g,
			   real lambda,real *dvdlambda,
			   const t_mdatoms *md,t_fcdata *fcd,
			   int *global_atom_index)
{
	int i,j,k,n,idx;
	int ai,aj,ak,al,am;
	int a1i,a1j,a1k,a1l,a2i,a2j,a2k,a2l;
	int type,cmapA;
	int t11,t21,t31,t12,t22,t32;
	int iphi1,ip1m1,ip1p1,ip1p2;
	int iphi2,ip2m1,ip2p1,ip2p2;
	int l1,l2,l3,l4;
	int pos1,pos2,pos3,pos4,tmp;
	
	real ty[4],ty1[4],ty2[4],ty12[4],tc[16],tx[16];
	real phi1,psi1,cos_phi1,sin_phi1,sign1,xphi1;
	real phi2,psi2,cos_phi2,sin_phi2,sign2,xphi2;
	real dx,xx,tt,tu,e,df1,df2,ddf1,ddf2,ddf12,vtot;
	real ra21,rb21,rg21,rg1,rgr1,ra2r1,rb2r1,rabr1;
	real ra22,rb22,rg22,rg2,rgr2,ra2r2,rb2r2,rabr2;
	real fg1,hg1,fga1,hgb1,gaa1,gbb1;
	real fg2,hg2,fga2,hgb2,gaa2,gbb2;
	real fac;
	
	rvec r1_ij, r1_kj, r1_kl,m1,n1;
	rvec r2_ij, r2_kj, r2_kl,m2,n2;
	rvec f1_i,f1_j,f1_k,f1_l;
	rvec f2_i,f2_j,f2_k,f2_l;
	rvec a1,b1,a2,b2;
	rvec f1,g1,h1,f2,g2,h2;
	rvec dtf1,dtg1,dth1,dtf2,dtg2,dth2;
	ivec jt1,dt1_ij,dt1_kj,dt1_lj;
	ivec jt2,dt2_ij,dt2_kj,dt2_lj;

    const real *cmapd;
	
	int loop_index[4][4] = {
		{0,4,8,12},
		{1,5,9,13},
		{2,6,10,14},
		{3,7,11,15}
	};
	
	/* Total CMAP energy */
	vtot = 0;
	
	for(n=0;n<nbonds; )
	{
		/* Five atoms are involved in the two torsions */
		type   = forceatoms[n++];
		ai     = forceatoms[n++];
		aj     = forceatoms[n++];
		ak     = forceatoms[n++];
		al     = forceatoms[n++];
		am     = forceatoms[n++];
		
		/* Which CMAP type is this */
		cmapA = forceparams[type].cmap.cmapA;
        cmapd = cmap_grid->cmapdata[cmapA].cmap;
		
		/* First torsion */
		a1i   = ai;
		a1j   = aj;
		a1k   = ak;
		a1l   = al;
		
		phi1  = dih_angle(x[a1i], x[a1j], x[a1k], x[a1l], pbc, r1_ij, r1_kj, r1_kl, m1, n1,
						   &sign1, &t11, &t21, &t31); /* 84 */
		
        cos_phi1 = cos(phi1);
        
		a1[0] = r1_ij[1]*r1_kj[2]-r1_ij[2]*r1_kj[1];
		a1[1] = r1_ij[2]*r1_kj[0]-r1_ij[0]*r1_kj[2];
		a1[2] = r1_ij[0]*r1_kj[1]-r1_ij[1]*r1_kj[0]; /* 9 */
		
		b1[0] = r1_kl[1]*r1_kj[2]-r1_kl[2]*r1_kj[1];
		b1[1] = r1_kl[2]*r1_kj[0]-r1_kl[0]*r1_kj[2];
		b1[2] = r1_kl[0]*r1_kj[1]-r1_kl[1]*r1_kj[0]; /* 9 */
		
		tmp = pbc_rvec_sub(pbc,x[a1l],x[a1k],h1);
		
		ra21  = iprod(a1,a1);       /* 5 */
		rb21  = iprod(b1,b1);       /* 5 */
		rg21  = iprod(r1_kj,r1_kj); /* 5 */
		rg1   = sqrt(rg21);
		
		rgr1  = 1.0/rg1;
		ra2r1 = 1.0/ra21;
		rb2r1 = 1.0/rb21;
		rabr1 = sqrt(ra2r1*rb2r1);
		
		sin_phi1 = rg1 * rabr1 * iprod(a1,h1) * (-1);
		
		if(cos_phi1 < -0.5 || cos_phi1 > 0.5)
		{
			phi1 = asin(sin_phi1);
			
			if(cos_phi1 < 0)
			{
				if(phi1 > 0)
				{
					phi1 = M_PI - phi1;
				}
				else
				{
					phi1 = -M_PI - phi1;
				}
			}
		}
		else
		{
			phi1 = acos(cos_phi1);
			
			if(sin_phi1 < 0)
			{
				phi1 = -phi1;
			}
		}
		
		xphi1 = phi1 + M_PI; /* 1 */
		
		/* Second torsion */
		a2i   = aj;
		a2j   = ak;
		a2k   = al;
		a2l   = am;
		
		phi2  = dih_angle(x[a2i], x[a2j], x[a2k], x[a2l], pbc, r2_ij, r2_kj, r2_kl, m2, n2,
						  &sign2, &t12, &t22, &t32); /* 84 */
		
        cos_phi2 = cos(phi2);

		a2[0] = r2_ij[1]*r2_kj[2]-r2_ij[2]*r2_kj[1];
		a2[1] = r2_ij[2]*r2_kj[0]-r2_ij[0]*r2_kj[2];
		a2[2] = r2_ij[0]*r2_kj[1]-r2_ij[1]*r2_kj[0]; /* 9 */
		
		b2[0] = r2_kl[1]*r2_kj[2]-r2_kl[2]*r2_kj[1];
		b2[1] = r2_kl[2]*r2_kj[0]-r2_kl[0]*r2_kj[2];
		b2[2] = r2_kl[0]*r2_kj[1]-r2_kl[1]*r2_kj[0]; /* 9 */
		
		tmp = pbc_rvec_sub(pbc,x[a2l],x[a2k],h2);
		
		ra22  = iprod(a2,a2);         /* 5 */
		rb22  = iprod(b2,b2);         /* 5 */
		rg22  = iprod(r2_kj,r2_kj);   /* 5 */
		rg2   = sqrt(rg22);
		
		rgr2  = 1.0/rg2;
		ra2r2 = 1.0/ra22;
		rb2r2 = 1.0/rb22;
		rabr2 = sqrt(ra2r2*rb2r2);
		
		sin_phi2 = rg2 * rabr2 * iprod(a2,h2) * (-1);
		
		if(cos_phi2 < -0.5 || cos_phi2 > 0.5)
		{
			phi2 = asin(sin_phi2);
			
			if(cos_phi2 < 0)
			{
				if(phi2 > 0)
				{
					phi2 = M_PI - phi2;
				}
				else
				{
					phi2 = -M_PI - phi2;
				}
			}
		}
		else
		{
			phi2 = acos(cos_phi2);
			
			if(sin_phi2 < 0)
			{
				phi2 = -phi2;
			}
		}
		
		xphi2 = phi2 + M_PI; /* 1 */
		
		/* Range mangling */
		if(xphi1<0)
		{
			xphi1 = xphi1 + 2*M_PI;
		}
		else if(xphi1>=2*M_PI)
		{
			xphi1 = xphi1 - 2*M_PI;
		}
		
		if(xphi2<0)
		{
			xphi2 = xphi2 + 2*M_PI;
		}
		else if(xphi2>=2*M_PI)
		{
			xphi2 = xphi2 - 2*M_PI;
		}
		
		/* Number of grid points */
		dx = 2*M_PI / cmap_grid->grid_spacing;
		
		/* Where on the grid are we */
		iphi1 = (int)(xphi1/dx);
		iphi2 = (int)(xphi2/dx);
		
		iphi1 = cmap_setup_grid_index(iphi1, cmap_grid->grid_spacing, &ip1m1,&ip1p1,&ip1p2);
		iphi2 = cmap_setup_grid_index(iphi2, cmap_grid->grid_spacing, &ip2m1,&ip2p1,&ip2p2);
		
		pos1    = iphi1*cmap_grid->grid_spacing+iphi2;
		pos2    = ip1p1*cmap_grid->grid_spacing+iphi2;
		pos3    = ip1p1*cmap_grid->grid_spacing+ip2p1;
		pos4    = iphi1*cmap_grid->grid_spacing+ip2p1;
		
		ty[0]   = cmapd[pos1*4];
		ty[1]   = cmapd[pos2*4];
		ty[2]   = cmapd[pos3*4];
		ty[3]   = cmapd[pos4*4];
		
		ty1[0]   = cmapd[pos1*4+1];
		ty1[1]   = cmapd[pos2*4+1];
		ty1[2]   = cmapd[pos3*4+1];
		ty1[3]   = cmapd[pos4*4+1];
		
		ty2[0]   = cmapd[pos1*4+2];
		ty2[1]   = cmapd[pos2*4+2];
		ty2[2]   = cmapd[pos3*4+2];
		ty2[3]   = cmapd[pos4*4+2];
		
		ty12[0]   = cmapd[pos1*4+3];
		ty12[1]   = cmapd[pos2*4+3];
		ty12[2]   = cmapd[pos3*4+3];
		ty12[3]   = cmapd[pos4*4+3];
		
		/* Switch to degrees */
		dx = 360.0 / cmap_grid->grid_spacing;
		xphi1 = xphi1 * RAD2DEG;
		xphi2 = xphi2 * RAD2DEG; 
		
		for(i=0;i<4;i++) /* 16 */
		{
			tx[i] = ty[i];
			tx[i+4] = ty1[i]*dx;
			tx[i+8] = ty2[i]*dx;
			tx[i+12] = ty12[i]*dx*dx;
		}
		
		idx=0;
		for(i=0;i<4;i++) /* 1056 */
		{
			for(j=0;j<4;j++)
			{
				xx = 0;
				for(k=0;k<16;k++)
				{
					xx = xx + cmap_coeff_matrix[k*16+idx]*tx[k];
				}
				
				idx++;
				tc[i*4+j]=xx;
			}
		}
		
		tt    = (xphi1-iphi1*dx)/dx;
		tu    = (xphi2-iphi2*dx)/dx;
		
		e     = 0;
		df1   = 0;
		df2   = 0;
		ddf1  = 0;
		ddf2  = 0;
		ddf12 = 0;
		
		for(i=3;i>=0;i--)
		{
			l1 = loop_index[i][3];
			l2 = loop_index[i][2];
			l3 = loop_index[i][1];
			
			e     = tt * e    + ((tc[i*4+3]*tu+tc[i*4+2])*tu + tc[i*4+1])*tu+tc[i*4];
			df1   = tu * df1  + (3.0*tc[l1]*tt+2.0*tc[l2])*tt+tc[l3];
			df2   = tt * df2  + (3.0*tc[i*4+3]*tu+2.0*tc[i*4+2])*tu+tc[i*4+1];
			ddf1  = tu * ddf1 + 2.0*3.0*tc[l1]*tt+2.0*tc[l2];
			ddf2  = tt * ddf2 + 2.0*3.0*tc[4*i+3]*tu+2.0*tc[4*i+2];
		}
		
		ddf12 = tc[5] + 2.0*tc[9]*tt + 3.0*tc[13]*tt*tt + 2.0*tu*(tc[6]+2.0*tc[10]*tt+3.0*tc[14]*tt*tt) +
		3.0*tu*tu*(tc[7]+2.0*tc[11]*tt+3.0*tc[15]*tt*tt);
		
		fac     = RAD2DEG/dx;
		df1     = df1   * fac;
		df2     = df2   * fac;
		ddf1    = ddf1  * fac * fac;
		ddf2    = ddf2  * fac * fac;
		ddf12   = ddf12 * fac * fac;
		
		/* CMAP energy */
		vtot += e;
		
		/* Do forces - first torsion */
		fg1       = iprod(r1_ij,r1_kj);
		hg1       = iprod(r1_kl,r1_kj);
		fga1      = fg1*ra2r1*rgr1;
		hgb1      = hg1*rb2r1*rgr1;
		gaa1      = -ra2r1*rg1;
		gbb1      = rb2r1*rg1;
		
		for(i=0;i<DIM;i++)
		{
			dtf1[i]   = gaa1 * a1[i];
			dtg1[i]   = fga1 * a1[i] - hgb1 * b1[i];
			dth1[i]   = gbb1 * b1[i];
			
			f1[i]     = df1  * dtf1[i];
			g1[i]     = df1  * dtg1[i];
			h1[i]     = df1  * dth1[i];
			
			f1_i[i]   =  f1[i];
			f1_j[i]   = -f1[i] - g1[i];
			f1_k[i]   =  h1[i] + g1[i];
			f1_l[i]   = -h1[i];
			
			f[a1i][i] = f[a1i][i] + f1_i[i];
			f[a1j][i] = f[a1j][i] + f1_j[i]; /* - f1[i] - g1[i] */                                                            
			f[a1k][i] = f[a1k][i] + f1_k[i]; /* h1[i] + g1[i] */                                                            
			f[a1l][i] = f[a1l][i] + f1_l[i]; /* h1[i] */                                                                       
		}
		
		/* Do forces - second torsion */
		fg2       = iprod(r2_ij,r2_kj);
		hg2       = iprod(r2_kl,r2_kj);
		fga2      = fg2*ra2r2*rgr2;
		hgb2      = hg2*rb2r2*rgr2;
		gaa2      = -ra2r2*rg2;
		gbb2      = rb2r2*rg2;
		
		for(i=0;i<DIM;i++)
		{
			dtf2[i]   = gaa2 * a2[i];
			dtg2[i]   = fga2 * a2[i] - hgb2 * b2[i];
			dth2[i]   = gbb2 * b2[i];
			
			f2[i]     = df2  * dtf2[i];
			g2[i]     = df2  * dtg2[i];
			h2[i]     = df2  * dth2[i];
			
			f2_i[i]   =  f2[i];
			f2_j[i]   = -f2[i] - g2[i];
			f2_k[i]   =  h2[i] + g2[i];
			f2_l[i]   = -h2[i];
			
			f[a2i][i] = f[a2i][i] + f2_i[i]; /* f2[i] */                                                                        
			f[a2j][i] = f[a2j][i] + f2_j[i]; /* - f2[i] - g2[i] */                                                              
			f[a2k][i] = f[a2k][i] + f2_k[i]; /* h2[i] + g2[i] */                   		
			f[a2l][i] = f[a2l][i] + f2_l[i]; /* - h2[i] */                                                                      
		}
		
		/* Shift forces */
		if(g)
		{
			copy_ivec(SHIFT_IVEC(g,a1j), jt1);
			ivec_sub(SHIFT_IVEC(g,a1i),  jt1,dt1_ij);
			ivec_sub(SHIFT_IVEC(g,a1k),  jt1,dt1_kj);
			ivec_sub(SHIFT_IVEC(g,a1l),  jt1,dt1_lj);
			t11 = IVEC2IS(dt1_ij);
			t21 = IVEC2IS(dt1_kj);
			t31 = IVEC2IS(dt1_lj);
			
			copy_ivec(SHIFT_IVEC(g,a2j), jt2);
			ivec_sub(SHIFT_IVEC(g,a2i),  jt2,dt2_ij);
			ivec_sub(SHIFT_IVEC(g,a2k),  jt2,dt2_kj);
			ivec_sub(SHIFT_IVEC(g,a2l),  jt2,dt2_lj);
			t12 = IVEC2IS(dt2_ij);
			t22 = IVEC2IS(dt2_kj);
			t32 = IVEC2IS(dt2_lj);
		}
		else if(pbc)
		{
			t31 = pbc_rvec_sub(pbc,x[a1l],x[a1j],h1);
			t32 = pbc_rvec_sub(pbc,x[a2l],x[a2j],h2);
		}
		else
		{
			t31 = CENTRAL;
			t32 = CENTRAL;
		}
		
		rvec_inc(fshift[t11],f1_i);
		rvec_inc(fshift[CENTRAL],f1_j);
		rvec_inc(fshift[t21],f1_k);
		rvec_inc(fshift[t31],f1_l);
		
		rvec_inc(fshift[t21],f2_i);
		rvec_inc(fshift[CENTRAL],f2_j);
		rvec_inc(fshift[t22],f2_k);
		rvec_inc(fshift[t32],f2_l);
	}	
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
	      const t_mdatoms *md,t_fcdata *fcd,
	      int *global_atom_index)
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

  costh=cos_angle(r_ij,r_kj); 		/* 25		*/
					/* 41 TOTAL	*/
  return costh;
}

real g96angles(int nbonds,
	       const t_iatom forceatoms[],const t_iparams forceparams[],
	       const rvec x[],rvec f[],rvec fshift[],
	       const t_pbc *pbc,const t_graph *g,
	       real lambda,real *dvdlambda,
	       const t_mdatoms *md,t_fcdata *fcd,
	       int *global_atom_index)
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
    
    rij_1    = gmx_invsqrt(iprod(r_ij,r_ij));
    rkj_1    = gmx_invsqrt(iprod(r_kj,r_kj));
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
		     const t_mdatoms *md,t_fcdata *fcd,
		     int *global_atom_index)
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
		      const t_mdatoms *md,t_fcdata *fcd,
		      int *global_atom_index)
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

static real bonded_tab(const char *type,int table_nr,
		       const bondedtable_t *table,real kA,real kB,real r,
		       real lambda,real *V,real *F)
{
  real k,tabscale,*VFtab,rt,eps,eps2,Yt,Ft,Geps,Heps2,Fp,VV,FF;
  int  n0,nnn;
  real v,f,dvdl;

  k = (1.0 - lambda)*kA + lambda*kB;

  tabscale = table->scale;
  VFtab    = table->tab;
  
  rt    = r*tabscale;
  n0    = rt;
  if (n0 >= table->n) {
    gmx_fatal(FARGS,"A tabulated %s interaction table number %d is out of the table range: r %f, between table indices %d and %d, table length %d",
	      type,table_nr,r,n0,n0+1,table->n);
  }
  eps   = rt - n0;
  eps2  = eps*eps;
  nnn   = 4*n0;
  Yt    = VFtab[nnn];
  Ft    = VFtab[nnn+1];
  Geps  = VFtab[nnn+2]*eps;
  Heps2 = VFtab[nnn+3]*eps2;
  Fp    = Ft + Geps + Heps2;
  VV    = Yt + Fp*eps;
  FF    = Fp + Geps + 2.0*Heps2;
  
  *F    = -k*FF*tabscale;
  *V    = k*VV;
  dvdl  = (kB - kA)*VV;
  
  return dvdl;
  
  /* That was 22 flops */
}

real tab_bonds(int nbonds,
	       const t_iatom forceatoms[],const t_iparams forceparams[],
	       const rvec x[],rvec f[],rvec fshift[],
	       const t_pbc *pbc,const t_graph *g,
	       real lambda,real *dvdlambda,
	       const t_mdatoms *md,t_fcdata *fcd,
	       int *global_atom_index)
{
  int  i,m,ki,ai,aj,type,table;
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
    dr   = dr2*gmx_invsqrt(dr2);		        /*  10		*/

    table = forceparams[type].tab.table;

    *dvdlambda += bonded_tab("bond",table,
			     &fcd->bondtab[table],
			     forceparams[type].tab.kA,
			     forceparams[type].tab.kB,
			     dr,lambda,&vbond,&fbond);  /*  22 */

    if (dr2 == 0.0)
      continue;

    
    vtot  += vbond;/* 1*/
    fbond *= gmx_invsqrt(dr2);			/*   6		*/
#ifdef DEBUG
    if (debug)
      fprintf(debug,"TABBONDS: dr = %10g  vbond = %10g  fbond = %10g\n",
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
  }					/* 62 TOTAL	*/
  return vtot;
}

real tab_angles(int nbonds,
		const t_iatom forceatoms[],const t_iparams forceparams[],
		const rvec x[],rvec f[],rvec fshift[],
		const t_pbc *pbc,const t_graph *g,
		real lambda,real *dvdlambda,
		const t_mdatoms *md,t_fcdata *fcd,
		int *global_atom_index)
{
  int  i,ai,aj,ak,t1,t2,type,table;
  rvec r_ij,r_kj;
  real cos_theta,cos_theta2,theta,dVdt,va,vtot;
  ivec jt,dt_ij,dt_kj;
  
  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    
    theta  = bond_angle(x[ai],x[aj],x[ak],pbc,
			r_ij,r_kj,&cos_theta,&t1,&t2);	/*  41		*/

    table = forceparams[type].tab.table;
  
    *dvdlambda += bonded_tab("angle",table,
			     &fcd->angletab[table],
			     forceparams[type].tab.kA,
			     forceparams[type].tab.kB,
			     theta,lambda,&va,&dVdt);  /*  22  */
    vtot += va;
    
    cos_theta2 = sqr(cos_theta);                /*   1		*/
    if (cos_theta2 < 1) {
      int  m;
      real snt,st,sth;
      real cik,cii,ckk;
      real nrkj2,nrij2;
      rvec f_i,f_j,f_k;
      
      st  = dVdt*gmx_invsqrt(1 - cos_theta2);	/*  12		*/
      sth = st*cos_theta;			/*   1		*/
#ifdef DEBUG
      if (debug)
	fprintf(debug,"ANGLES: theta = %10g  vth = %10g  dV/dtheta = %10g\n",
		theta*RAD2DEG,va,dVdt);
#endif
      nrkj2=iprod(r_kj,r_kj);			/*   5		*/
      nrij2=iprod(r_ij,r_ij);
      
      cik=st*gmx_invsqrt(nrkj2*nrij2);		/*  12		*/ 
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
    }                                           /* 169 TOTAL	*/
  }
  return vtot;
}

real tab_dihs(int nbonds,
	      const t_iatom forceatoms[],const t_iparams forceparams[],
	      const rvec x[],rvec f[],rvec fshift[],
	      const t_pbc *pbc,const t_graph *g,
	      real lambda,real *dvdlambda,
	      const t_mdatoms *md,t_fcdata *fcd,
	      int *global_atom_index)
{
  int  i,type,ai,aj,ak,al,table;
  int  t1,t2,t3;
  rvec r_ij,r_kj,r_kl,m,n;
  real phi,sign,ddphi,vpd,vtot;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    al   = forceatoms[i++];
    
    phi=dih_angle(x[ai],x[aj],x[ak],x[al],pbc,r_ij,r_kj,r_kl,m,n,
                  &sign,&t1,&t2,&t3);			/*  84  */

    table = forceparams[type].tab.table;

    /* Hopefully phi+M_PI never results in values < 0 */
    *dvdlambda += bonded_tab("dihedral",table,
			     &fcd->dihtab[table],
			     forceparams[type].tab.kA,
			     forceparams[type].tab.kB,
			     phi+M_PI,lambda,&vpd,&ddphi);
		       
    vtot += vpd;
    do_dih_fup(ai,aj,ak,al,-ddphi,r_ij,r_kj,r_kl,m,n,
	       f,fshift,pbc,g,x,t1,t2,t3);			/* 112	*/

#ifdef DEBUG
    fprintf(debug,"pdih: (%d,%d,%d,%d) phi=%g\n",
	    ai,aj,ak,al,phi);
#endif
  } /* 227 TOTAL 	*/

  return vtot;
}

void calc_bonds(FILE *fplog,const gmx_multisim_t *ms,
		const t_idef *idef,
		rvec x[],history_t *hist,
		rvec f[],t_forcerec *fr,
		const t_pbc *pbc,const t_graph *g,
		gmx_enerdata_t *enerd,t_nrnb *nrnb,
		real lambda,
		const t_mdatoms *md,
		t_fcdata *fcd,int *global_atom_index,
		t_atomtypes *atype, gmx_genborn_t *born,
		gmx_bool bPrintSepPot,gmx_large_int_t step)
{
  int    ftype,nbonds,ind,nat1;
  real   *epot,v,dvdl;
  const  t_pbc *pbc_null;
  char   buf[22];

  if (fr->bMolPBC)
    pbc_null = pbc;
  else
    pbc_null = NULL;

  if (bPrintSepPot)
    fprintf(fplog,"Step %s: bonded V and dVdl for this node\n",
	    gmx_step_str(step,buf));

#ifdef DEBUG
  if (g && debug)
    p_graph(debug,"Bondage is fun",g);
#endif
  
  epot = enerd->term;

  /* Do pre force calculation stuff which might require communication */
  if (idef->il[F_ORIRES].nr) {
    epot[F_ORIRESDEV] = calc_orires_dev(ms,idef->il[F_ORIRES].nr,
					idef->il[F_ORIRES].iatoms,
					idef->iparams,md,(const rvec*)x,
					pbc_null,fcd,hist);
  }
  if (idef->il[F_DISRES].nr) {
    calc_disres_R_6(ms,idef->il[F_DISRES].nr,
		    idef->il[F_DISRES].iatoms,
		    idef->iparams,(const rvec*)x,pbc_null,
		    fcd,hist);
  }
  
  /* Loop over all bonded force types to calculate the bonded forces */
  for(ftype=0; (ftype<F_NRE); ftype++) {
	  if(ftype<F_GB12 || ftype>F_GB14) {
    if (interaction_function[ftype].flags & IF_BOND &&
	!(ftype == F_CONNBONDS || ftype == F_POSRES)) {
      nbonds=idef->il[ftype].nr;
      if (nbonds > 0) {
	ind  = interaction_function[ftype].nrnb_ind;
	nat1 = interaction_function[ftype].nratoms + 1;
	dvdl = 0;
	if (ftype < F_LJ14 || ftype > F_LJC_PAIRS_NB) {
		if(ftype==F_CMAP)
		{
			v = cmap_dihs(nbonds,idef->il[ftype].iatoms,
						  idef->iparams,&idef->cmap_grid,
						  (const rvec*)x,f,fr->fshift,
						  pbc_null,g,lambda,&dvdl,md,fcd,
						  global_atom_index);
		}
		else
		{
			v =
	    interaction_function[ftype].ifunc(nbonds,idef->il[ftype].iatoms,
					      idef->iparams,
					      (const rvec*)x,f,fr->fshift,
					      pbc_null,g,lambda,&dvdl,md,fcd,
					      global_atom_index);
		}

	  if (bPrintSepPot) {
	    fprintf(fplog,"  %-23s #%4d  V %12.5e  dVdl %12.5e\n",
		    interaction_function[ftype].longname,nbonds/nat1,v,dvdl);
	  }
	} else {
	  v = do_listed_vdw_q(ftype,nbonds,idef->il[ftype].iatoms,
			      idef->iparams,
			      (const rvec*)x,f,fr->fshift,
			      pbc_null,g,
			      lambda,&dvdl,
			      md,fr,&enerd->grpp,global_atom_index);
	  if (bPrintSepPot) {
	    fprintf(fplog,"  %-5s + %-15s #%4d                  dVdl %12.5e\n",
		    interaction_function[ftype].longname,
		    interaction_function[F_COUL14].longname,nbonds/nat1,dvdl);
	  }
	}
	if (ind != -1)
	  inc_nrnb(nrnb,ind,nbonds/nat1);
	epot[ftype]        += v;
	enerd->dvdl_nonlin += dvdl;
      }
    }
  }
  }
  /* Copy the sum of violations for the distance restraints from fcd */
  if (fcd)
    epot[F_DISRESVIOL] = fcd->disres.sumviol;
}

void calc_bonds_lambda(FILE *fplog,
		       const t_idef *idef,
		       rvec x[],
		       t_forcerec *fr,
		       const t_pbc *pbc,const t_graph *g,
		       gmx_enerdata_t *enerd,t_nrnb *nrnb,
		       real lambda,
		       const t_mdatoms *md,
		       t_fcdata *fcd,int *global_atom_index)
{
    int    ftype,nbonds_np,nbonds,ind, nat1;
  real   *epot,v,dvdl;
  rvec   *f,*fshift_orig;
  const  t_pbc *pbc_null;
  t_iatom *iatom_fe;

  if (fr->bMolPBC)
    pbc_null = pbc;
  else
    pbc_null = NULL;
  
  epot = enerd->term;
  
  snew(f,fr->natoms_force);
  /* We want to preserve the fshift array in forcerec */
  fshift_orig = fr->fshift;
  snew(fr->fshift,SHIFTS);

  /* Loop over all bonded force types to calculate the bonded forces */
  for(ftype=0; (ftype<F_NRE); ftype++) {
      if(ftype<F_GB12 || ftype>F_GB14) {
          
          if (interaction_function[ftype].flags & IF_BOND &&
              !(ftype == F_CONNBONDS || ftype == F_POSRES)) 
          {
              nbonds_np = idef->il[ftype].nr_nonperturbed;
              nbonds    = idef->il[ftype].nr - nbonds_np;
              nat1 = interaction_function[ftype].nratoms + 1;
              if (nbonds > 0) {
                  ind  = interaction_function[ftype].nrnb_ind;
                  iatom_fe = idef->il[ftype].iatoms + nbonds_np;
                  dvdl = 0;
                  if (ftype < F_LJ14 || ftype > F_LJC_PAIRS_NB) {
                      v =
                          interaction_function[ftype].ifunc(nbonds,iatom_fe,
                                                            idef->iparams,
                                                            (const rvec*)x,f,fr->fshift,
                                                            pbc_null,g,lambda,&dvdl,md,fcd,
                                                            global_atom_index);
                  } else {
                      v = do_listed_vdw_q(ftype,nbonds,iatom_fe,
                                          idef->iparams,
                                          (const rvec*)x,f,fr->fshift,
                                          pbc_null,g,
                                          lambda,&dvdl,
                                          md,fr,&enerd->grpp,global_atom_index);
                  }
                  if (ind != -1)
                      inc_nrnb(nrnb,ind,nbonds/nat1);
                  epot[ftype] += v;
              }
          }
      }
  }

  sfree(fr->fshift);
  fr->fshift = fshift_orig;
  sfree(f);
}
