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

#include <stdio.h>
#include <math.h>
#include "assert.h"
#include "typedefs.h"
#include "vec.h"
#include "ewald_util.h"
#include "smalloc.h"
#include "physics.h"
#include "txtdump.h"
#include "futil.h"
#include "names.h"
#include "fftgrid.h"
#include "writeps.h"
#include "macros.h"
#include "xvgr.h"


real calc_ewaldcoeff(real rc,real dtol)
{
  real x=5,low,high;
  int n,i=0;
  
  
  do {
    i++;
    x*=2;
  } while (erfc(x*rc) > dtol);

  n=i+60; /* search tolerance is 2^-60 */
  low=0;
  high=x;
  for(i=0;i<n;i++) {
    x=(low+high)/2;
    if (erfc(x*rc) > dtol)
      low=x;
    else 
      high=x;
  }
  return x;
}



real ewald_LRcorrection(FILE *fp,t_nsborder *nsb,t_commrec *cr,t_forcerec *fr,
			real charge[],t_block *excl,rvec x[],
			matrix box,rvec mu_tot,real qsum, int ewald_geometry,
			real epsilon_surface,matrix lr_vir)
{
  real Vself;
  int     i,i1,i2,j,k,m,iv,jv;
  atom_id *AA;
  double  q2sum,Vexcl; /* Necessary for precision */
  real    vc,qi,dr,dr2,rinv,fscal,Vcharge,Vdipole,rinv2,ewc=fr->ewaldcoeff;
  rvec    df,dx,mutot,dipcorr;
  rvec    *f=fr->f_el_recip;
  real    vol = box[XX][XX]*box[YY][YY]*box[ZZ][ZZ];
  real    dipole_coeff,qq;
  /*#define TABLES*/
#ifdef TABLES
  real    tabscale=fr->tabscale;
  real    eps,eps2,VV,FF,F,Y,Geps,Heps2,Fp,fijC,r1t;
  real    *VFtab=fr->coulvdwtab;
  int     n0,n1,nnn;
#else
  double  isp=0.564189583547756;
#endif
  int     start = START(nsb);
  int     end   = start+HOMENR(nsb);
  
  AA      = excl->a;
  Vexcl   = 0;
  q2sum   = 0; 
  Vdipole = 0;
  Vcharge = 0;

  if (epsilon_surface == 0)
    dipole_coeff=0;
  else { 
    switch (ewald_geometry) {
    case eewg3D:
      dipole_coeff = 2*M_PI*ONE_4PI_EPS0/((2*epsilon_surface+1)*vol);
      break;
    case eewg3DC:
      dipole_coeff = 2*M_PI*ONE_4PI_EPS0/vol;
      break;
    default:
      fatal_error(0,"Unsupported Ewald geometry");
      dipole_coeff=0;
      break;
    }
  }
  /* Note that we have to transform back to gromacs units, since
   * mu_tot contains the dipole in debye units (for output).
   */
  for(i=0; (i<DIM); i++) {
    mutot[i]   = mu_tot[i]*DEBYE2ENM;
    dipcorr[i] = 2.0*dipole_coeff*mutot[i];
  }
  if (debug) {
    fprintf(debug,"dipcorr = %8.3f  %8.3f  %8.3f\n",
	    dipcorr[XX],dipcorr[YY],dipcorr[ZZ]);
    fprintf(debug,"mutot   = %8.3f  %8.3f  %8.3f\n",
	    mutot[XX],mutot[YY],mutot[ZZ]);
  }
  for(i=start; (i<end); i++) {
    /* Initiate local variables (for this i-particle) to 0 */
    qi  = charge[i]*ONE_4PI_EPS0;
    i1  = excl->index[i];
    i2  = excl->index[i+1];
    q2sum += charge[i]*charge[i];
    
    /* Loop over excluded neighbours */
    for(j=i1; (j<i2); j++) {
      k = AA[j];
      /* 
       * First we must test whether k <> i, and then, because the
       * exclusions are all listed twice i->k and k->i we must select
       * just one of the two.
       * As a minor optimization we only compute forces when the charges
       * are non-zero.
       */
      if (k > i) {
	qq = qi*charge[k];
	if (qq != 0.0) {
	  dr2 = 0;
	  rvec_sub(x[i],x[k],dx);
	  for(m=DIM-1; (m>=0); m--) {
	    if (dx[m] > 0.5*box[m][m])
	      rvec_dec(dx,box[m]);
	    else if (dx[m] < -0.5*box[m][m])
	      rvec_inc(dx,box[m]);
	    
	    dr2  += dx[m]*dx[m];
	  }

	  /* It might be possible to optimize this slightly 
	   * when using spc or similar water models:
	   * Potential could be stored once, at the beginning,
	   * and forces could use that bonds are constant 
	   */
	  /* The Ewald interactions are tabulated. If you 
	   * remove the tabulation you must replace this
	   * part too
	   */

	  rinv              = invsqrt(dr2);
	  rinv2             = rinv*rinv;
	  dr                = 1.0/rinv;      
#ifdef TABLES
	  r1t               = tabscale*dr;   
	  n0                = r1t;
	  assert(n0 >= 3);
	  n1                = 12*n0;
	  eps               = r1t-n0;
	  eps2              = eps*eps;
	  nnn               = n1;
	  Y                 = VFtab[nnn];
	  F                 = VFtab[nnn+1];
	  Geps              = eps*VFtab[nnn+2];
	  Heps2             = eps2*VFtab[nnn+3];
	  Fp                = F+Geps+Heps2;
	  VV                = Y+eps*Fp;
	  FF                = Fp+Geps+2.0*Heps2;
	  vc                = qq*(rinv-VV);
	  fijC              = qq*FF;
	  Vexcl            += vc;
	  
	  fscal             = vc*rinv2+fijC*tabscale*rinv;
	  /* End of tabulated interaction part */
#else

	  /* This is the code you would want instead if not using
	   * tables: 
	   */
	  vc      = qq*erf(ewc*dr)*rinv;
	  Vexcl  += vc;
	  fscal   = rinv2*(vc-2.0*qq*ewc*isp*exp(-ewc*ewc*dr2));
#endif
	  
	  /* The force vector is obtained by multiplication with the 
	   * distance vector 
	   */
	  svmul(fscal,dx,df);
	  if (debug)
	    fprintf(debug,"dr=%8.4f, fscal=%8.0f, df=%10.0f,%10.0f,%10.0f\n",
		    dr,fscal,df[XX],df[YY],df[ZZ]);
	  rvec_inc(f[k],df);
	  rvec_dec(f[i],df);
	  for(iv=0; (iv<DIM); iv++)
	    for(jv=0; (jv<DIM); jv++)
	      lr_vir[iv][jv]+=0.5*dx[iv]*df[jv];
	}
      }
    }
    /* Dipole correction on force */
    if (dipole_coeff != 0) {
      if (ewald_geometry == eewg3D) {
	for(j=0; (j<DIM); j++)
	  f[i][j] -= dipcorr[j]*charge[i];
      } 
      else if (ewald_geometry == eewg3DC) {
	f[i][ZZ] -= dipcorr[ZZ]*charge[i];
      }
    }
  }
  
  /* Global corrections only on master process */
  if (MASTER(cr)) {
    /* Apply charge correction */
    /* use vc as a dummy variable */
    vc = qsum*qsum*M_PI*ONE_4PI_EPS0/(2.0*vol*vol*ewc*ewc);
    for(iv=0; (iv<DIM); iv++)
      lr_vir[iv][iv] += vc;
    Vcharge = -vol*vc;
    
    /* Apply surface dipole correction:
     * correction = dipole_coeff * (dipole)^2
     */
    if (dipole_coeff != 0) {
      if (ewald_geometry == eewg3D) 
	Vdipole = dipole_coeff*iprod(mutot,mutot);
      else if (ewald_geometry == eewg3DC) 
	Vdipole = dipole_coeff*mutot[ZZ]*mutot[ZZ];
    }
  }
  
  Vself=ewc*ONE_4PI_EPS0*q2sum/sqrt(M_PI);
  
  if (debug) {
    fprintf(debug,"Long Range corrections for Ewald interactions:\n");
    fprintf(debug,"start=%d,natoms=%d\n",start,end-start);
    fprintf(debug,"q2sum = %g, Vself=%g\n",q2sum, Vself);
    fprintf(debug,"Long Range correction: Vexcl=%g\n",Vexcl);
    if (MASTER(cr)) {
      fprintf(debug,"Total charge correction: Vcharge=%g\n",Vcharge);
      if (epsilon_surface>0)
	fprintf(debug,"Total dipole correction: Vdipole=%g\n",Vdipole);
    }
  }

  /* Return the correction to the energy */
  return (Vdipole+Vcharge-Vself-Vexcl);
}
  

