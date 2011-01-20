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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include "maths.h"
#include "typedefs.h"
#include "vec.h"
#include "coulomb.h"
#include "smalloc.h"
#include "physics.h"
#include "txtdump.h"
#include "futil.h"
#include "names.h"
#include "writeps.h"
#include "macros.h"

real calc_ewaldcoeff(real rc,real dtol)
{
  real x=5,low,high;
  int n,i=0;
  
  
  do {
    i++;
    x*=2;
  } while (gmx_erfc(x*rc) > dtol);

  n=i+60; /* search tolerance is 2^-60 */
  low=0;
  high=x;
  for(i=0;i<n;i++) {
    x=(low+high)/2;
    if (gmx_erfc(x*rc) > dtol)
      low=x;
    else 
      high=x;
  }
  return x;
}



real ewald_LRcorrection(FILE *fplog,
			int start,int end,
			t_commrec *cr,t_forcerec *fr,
			real *chargeA,real *chargeB,
			t_blocka *excl,rvec x[],
			matrix box,rvec mu_tot[],
			int ewald_geometry,real epsilon_surface,
			real lambda,real *dvdlambda,
			real *vdip,real *vcharge)
{
  int     i,i1,i2,j,k,m,iv,jv,q;
  atom_id *AA;
  double  q2sumA,q2sumB,Vexcl,dvdl_excl; /* Necessary for precision */
  real    one_4pi_eps;
  real    v,vc,qiA,qiB,dr,dr2,rinv,fscal,enercorr;
  real    VselfA,VselfB=0,Vcharge[2],Vdipole[2],rinv2,ewc=fr->ewaldcoeff,ewcdr;
  rvec    df,dx,mutot[2],dipcorrA,dipcorrB;
  rvec    *f=fr->f_novirsum;
  tensor  dxdf;
  real    vol = box[XX][XX]*box[YY][YY]*box[ZZ][ZZ];
  real    L1,dipole_coeff,qqA,qqB,qqL,vr0;
  /*#define TABLES*/
#ifdef TABLES
  real    tabscale=fr->tabscale;
  real    eps,eps2,VV,FF,F,Y,Geps,Heps2,Fp,fijC,r1t;
  real    *VFtab=fr->coulvdwtab;
  int     n0,n1,nnn;
#else
  double  isp=0.564189583547756;
#endif
  int     niat;
  gmx_bool    bFreeEnergy = (chargeB != NULL);
  gmx_bool    bMolPBC = fr->bMolPBC;

  one_4pi_eps = ONE_4PI_EPS0/fr->epsilon_r;
  vr0 = ewc*2/sqrt(M_PI);

  AA         = excl->a;
  Vexcl      = 0;
  dvdl_excl  = 0;
  q2sumA     = 0;
  q2sumB     = 0;
  Vdipole[0] = 0;
  Vdipole[1] = 0;
  Vcharge[0] = 0;
  Vcharge[1] = 0;
  L1         = 1.0-lambda;

  /* Note that we have to transform back to gromacs units, since
   * mu_tot contains the dipole in debye units (for output).
   */
  for(i=0; (i<DIM); i++) {
    mutot[0][i] = mu_tot[0][i]*DEBYE2ENM;
    mutot[1][i] = mu_tot[1][i]*DEBYE2ENM;
    dipcorrA[i] = 0;
    dipcorrB[i] = 0;
  }
  dipole_coeff=0;
  switch (ewald_geometry) {
  case eewg3D:
    if (epsilon_surface != 0) {
      dipole_coeff =
	2*M_PI*ONE_4PI_EPS0/((2*epsilon_surface + fr->epsilon_r)*vol);
      for(i=0; (i<DIM); i++) {
	dipcorrA[i] = 2*dipole_coeff*mutot[0][i];
	dipcorrB[i] = 2*dipole_coeff*mutot[1][i];
      }
    }
    break;
  case eewg3DC:
    dipole_coeff = 2*M_PI*one_4pi_eps/vol;
    dipcorrA[ZZ] = 2*dipole_coeff*mutot[0][ZZ];
    dipcorrB[ZZ] = 2*dipole_coeff*mutot[1][ZZ];
    break;
  default:
    gmx_incons("Unsupported Ewald geometry");
    break;
  }
  if (debug) {
    fprintf(debug,"dipcorr = %8.3f  %8.3f  %8.3f\n",
	    dipcorrA[XX],dipcorrA[YY],dipcorrA[ZZ]);
    fprintf(debug,"mutot   = %8.3f  %8.3f  %8.3f\n",
	    mutot[0][XX],mutot[0][YY],mutot[0][ZZ]);
  }

  if (DOMAINDECOMP(cr))
    niat = excl->nr;
  else
    niat = end; 
      
  clear_mat(dxdf);
  if (!bFreeEnergy) {
    for(i=start; (i<niat); i++) {
      /* Initiate local variables (for this i-particle) to 0 */
      qiA = chargeA[i]*one_4pi_eps;
      i1  = excl->index[i];
      i2  = excl->index[i+1];
      if (i < end)
	q2sumA += chargeA[i]*chargeA[i];
      
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
	  qqA = qiA*chargeA[k];
	  if (qqA != 0.0) {
	    rvec_sub(x[i],x[k],dx);
	    if (bMolPBC) {
	      /* Cheap pbc_dx, assume excluded pairs are at short distance. */
	      for(m=DIM-1; (m>=0); m--) {
		if (dx[m] > 0.5*box[m][m])
		  rvec_dec(dx,box[m]);
		else if (dx[m] < -0.5*box[m][m])
		  rvec_inc(dx,box[m]);
	      }
	    }
	    dr2 = norm2(dx);
	    /* Distance between two excluded particles may be zero in the
	     * case of shells
	     */
	    if (dr2 != 0) {
	      rinv              = gmx_invsqrt(dr2);
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
	      vc                = qqA*(rinv-VV);
	      fijC              = qqA*FF;
	      Vexcl            += vc;
	      
	      fscal             = vc*rinv2+fijC*tabscale*rinv;
	    /* End of tabulated interaction part */
#else
	    
	    /* This is the code you would want instead if not using
	     * tables: 
	     */
	      ewcdr   = ewc*dr;
	      vc      = qqA*gmx_erf(ewcdr)*rinv;
	      Vexcl  += vc;
#ifdef GMX_DOUBLE
	      /* Relative accuracy at R_ERF_R_INACC of 3e-10 */
#define	      R_ERF_R_INACC 0.006
#else
              /* Relative accuracy at R_ERF_R_INACC of 2e-5 */
#define	      R_ERF_R_INACC 0.1
#endif
	      if (ewcdr > R_ERF_R_INACC) {
		fscal = rinv2*(vc - 2.0*qqA*ewc*isp*exp(-ewcdr*ewcdr));
	      } else {
		/* Use a fourth order series expansion for small ewcdr */
		fscal = ewc*ewc*qqA*vr0*(2.0/3.0 - 0.4*ewcdr*ewcdr);
	      }
#endif
	      /* The force vector is obtained by multiplication with the 
	       * distance vector 
	       */
	      svmul(fscal,dx,df);
	      rvec_inc(f[k],df);
	      rvec_dec(f[i],df);
	      for(iv=0; (iv<DIM); iv++)
		for(jv=0; (jv<DIM); jv++)
		  dxdf[iv][jv] += dx[iv]*df[jv];
	    } else {
	      Vexcl += qqA*vr0;
	    }
	  }
	}
      }
      /* Dipole correction on force */
      if (dipole_coeff != 0) {
	for(j=0; (j<DIM); j++)
	  f[i][j] -= dipcorrA[j]*chargeA[i];
      }
    }
  } else {
    for(i=start; (i<niat); i++) {
      /* Initiate local variables (for this i-particle) to 0 */
      qiA = chargeA[i]*one_4pi_eps;
      qiB = chargeB[i]*one_4pi_eps;
      i1  = excl->index[i];
      i2  = excl->index[i+1];
      if (i < end) {
	q2sumA += chargeA[i]*chargeA[i];
	q2sumB += chargeB[i]*chargeB[i];
      }
      
      /* Loop over excluded neighbours */
      for(j=i1; (j<i2); j++) {
	k = AA[j];
	if (k > i) {
	  qqA = qiA*chargeA[k];
	  qqB = qiB*chargeB[k];
	  if (qqA != 0.0 || qqB != 0.0) {
	    qqL = L1*qqA + lambda*qqB;
	    rvec_sub(x[i],x[k],dx);
	    if (bMolPBC) {
	      /* Cheap pbc_dx, assume excluded pairs are at short distance. */
	      for(m=DIM-1; (m>=0); m--) {
		if (dx[m] > 0.5*box[m][m])
		  rvec_dec(dx,box[m]);
		else if (dx[m] < -0.5*box[m][m])
		  rvec_inc(dx,box[m]);
	      }
	    }
	    dr2 = norm2(dx);
	    if (dr2 != 0) {
	      rinv   = gmx_invsqrt(dr2);
	      rinv2  = rinv*rinv;
	      dr     = 1.0/rinv;      
	      v      = gmx_erf(ewc*dr)*rinv;
	      vc     = qqL*v;
	      Vexcl += vc;
	      fscal  = rinv2*(vc-2.0*qqL*ewc*isp*exp(-ewc*ewc*dr2));
	      svmul(fscal,dx,df);
	      rvec_inc(f[k],df);
	      rvec_dec(f[i],df);
	      for(iv=0; (iv<DIM); iv++)
		for(jv=0; (jv<DIM); jv++)
		  dxdf[iv][jv] += dx[iv]*df[jv];
	      dvdl_excl += (qqB - qqA)*v;
	    } else {
	      Vexcl     +=         qqL*vr0;
	      dvdl_excl += (qqB - qqA)*vr0;
	    }
	  }
	}
      }
      /* Dipole correction on force */
      if (dipole_coeff != 0) {
	for(j=0; (j<DIM); j++)
	  f[i][j] -= L1*dipcorrA[j]*chargeA[i]
	    + lambda*dipcorrB[j]*chargeB[i];
      } 
    }
  }
  for(iv=0; (iv<DIM); iv++)
    for(jv=0; (jv<DIM); jv++)
      fr->vir_el_recip[iv][jv] += 0.5*dxdf[iv][jv];
      
  /* Global corrections only on master process */
  if (MASTER(cr)) {
    for(q=0; q<(bFreeEnergy ? 2 : 1); q++) {
      /* Apply charge correction */
      /* use vc as a dummy variable */
      vc = fr->qsum[q]*fr->qsum[q]*M_PI*one_4pi_eps/(2.0*vol*vol*ewc*ewc);
      for(iv=0; (iv<DIM); iv++)
	fr->vir_el_recip[iv][iv] +=
	  (bFreeEnergy ? (q==0 ? L1*vc : lambda*vc) : vc);
      Vcharge[q] = -vol*vc;
      
      /* Apply surface dipole correction:
       * correction = dipole_coeff * (dipole)^2
       */
      if (dipole_coeff != 0) {
	if (ewald_geometry == eewg3D)
	  Vdipole[q] = dipole_coeff*iprod(mutot[q],mutot[q]);
	else if (ewald_geometry == eewg3DC)
	  Vdipole[q] = dipole_coeff*mutot[q][ZZ]*mutot[q][ZZ];
      }
    }
  }    
  
  VselfA = ewc*one_4pi_eps*q2sumA/sqrt(M_PI);

  if (!bFreeEnergy) {
    *vcharge = Vcharge[0];
    *vdip    = Vdipole[0];
    enercorr = *vcharge + *vdip - VselfA - Vexcl;
   } else {
    VselfB = ewc*one_4pi_eps*q2sumB/sqrt(M_PI);
    *vcharge = L1*Vcharge[0] + lambda*Vcharge[1];
    *vdip    = L1*Vdipole[0] + lambda*Vdipole[1];
    enercorr = *vcharge + *vdip - (L1*VselfA + lambda*VselfB) - Vexcl;
    *dvdlambda += Vdipole[1] + Vcharge[1] - VselfB
      - (Vdipole[0] + Vcharge[0] - VselfA) - dvdl_excl;
  }

  if (debug) {
    fprintf(debug,"Long Range corrections for Ewald interactions:\n");
    fprintf(debug,"start=%d,natoms=%d\n",start,end-start);
    fprintf(debug,"q2sum = %g, Vself=%g\n",
	    L1*q2sumA+lambda*q2sumB,L1*VselfA+lambda*VselfB);
    fprintf(debug,"Long Range correction: Vexcl=%g\n",Vexcl);
    if (MASTER(cr)) {
      fprintf(debug,"Total charge correction: Vcharge=%g\n",
	      L1*Vcharge[0]+lambda*Vcharge[1]);
      if (epsilon_surface > 0 || ewald_geometry == eewg3DC) {
	fprintf(debug,"Total dipole correction: Vdipole=%g\n",
		L1*Vdipole[0]+lambda*Vdipole[1]);
      }
    }
  }
    
  /* Return the correction to the energy */
  return enercorr;
}
