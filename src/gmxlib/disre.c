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
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "typedefs.h"
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "futil.h"
#include "xvgr.h"
#include "fatal.h"
#include "bondf.h"
#include "copyrite.h"
#include "disre.h"
#include "main.h"

void init_disres(FILE *log,int nfa,t_iatom forceatoms[],t_iparams ip[],
		 t_inputrec *ir,t_commrec *mcr,t_fcdata *fcd)
{
  int          fa;
  t_disresdata *dd;

  dd = &(fcd->disres);

  dd->dr_weighting = ir->eDisreWeighting;
  dd->dr_fc        = ir->dr_fc;
  dd->dr_tau       = ir->dr_tau;
  if (dd->dr_tau == 0.0) {
    dd->dr_bMixed = FALSE;
    dd->ETerm = 0.0;
  } else {
    dd->dr_bMixed = ir->bDisreMixed;
    dd->ETerm = exp(-(ir->delta_t/ir->dr_tau));
  }
  dd->ETerm1 = 1.0 - dd->ETerm;
  dd->exp_min_t_tau = 1.0;
  
  dd->nr = 0;
  for(fa=0; fa<nfa; fa+=3)
    if (fa==0 || 
	ip[forceatoms[fa-3]].disres.label != ip[forceatoms[fa]].disres.label)
      dd->nr++;
  dd->npr = nfa/(interaction_function[F_DISRES].nratoms+1);
  snew(dd->rt,dd->npr);
  snew(dd->rav,dd->npr);
  /* Allocate Rt_6 and Rav_6 consecutively in memory so they can be
   * averaged over the processors in one call (in calc_disre_R_6)
   */
  snew(dd->Rt_6,2*dd->nr);
  dd->Rav_6 = &(dd->Rt_6[dd->nr]);
  if (mcr)
    snew(dd->Rtl_6,dd->nr);
  else
    dd->Rtl_6 = dd->Rt_6;
  
  if (dd->npr > 0) {
    fprintf(log,"There are %d distance restraints involving %d atom pairs\n",
	    dd->nr,dd->npr);
    if (mcr)
      check_multi_int(log,mcr,fcd->disres.nr,
		      "the number of distance restraints");
    please_cite(log,"Torda89a");
  }
}

void calc_disres_R_6(t_commrec *mcr,
		     int nfa,t_iatom forceatoms[],t_iparams ip[],
		     rvec x[],bool bFullPBC,t_fcdata *fcd)
{
  atom_id     ai,aj;
  int         fa,res,i,pair,ki,kj,m;
  int         type,label;
  rvec        dx;
  real        rt_1,rt_3,rt2,rav_3,*rt,*rav,*Rtl_6,*Rt_6,*Rav_6;
  ivec        it,jt,dt;
  t_disresdata *dd;
  real        ETerm,ETerm1,cf1,cf2,invn=0;

  dd = &(fcd->disres);
  ETerm        = dd->ETerm;
  ETerm1       = dd->ETerm1;
  rt           = dd->rt;
  rav          = dd->rav;
  Rtl_6        = dd->Rtl_6;
  Rt_6         = dd->Rt_6;
  Rav_6        = dd->Rav_6;

  /* scaling factor to smoothly turn on the restraint forces *
   * when using time averaging                               */
  dd->exp_min_t_tau *= ETerm;

  cf1 = dd->exp_min_t_tau;
  cf2 = 1.0/(1.0 - dd->exp_min_t_tau);
  
  if (mcr)
    invn = 1.0/mcr->nnodes;

  /* 'loop' over all atom pairs (pair_nr=fa/3) involved in restraints, *
   * the total number of atoms pairs is nfa/3                          */
  res = 0;
  fa  = 0;
  while (fa < nfa) {
    type  = forceatoms[fa];
    label = ip[type].disres.label;

    Rav_6[res] = 0.0;
    Rt_6[res]  = 0.0;
  
    /* Loop over the atom pairs of 'this' restraint */
    while (fa<nfa && ip[forceatoms[fa]].disres.label==label) {
      pair = fa/3;
      ai   = forceatoms[fa+1];
      aj   = forceatoms[fa+2];

      if (bFullPBC)
	pbc_dx(x[ai],x[aj],dx);
      else
	rvec_sub(x[ai],x[aj],dx);
      rt2  = iprod(dx,dx);
      rt_1 = invsqrt(rt2);
      rt_3 = rt_1*rt_1*rt_1;

      rt[pair]  = sqrt(rt2);
      rav_3     = cf2*((ETerm - cf1)*rav[pair] + ETerm1*rt_3);
      rav[pair] = rav_3;
      
      Rt_6[res]  += rt_3*rt_3;
      Rav_6[res] += rav_3*rav_3;

      fa += 3;
    }
    if (mcr) {
      Rtl_6[res]  = Rt_6[res];
      Rt_6[res]  *= invn;
      Rav_6[res] *= invn;
    }
    res++;
  }
  
  if ((mcr) && PAR(mcr))
    gmx_sum(2*dd->nr,Rt_6,mcr);
}

real ta_disres(int nfa,t_iatom forceatoms[],t_iparams ip[],
	       rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	       matrix box,real lambda,real *dvdlambda,
	       t_mdatoms *md,int ngrp,real egnb[],real egcoul[],
	       t_fcdata *fcd)
{
  const real sixth=1.0/6.0;
  const real seven_three=7.0/3.0;
  
  atom_id     ai,aj;
  int         fa,res,npairs,p,pair,ki=CENTRAL,m;
  int         type,label;
  rvec        dx;
  real        weight_rt_1;
  rvec        *fshift;
  real        smooth_fc,Rt,Rav,rt2,*Rtl_6,*Rt_6,*Rav_6;
  real        k0,f_scal=0,fmax_scal,fk_scal,fij;
  real        tav_viol,instant_viol,mixed_viol,violtot,vtot;
  real        tav_viol_Rav7,instant_viol_Rav7;
  real        up1,up2,low;
  bool        bConservative,bMixed,bViolation;
  ivec        it,jt,dt;
  t_disresdata *dd;
  int         dr_weighting;
  bool        dr_bMixed,bFullPBC;
  real        dr_fc;

  bFullPBC = (fr->ePBC == epbcFULL);

  dd = &(fcd->disres);
  dr_weighting = dd->dr_weighting;
  dr_bMixed    = dd->dr_bMixed;
  dr_fc        = dd->dr_fc;
  Rtl_6        = dd->Rtl_6;
  Rt_6         = dd->Rt_6;
  Rav_6        = dd->Rav_6;

  tav_viol=instant_viol=mixed_viol=tav_viol_Rav7=instant_viol_Rav7=0;
  /* scaling factor to smoothly turn on the restraint forces *
   * when using time averaging                               */
  smooth_fc = dr_fc * (1.0 - dd->exp_min_t_tau); 
  
  fshift  = fr->fshift; 
  violtot = 0;
  vtot    = 0;
  
  /* 'loop' over all atom pairs (pair_nr=fa/3) involved in restraints, *
   * the total number of atoms pairs is nfa/3                          */
  res  = 0;
  fa   = 0;
  while (fa < nfa) {
    type  = forceatoms[fa];
    /* Take action depending on restraint, calculate scalar force */
    label = ip[type].disres.label;
    up1   = ip[type].disres.up1;
    up2   = ip[type].disres.up2;
    low   = ip[type].disres.low;
    k0    = smooth_fc*ip[type].disres.kfac;

    /* Count the number of pairs in this restraint */
    npairs=1;
    while(fa + 3*npairs < nfa &&
	  ip[forceatoms[fa + 3*npairs]].disres.label == label)
      npairs++;

    /* save some flops when there is only one pair */
    if (ip[type].disres.type != 2) {
      bConservative = (dr_weighting == edrwConservative) && (npairs>1);
      bMixed        = dr_bMixed;
      Rt  = pow(Rt_6[res],-sixth);
      Rav = pow(Rav_6[res],-sixth);
    } else {
      /* When rtype=2 use the instantaneous not ensemble avereged distance */
      bConservative = npairs>1;
      bMixed        = FALSE;
      Rt  = pow(Rtl_6[res],-sixth);
      Rav = Rt;
    }
    
    if (Rav > up1) {
      bViolation = TRUE;
      tav_viol = Rav - up1;
    } else if (Rav < low) {
      bViolation = TRUE;
      tav_viol = Rav - low;
    } else
      bViolation = FALSE;

    if (bViolation) {
      /* NOTE: there is no real potential when time averaging is applied */
      vtot += 0.5*k0*sqr(tav_viol);
      if (1/vtot == 0)
	printf("vtot is inf: %f\n",vtot);
      if (!bMixed) {
	f_scal   = -k0*tav_viol;
	violtot += fabs(tav_viol);
      } else {
	if (Rt > up1) {
	  if (tav_viol > 0)
	    instant_viol = Rt - up1;
	  else
	    bViolation = FALSE;
	} else if (Rt < low) {
	  if (tav_viol < 0)
	    instant_viol = Rt - low;
	  else
	    bViolation = FALSE;
	} else
	  bViolation = FALSE;
	if (bViolation) {
	  mixed_viol = sqrt(tav_viol*instant_viol);
	  f_scal     = -k0*mixed_viol;
	  violtot   += mixed_viol;
	}
      }
    }
    /* save some flops and don't use variables when they are not calculated */
    if (bViolation) {
      fmax_scal = -k0*(up2-up1);
      /* Correct the force for the number of restraints */
      if (bConservative) {
	f_scal  = max(f_scal,fmax_scal);
	if (!bMixed)
	  f_scal *= Rav/Rav_6[res];
	else {
	  f_scal /= 2*mixed_viol;
	  tav_viol_Rav7     = tav_viol*Rav/Rav_6[res];
	  instant_viol_Rav7 = instant_viol*Rt/Rt_6[res];
	}
      } else {
	f_scal /= (real)npairs;
	f_scal  = max(f_scal,fmax_scal);
      }    
      
      /* Exert the force ... */
      
      /* Loop over the atom pairs of 'this' restraint */
      for(p=0; p<npairs; p++) {
	pair = fa/3;
	ai   = forceatoms[fa+1];
	aj   = forceatoms[fa+2];

	if (bFullPBC) 
	  ki = pbc_dx(x[ai],x[aj],dx);
	else
	  rvec_sub(x[ai],x[aj],dx);
	rt2 = iprod(dx,dx);
	
	weight_rt_1 = invsqrt(rt2);
	
	if (bConservative) {
	  if (!dr_bMixed)
	    weight_rt_1 *= pow(dd->rav[pair],seven_three);
	  else
	    weight_rt_1 *= tav_viol_Rav7*pow(dd->rav[pair],seven_three)+
	    instant_viol_Rav7*pow(dd->rt[pair],-7);
	}
	
	fk_scal  = f_scal*weight_rt_1;

	if (g) {
	  ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
	  ki=IVEC2IS(dt);
	}

	for(m=0; m<DIM; m++) {
	  fij            = fk_scal*dx[m];
	  
	  f[ai][m]      += fij;
	  f[aj][m]      -= fij;
	  fshift[ki][m] += fij;
	  fshift[CENTRAL][m] -= fij;
	}
	fa += 3;
      }
    } else
      fa += 3*npairs;
    res++;
  }
  
  dd->sumviol = violtot;

  /* Return energy */
  return vtot;
}
