/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_disre_c = "$Id$";

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

/* Some local variables */
static int  dr_weighting; /* Weighting of pairs in one restraint        */
static bool dr_bMixed;    /* Use sqrt of the instantaneous times        *
			   * the time averaged violation                */
static real dr_fc;	  /* Force constant for disres,                 *
			   * which is multiplied by a (possibly)        *
			   * different factor for each restraint        */
static real dr_tau;	  /* Time constant for disres			*/
static real ETerm,ETerm1; /* Exponential constants			*/

static t_drblock   drblock;

void init_disres(FILE *log,int nfa,t_inputrec *ir)
{
  dr_weighting = ir->eDisreWeighting;
  dr_fc        = ir->dr_fc;
  dr_tau       = ir->dr_tau;
  if (dr_tau == 0.0) {
    dr_bMixed = FALSE;
    ETerm = 0.0;
  } else {
    dr_bMixed = ir->bDisreMixed;
    ETerm = exp(-(ir->delta_t/dr_tau));
  }
  ETerm1 = 1.0 - ETerm;
  
  drblock.ndr = nfa/(interaction_function[F_DISRES].nratoms+1);
  snew(drblock.rav, drblock.ndr);
  snew(drblock.rt,  drblock.ndr);
  
  if (drblock.ndr > 0) {
    fprintf(log,"Done init_disre, ndr=%d\n",drblock.ndr);
    please_cite(log,"Torda89a");
  }
}

t_drblock *get_drblock(void)
{
  return &drblock;
}

real ta_disres(int nfa,t_iatom forceatoms[],t_iparams ip[],
	       rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	       matrix box,real lambda,real *dvdlambda,
	       t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
#define MAX_DRPAIRS 1000
  static real sixth=1.0/6.0;
  static real seven_three=7.0/3.0;
  
  atom_id     ai,aj;
  int         fa,fa_start,i,pair,ki,kj,m;
  int         restraint,pair_nr,pair_nr_start,rindex,npairs,type;
  rvec        dx[MAX_DRPAIRS];
  real        weight_rt_1[MAX_DRPAIRS];
  rvec        *fshift;
  real        smooth_fc,rt,Rt,Rav,rav_3,rt_1,rt_3,Rav_6,Rt_6,rt2;
  real        k0,f_scal,fmax_scal,fk_scal,fij;
  real        tav_viol,instant_viol,mixed_viol,violtot;
  real        tav_viol_Rav7,instant_viol_Rav7;
  real        up1,up2,low;
  bool        bConservative,bMixed,bViolation;
  static real exp_min_t_tau=1.0;
 
  ivec        it,jt,dt;
  
  tav_viol=instant_viol=mixed_viol=tav_viol_Rav7=instant_viol_Rav7=0;
  /* scaling factor to smoothly turn on the restraint forces *
   * when using time averaging                               */
  exp_min_t_tau *= ETerm;
  smooth_fc      = dr_fc * (1.0 - exp_min_t_tau); 

  fshift  = fr->fshift; 
  violtot = 0;
  pair_nr = 0;
  
  /* 'loop' over all atom pairs (pair_nr=fa/3) involved in restraints, *
   * the total number of atoms pairs is nfa/3                          */
  fa=0;
  while (fa < nfa) {
    restraint = forceatoms[fa];
    rindex    = ip[restraint].disres.index;

    fa_start      = fa;
    pair_nr_start = pair_nr;

    Rav_6  = 0.0;
    Rt_6   = 0.0;
  
    /* Loop over the atom pairs of 'this' restraint */
    for(pair=0; (fa<nfa) && (ip[forceatoms[fa]].disres.index == rindex); 
	pair++,pair_nr++,fa+=3) {
      if (pair >= MAX_DRPAIRS)
	fatal_error(0,"Too many distance restraint pairs");
      
      ai   = forceatoms[fa+1];
      aj   = forceatoms[fa+2];
      
      rvec_sub(x[ai],x[aj],dx[pair]);
      rt2  = iprod(dx[pair],dx[pair]);
      rt   = sqrt(rt2);
      rt_1 = invsqrt(rt2);
      rt_3 = rt_1*rt_1*rt_1;
      if (drblock.rav[pair_nr] == 0)
	rav_3 = rt_3;
      else
	rav_3 = ETerm*drblock.rav[pair_nr] + ETerm1*rt_3;

      weight_rt_1[pair] = rt_1;

      drblock.rt[pair_nr]  = rt;
      drblock.rav[pair_nr] = rav_3;
	
      Rav_6  += rav_3*rav_3;
      Rt_6   += rt_3*rt_3;
    }
    npairs=pair;
    
    /* Take action depending on restraint, calculate scalar force */
    type = ip[restraint].disres.type;
    up1  = ip[restraint].disres.up1;
    up2  = ip[restraint].disres.up2;
    low  = ip[restraint].disres.low;
    k0   = smooth_fc*ip[restraint].disres.fac;

    /* save some flops when there is only one pair */
    if (type != 2) {
      bConservative = (dr_weighting == edrwConservative) && (npairs>1);
      bMixed        = dr_bMixed;
    } else {
      bConservative = npairs>1;
      bMixed        = FALSE;
    }
    
    fmax_scal = -k0*(up2-up1); 

    Rt   = pow(Rt_6,-sixth);
    /* Compute rav from ensemble / restraint averaging */
    if (type != 2)
      Rav  = pow(Rav_6,-sixth);
    else
      /* Dirty trick to not use time averaging when type=2 */
      Rav  = Rt; 
    
    if (Rav > up1) {
      bViolation = TRUE;
      tav_viol = Rav - up1;
    } else if (Rav < low) {
      bViolation = TRUE;
      tav_viol = Rav - low;
    } else
      bViolation = FALSE;
    if (bViolation) {
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
	  mixed_viol   = sqrt(tav_viol*instant_viol);
	  f_scal       = -k0*mixed_viol;
	  violtot     += mixed_viol;
	}
      }
    }
    /* save some flops and don't use variables when they are not calculated */
    if (bViolation) {
      /* Correct the force for the number of restraints */
      if (bConservative) {
	f_scal  = max(f_scal,fmax_scal);
	if (!bMixed)
	  f_scal *= Rav/Rav_6;
	else {
	  f_scal /= 2*mixed_viol;
	  tav_viol_Rav7     = tav_viol*Rav/Rav_6;
	  instant_viol_Rav7 = instant_viol*Rt/Rt_6;
	}
      } else {
	f_scal /= (real)npairs;
	f_scal  = max(f_scal,fmax_scal);
      }    
      
      /* Exert the force ... */
      i = fa_start;
      
      for(pair=0; (pair<npairs); pair++) {
	restraint= forceatoms[i++];
	ai       = forceatoms[i++];
	aj       = forceatoms[i++];

	ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
	ki=IVEC2IS(dt);
	
	if (bConservative) {
	  rav_3 = drblock.rav[pair_nr_start+pair];
	  if (!dr_bMixed)
	    weight_rt_1[pair] *= pow(rav_3,seven_three);
	  else {
	    rt = drblock.rt[pair_nr_start+pair];
	    weight_rt_1[pair] *= tav_viol_Rav7*pow(rav_3,seven_three)+
	      instant_viol_Rav7*pow(rt,-7);
	  }
	}
	
	fk_scal  = f_scal*weight_rt_1[pair];
	
	for(m=0; (m<DIM); m++) {
	  fij            = fk_scal*dx[pair][m];
	  
	  f[ai][m]      += fij;
	  f[aj][m]      -= fij;
	  fshift[ki][m] += fij;
	  fshift[CENTRAL][m] -= fij;
	}
      }
    }
  }
  /* Return violation */
  return violtot;
}

