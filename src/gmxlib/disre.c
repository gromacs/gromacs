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
 * Good gRace! Old Maple Actually Chews Slate
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

void init_disres(FILE *log,int nrestraints,t_inputrec *ir)
{
  dr_weighting = ir->eDisreWeighting;
  dr_fc  = ir->dr_fc;
  dr_tau = ir->dr_tau;
  if (dr_tau == 0.0)
    ETerm=0.0;
  else
    ETerm=exp(-(ir->delta_t/dr_tau));

  ETerm1=1.0-ETerm;
  
  drblock.ndr    = nrestraints/(interaction_function[F_DISRES].nratoms+1);
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

real ta_disres(FILE *log,int nrestraints,t_iatom forceatoms[],t_iparams ip[],
	       rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	       matrix box,real lambda,real *dvdlambda,
	       t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
#define MAX_DRPAIRS 100
  static real sixth=1.0/6.0;
  
  atom_id     ai,aj;
  int         i,k,ki,kj,m,restraint,index,rindex,i0,npairs;
  rvec        dx[MAX_DRPAIRS];
  real        weight_rt_1[MAX_DRPAIRS];
  rvec        *fshift;
  real        smooth_fc,rt,rav,rav_3,rt_1,rt_3,rav_6,rt_6,rt2;
  real        k0,f_scal,fmax_scal,fk_scal,fij,viol,pseudo_viol,violtot;
  real        up1,up2;
  bool        bConservative;
  static real exp_min_t_tau=1.0;

  /* scaling factor to smoothly turn on the restraint forces *
   * when using time averaging                               */
  exp_min_t_tau *= ETerm;
  smooth_fc      = dr_fc * (1.0 - exp_min_t_tau); 

  bConservative = (dr_weighting == edrwConservative);
  fshift  = fr->fshift; 
  violtot = 0;
  index   = 0;
  
  for(i=0; (i<nrestraints); ) {  
    i0        = i;
    restraint = forceatoms[i];
    rindex    = ip[restraint].disres.index;

    rav_6  = 0.0;
    rt_6   = 0.0;
    k      = 0;
  
    /* Loop over the distance restraints that should be treated 
     * simultaneously
     */
    for(k=0; (ip[forceatoms[i]].disres.index == rindex) && (i<nrestraints); 
	k++,index++,i+=3) {
      if (k >= MAX_DRPAIRS)
	fatal_error(0,"Too many distance restraint pairs");
      
      ai   = forceatoms[i+1];
      aj   = forceatoms[i+2];
      
      rvec_sub(x[ai],x[aj],dx[k]);
      rt2  = iprod(dx[k],dx[k]);
      rt   = sqrt(rt2);
      rt_1 = invsqrt(rt2);
      rt_3 = rt_1*rt_1*rt_1;
      if (bConservative)
	weight_rt_1[k] = rt_3*rt_3*rt_1*rt_1;
      else
	weight_rt_1[k] = rt_1;
      if (drblock.rav[index] == 0)
	rav_3 = rt_3;
      else
	rav_3 = ETerm*drblock.rav[index] + ETerm1*rt_3;
      
      drblock.rt[index]  = rt;
      drblock.rav[index] = rav_3;
	
      rav_6  += rav_3*rav_3;
      rt_6   += rt_3*rt_3;
    }
    npairs=k;
    
    /* Take action depending on restraint, calculate scalar force */
    up1       = ip[restraint].disres.up1;
    up2       = ip[restraint].disres.up2;
    k0        = smooth_fc*ip[restraint].disres.fac;

    fmax_scal = -k0*(up2-up1); 

    /* Compute rav from ensemble / restraint averaging */
    rav  = pow(rav_6,-sixth);
    rt   = pow(rt_6,-sixth);

    if (rav > up1)
      viol = rav-up1;
    else
      viol = 0.0;

    if (!dr_bMixed) {
      f_scal = -k0*viol;
    } else
      if (rt > up1) {
	pseudo_viol = sqrt(viol*(rt - up1));
	f_scal      = -k0*pseudo_viol;
      } else
	f_scal = 0.0;
    violtot += viol;
    
    /* Correct the force for the number of restraints */
    if (bConservative) {
      f_scal  = max(f_scal,fmax_scal);
      f_scal *= pow(rav,7);
    } else {
      f_scal /= (real)npairs;
      f_scal  = max(f_scal,fmax_scal);
    }    

    /* Exert the force ... */
    i    = i0;

    for(k=0; (k<npairs); k++) {
      restraint= forceatoms[i++];
      ai       = forceatoms[i++];
      aj       = forceatoms[i++];
      ki       = SHIFT_INDEX(g,ai);
      kj       = SHIFT_INDEX(g,aj);
      
      fk_scal  = f_scal*weight_rt_1[k];
      
      for(m=0; (m<DIM); m++) {
	fij            = fk_scal*dx[k][m];
	
	f[ai][m]      += fij;
	f[aj][m]      -= fij;
	fshift[ki][m] += fij;
	fshift[kj][m] -= fij;
      }
    }
  }
  /* Return violation */
  return violtot;
}

