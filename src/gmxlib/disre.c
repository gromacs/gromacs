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

/* Some local variables */
static real dr_fc;	  /* Force constant for disres			*/
static real dr_tau;	  /* Time constant for disres			*/
static real ETerm,ETerm1; /* Exponential constants			*/

static t_drblock   drblock;

void init_disres(FILE *log,int nbonds,t_inputrec *ir)
{
  int i,type,ftype,ftind;

  dr_fc  = ir->dr_fc;
  dr_tau = ir->dr_tau;
  if (dr_tau == 0.0)
    ETerm=0.0;
  else
    ETerm=exp(-(ir->delta_t/dr_tau));

  ETerm1=1.0-ETerm;
  
  drblock.ndr    = nbonds/(interaction_function[F_DISRES].nratoms+1);
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

real ta_disres(FILE *log,int nbonds,t_iatom forceatoms[],t_iparams ip[],
	       rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	       matrix box,real lambda,real *dvdlambda,
	       t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
#define MAX_DRPAIRS 100
  static real third=1.0/3.0;
  static real sixth=1.0/6.0;
  
  atom_id     ai,aj;
  int         i,j,k,kmax,ki,kj,m,type,index,rindex,i0;
  rvec        dx[MAX_DRPAIRS];
  real        rt_1[MAX_DRPAIRS];
  rvec        *fshift;
  real        rt,rN,rav,rav_3,rt_11,rt_3,rav_6,rt_6,rt222,k_1;
  real        k0,f_scal,fmax_scal,fk_scal,fij,viol,viol2,violtot;
  real        rx0,rx1;
  
  fshift  = fr->fshift; 
  violtot = 0;
  index   = 0;
  
  for(i=0; (i<nbonds); ) {  
    i0     = i;
    type   = forceatoms[i];
    rindex = ip[type].disres.index;
    
    rav_6  = 0.0;
    rt_6   = 0.0;
    k      = 0;
  
    /* Loop over the distance restraints that should be treated 
     * simultaneously
     */
    for(k=0; (ip[forceatoms[i]].disres.index == rindex) && (i<nbonds); 
	k++,index++,i+=3) {
      if (k >= MAX_DRPAIRS)
	fatal_error(0,"Too many distance restraint pairs");
      
      ai   = forceatoms[i+1];
      aj   = forceatoms[i+2];
      
      rvec_sub(x[ai],x[aj],dx[k]);
      rt222   = iprod(dx[k],dx[k]);
      rt      = sqrt(rt222);
      rt_11   = invsqrt(rt222);
      rt_3    = rt_11*rt_11*rt_11;
      rt_1[k] = rt_11;
      if (drblock.rav[index] == 0)
	rav_3 = rt_3;
      else
	rav_3 = ETerm*drblock.rav[index] + ETerm1*rt_3;
      
      drblock.rt[index]  = rt;
      drblock.rav[index] = rav_3;
	
      rav_6  += rav_3*rav_3;
      rt_6   += rt_3*rt_3;
    }
    
    /* Take action depending on type, calculate scalar force */
    rx0       = ip[type].disres.rx0;
    rx1       = ip[type].disres.rx1;
    
    k0        = dr_fc;
    fmax_scal = -k0*(rx1-rx0); 
    f_scal    = 0.0;
    viol      = 0.0;
    
    /* Compute rav from ensemble / restraint averaging */
    rav  = pow(rav_6,-sixth);
    rt   = pow(rt_6,-sixth);
    
    if (rav > rx0) {
      if (ip[type].disres.type == 1) {
	if (rt > rx0) {
	  viol   = sqrt((rav - rx0)*(rt - rx0));
	  f_scal = -k0*viol;
	}
      }
      else if (ip[type].disres.type == 2) {
	viol   = rav-rx0;
	f_scal = -k0*viol;
      }
      else
	fatal_error(0,"No such disres type: %d\n",
		    ip[type].disres.type);
    }
    violtot += viol;
    
    /* Divide the force by number of restraints */
    k_1     = 1.0/(real) k;
    f_scal *= k_1;
    f_scal  = max(f_scal,fmax_scal);
    
    /* Exert the force ... */
    i    = i0;
    kmax = k;
    k    = 0;
    
    for(k=0; (k<kmax); k++) {
      type     = forceatoms[i++];
      ai       = forceatoms[i++];
      aj       = forceatoms[i++];
      ki       = SHIFT_INDEX(g,ai);
      kj       = SHIFT_INDEX(g,aj);
      
      fk_scal  = f_scal*rt_1[k];
      
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

