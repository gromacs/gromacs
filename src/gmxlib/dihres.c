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
#include "typedefs.h"
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "physics.h"
#include "vec.h"
#include "futil.h"
#include "xvgr.h"
#include "fatal.h"
#include "bondf.h"
#include "copyrite.h"
#include "disre.h"
#include "main.h"

void init_dihres(FILE *fplog,int nfa,const t_iatom forceatoms[],
		 const t_iparams ip[],
		 const t_inputrec *ir,t_fcdata *fcd)
{
  int fa;
  t_dihresdata *dd;

  
  dd = &(fcd->dihres);
  
  dd->dihre_fc        = ir->dihre_fc;
  dd->dihre_tau       = ir->dihre_tau;
  if (dd->dihre_tau == 0.0) {
    dd->ETerm = 0.0;
  } else {
    dd->ETerm = exp(-(ir->delta_t/ir->dihre_tau));
  }
  dd->ETerm1 = 1.0 - dd->ETerm;
  dd->exp_min_t_tau = 1.0;
  
  dd->nr = 0;
  for(fa=0; fa<nfa; fa+=5)
    if (fa==0 || 
	ip[forceatoms[fa-5]].dihres.label != ip[forceatoms[fa]].dihres.label)
      dd->nr++;
  dd->ndih = nfa/(interaction_function[F_DIHRES].nratoms+1);
  snew(dd->diht,dd->ndih);
  snew(dd->dihav,dd->ndih);
  
  if (dd->ndih > 0) {
    fprintf(fplog,
	    "There are %d dihedral restraints involving %d atom quartets\n",
	    dd->nr,dd->ndih);
  }
}

real ta_dihres(int nfa,const t_iatom forceatoms[],const t_iparams ip[],
	       const rvec x[],rvec f[],rvec fshift[],
	       const t_pbc *pbc,const t_graph *g,
	       real lambda,real *dvdlambda,
	       const t_mdatoms *md,t_fcdata *fcd)
{
  real vtot = 0;
  int  ai,aj,ak,al,i,k,type,typep,label,power,t1,t2,t3;
  real phi0,phi,ddphi,ddp,dp,dp2,dphi,kfac,cos_phi,sign,d2r,fc;
  rvec r_ij,r_kj,r_kl,m,n;
  
  fc  = fcd->dihres.dihre_fc;
  d2r = DEG2RAD;
  k   = 0;
  for(i=0; (i<nfa); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    al   = forceatoms[i++];
    
    phi0  = ip[type].dihres.phi*d2r;
    dphi  = ip[type].dihres.dphi*d2r;
    kfac  = ip[type].dihres.kfac*fc; 
    power = ip[type].dihres.power;
    label = ip[type].dihres.label;
    
    phi = dih_angle(x[ai],x[aj],x[ak],x[al],pbc,r_ij,r_kj,r_kl,m,n,
                    &cos_phi,&sign,&t1,&t2,&t3);	  
    /* 84 flops */
    
    if (debug)
        fprintf(debug,"dihres[%d]: %d %d %d %d : phi=%f, dphi=%f, kfac=%f, power=%d, label=%d\n",
                k++,ai,aj,ak,al,phi,dphi,kfac,power,label);

    /* phi can jump if phi0 is close to Pi/-Pi, which will cause huge
     * force changes if we just apply a normal harmonic.
     * Instead, we first calculate phi-phi0 and take it modulo (-Pi,Pi).
     * This means we will never have the periodicity problem, unless
     * the dihedral is Pi away from phiO, which is very unlikely due to
     * the potential.
     */
    dp = phi-phi0;  
    if (fabs(dp) > dphi) {
      /* dp cannot be outside (-2*pi,2*pi) */
      if (dp >= M_PI)
	dp -= 2*M_PI;
      else if(dp < -M_PI)
	dp += 2*M_PI;

      ddp = dp-dphi;      
      vtot += 0.5*kfac*ddp*ddp;
      if (phi < phi0)
	ddphi = kfac*ddp;
      else
	ddphi = -kfac*ddp;
      
      do_dih_fup(ai,aj,ak,al,(real)(-ddphi),r_ij,r_kj,r_kl,m,n,
		 f,fshift,pbc,g,x,t1,t2,t3);			/* 112		*/
    }
  }
  return vtot;
}
