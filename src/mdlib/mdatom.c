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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "mdatoms.h"
#include "smalloc.h"
#include "main.h"
#include "qmmm.h"

#define ALMOST_ZERO 1e-30

t_mdatoms *init_mdatoms(FILE *fp,t_atoms *atoms,bool bFreeEnergy)
{
  int    i,g;
  double tmA,tmB;
  t_atom *atom;
  t_mdatoms *md;
  
  snew(md,1);

  md->bVCMgrps = FALSE;
  tmA = 0.0;
  tmB = 0.0;
  for(i=0; i<atoms->nr; i++) {
    atom = &atoms->atom[i];

    if (atom->grpnr[egcVCM] > 0)
      md->bVCMgrps = TRUE;
    
    if (bFreeEnergy && PERTURBED(*atom)) {
      md->nPerturbed++;
      if (atom->mB != atom->m)
	md->nMassPerturbed++;
      if (atom->qB != atom->q)
	md->nChargePerturbed++;
    }
    
    tmA += atom->m;
    tmB += atom->mB;
  }

  md->tmassA = tmA;
  md->tmassB = tmB;
  
  if (bFreeEnergy && fp)
    fprintf(fp,
	    "There are %d atoms and %d charges for free energy perturbation\n",
	    md->nPerturbed,md->nChargePerturbed);

  return md;
}

void atoms2md(t_atoms *atoms,t_inputrec *ir,int norires,
	      int nindex,int *index,
	      int start,int homenr,
	      t_mdatoms *md)
{
  int       i,g;
  real      mA,mB,fac;
  t_atom    *atom;
  t_grpopts *opts;

  opts = &ir->opts;

  if (index == NULL) {
    md->nr = atoms->nr;
  } else {
    md->nr = nindex;
  }

  if (md->nr > md->nalloc) {
    md->nalloc = over_alloc_dd(md->nr);

    if (md->nMassPerturbed) {
      srenew(md->massA,md->nalloc);
      srenew(md->massB,md->nalloc);
    }
    srenew(md->massT,md->nalloc);
    srenew(md->invmass,md->nalloc);
    srenew(md->chargeA,md->nalloc);
    if (md->nPerturbed) {
      srenew(md->chargeB,md->nalloc);
    }
    srenew(md->typeA,md->nalloc);
    if (md->nPerturbed) {
      srenew(md->typeB,md->nalloc);
    }
    srenew(md->ptype,md->nalloc);
    if (opts->ngtc > 1) {
      srenew(md->cTC,md->nalloc);
      /* We always copy cTC with domain decomposition */
    }
    srenew(md->cENER,md->nalloc);
    if (opts->ngacc > 1)
      srenew(md->cACC,md->nalloc);
    if (opts->nFreeze &&
	(opts->ngfrz > 1 ||
	 opts->nFreeze[0][XX] || opts->nFreeze[0][YY] || opts->nFreeze[0][ZZ]))
      srenew(md->cFREEZE,md->nalloc);
    if (md->bVCMgrps)
      srenew(md->cVCM,md->nalloc);
    if (norires)
      srenew(md->cORF,md->nalloc);
    if (md->nPerturbed)
      srenew(md->bPerturbed,md->nalloc);
    
    /* The user should fix this */
    if (FALSE)
      srenew(md->cU1,md->nalloc);
    if (FALSE)
      srenew(md->cU2,md->nalloc);
    
    if (ir->bQMMM)
      srenew(md->bQM,md->nalloc);
  }

  for(i=0; (i<md->nr); i++) {
    if (index == NULL)
      atom = &atoms->atom[i];
    else
      atom = &atoms->atom[index[i]];

    if (md->cFREEZE)
      md->cFREEZE[i]	= atom->grpnr[egcFREEZE];
    if (EI_ENERGY_MINIMIZATION(ir->eI)) {
      mA = 1.0;
      mB = 1.0;
    } else if (ir->eI == eiBD) {
      /* Make the mass proportional to the friction coefficient for BD.
       * This is necessary for the constraint algorithms.
       */
      if (ir->bd_fric) {
	mA = ir->bd_fric*ir->delta_t;
	mB = ir->bd_fric*ir->delta_t;
      } else {
	fac = ir->delta_t/opts->tau_t[atom->grpnr[egcTC]];
	mA = atom->m*fac;
	mB = atom->mB*fac;
      }
    } else {
      mA = atom->m;
      mB = atom->mB;
    }
    if (md->nMassPerturbed) {
      md->massA[i]	= mA;
      md->massB[i]	= mB;
    }
    md->massT[i]	= mA;
    if (mA == 0.0) {
      md->invmass[i] = 0;
    } else if (md->cFREEZE) {
      g = md->cFREEZE[i];
      if (opts->nFreeze[g][XX] && opts->nFreeze[g][YY] && opts->nFreeze[g][ZZ])
	/* Set the mass of completely frozen particles to ALMOST_ZERO iso 0
	 * to avoid div by zero in lincs or shake.
	 * Note that constraints can still move a partially frozen particle.
	 */
	md->invmass[i]	= ALMOST_ZERO;
      else
	md->invmass[i]	= 1.0/mA;
    } else {
      md->invmass[i]	= 1.0/mA;
    }
    md->chargeA[i]	= atom->q;
    md->typeA[i]	= atom->type;
    if (md->nPerturbed) {
      md->chargeB[i]	= atom->qB;
      md->typeB[i]	= atom->typeB;
      md->bPerturbed[i] = PERTURBED(*atom);
    }
    md->ptype[i]	= atom->ptype;
    if (md->cTC)
      md->cTC[i]	= atom->grpnr[egcTC];
    md->cENER[i]	= atom->grpnr[egcENER];
    if (md->cACC)
      md->cACC[i]	= atom->grpnr[egcACC];
    if (md->cVCM)
      md->cVCM[i]      	= atom->grpnr[egcVCM];
    if (md->cORF)
      md->cORF[i]      	= atom->grpnr[egcORFIT];

    if (md->cU1)
      md->cU1[i]      	= atom->grpnr[egcUser1];
    if (md->cU2)
      md->cU2[i]      	= atom->grpnr[egcUser2];

    if (ir->bQMMM) {
      if (atom->grpnr[egcQMMM] < atoms->grps[egcQMMM].nr-1) {
	md->bQM[i]      = TRUE;
      } else {
	md->bQM[i]      = FALSE;
      }
    }
  }

  md->start  = start;
  md->homenr = homenr;
  md->lambda = 0;
}

void update_mdatoms(t_mdatoms *md,real lambda)
{
  int    i,end;
  real   L1=1.0-lambda;
  
  end=md->nr;

  if (md->nMassPerturbed) {
    for(i=0; (i<end); i++) {
      if (md->bPerturbed[i]) {
	md->massT[i] = L1*md->massA[i] + lambda*md->massB[i];
	if (md->invmass[i] > 1.1*ALMOST_ZERO)
	  md->invmass[i] = 1.0/md->massT[i];
      }
    }
    md->tmass = L1*md->tmassA + lambda*md->tmassB;
  } else {
    md->tmass = md->tmassA;
  }
  md->lambda = lambda;
}
