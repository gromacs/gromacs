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
#include "mtop_util.h"
#include "gmx_omp_nthreads.h"

#define ALMOST_ZERO 1e-30

t_mdatoms *init_mdatoms(FILE *fp,gmx_mtop_t *mtop,gmx_bool bFreeEnergy)
{
  int    mb,a,g,nmol;
  double tmA,tmB;
  t_atom *atom;
  t_mdatoms *md;
  gmx_mtop_atomloop_all_t aloop;
  t_ilist *ilist;

  snew(md,1);

  md->nenergrp = mtop->groups.grps[egcENER].nr;
  md->bVCMgrps = FALSE;
  tmA = 0.0;
  tmB = 0.0;

  aloop = gmx_mtop_atomloop_all_init(mtop);
  while(gmx_mtop_atomloop_all_next(aloop,&a,&atom)) {
    if (ggrpnr(&mtop->groups,egcVCM,a) > 0)
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

  md->bOrires = gmx_mtop_ftype_count(mtop,F_ORIRES);

  return md;
}

void atoms2md(gmx_mtop_t *mtop,t_inputrec *ir,
	      int nindex,int *index,
	      int start,int homenr,
	      t_mdatoms *md)
{
  gmx_mtop_atomlookup_t alook;
  int       i;
  t_grpopts *opts;
  gmx_groups_t *groups;
  gmx_molblock_t *molblock;

  opts = &ir->opts;

  groups = &mtop->groups;

  molblock = mtop->molblock;

  /* Index==NULL indicates particle decomposition,
   * unless we have an empty DD node, so also check for homenr and start.
   * This should be signaled properly with an extra parameter or nindex==-1.
   */
  if (index == NULL && (homenr > 0 || start > 0)) {
    md->nr = mtop->natoms;
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
    if (md->bOrires)
      srenew(md->cORF,md->nalloc);
    if (md->nPerturbed)
      srenew(md->bPerturbed,md->nalloc);
    
    /* Note that these user t_mdatoms array pointers are NULL
     * when there is only one group present.
     * Therefore, when adding code, the user should use something like:
     * gprnrU1 = (md->cU1==NULL ? 0 : md->cU1[localatindex])
     */
    if (mtop->groups.grpnr[egcUser1] != NULL)
      srenew(md->cU1,md->nalloc);
    if (mtop->groups.grpnr[egcUser2] != NULL)
      srenew(md->cU2,md->nalloc);
    
    if (ir->bQMMM)
      srenew(md->bQM,md->nalloc);
    if (ir->bAdress) {
      srenew(md->wf,md->nalloc);
      srenew(md->tf_table_index,md->nalloc);
    }
  }

  alook = gmx_mtop_atomlookup_init(mtop);

#pragma omp parallel for num_threads(gmx_omp_nthreads_get(emntDefault)) schedule(static)
  for(i=0; i<md->nr; i++) {
    int     g,ag,molb;
    real    mA,mB,fac;
    t_atom  *atom;

    if (index == NULL) {
      ag = i;
    } else {
      ag   = index[i];
    }
    gmx_mtop_atomnr_to_atom(alook,ag,&atom);

    if (md->cFREEZE) {
      md->cFREEZE[i] = ggrpnr(groups,egcFREEZE,ag);
    }
        if (EI_ENERGY_MINIMIZATION(ir->eI))
        {
            /* Displacement is proportional to F, masses used for constraints */
            mA = 1.0;
            mB = 1.0;
        }
        else if (ir->eI == eiBD)
        {
            /* With BD the physical masses are irrelevant.
             * To keep the code simple we use most of the normal MD code path
             * for BD. Thus for constraining the masses should be proportional
             * to the friction coefficient. We set the absolute value such that
             * m/2<(dx/dt)^2> = m/2*2kT/fric*dt = kT/2 => m=fric*dt/2
             * Then if we set the (meaningless) velocity to v=dx/dt, we get the
             * correct kinetic energy and temperature using the usual code path.
             * Thus with BD v*dt will give the displacement and the reported
             * temperature can signal bad integration (too large time step).
             */
            if (ir->bd_fric > 0)
            {
                mA = 0.5*ir->bd_fric*ir->delta_t;
                mB = 0.5*ir->bd_fric*ir->delta_t;
            }
            else
            {
                /* The friction coefficient is mass/tau_t */
                fac = ir->delta_t/opts->tau_t[md->cTC ? groups->grpnr[egcTC][ag] : 0];
                mA = 0.5*atom->m*fac;
                mB = 0.5*atom->mB*fac;
            }
        }
        else
        {
            mA = atom->m;
            mB = atom->mB;
        }
    if (md->nMassPerturbed) {
      md->massA[i]	= mA;
      md->massB[i]	= mB;
    }
    md->massT[i]	= mA;
    if (mA == 0.0) {
      md->invmass[i]    = 0;
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
      md->cTC[i]	= groups->grpnr[egcTC][ag];
    md->cENER[i]	=
      (groups->grpnr[egcENER] ? groups->grpnr[egcENER][ag] : 0);
    if (md->cACC)
      md->cACC[i]	= groups->grpnr[egcACC][ag];
    if (md->cVCM)
      md->cVCM[i]     	= groups->grpnr[egcVCM][ag];
    if (md->cORF)
      md->cORF[i]     	= groups->grpnr[egcORFIT][ag];

    if (md->cU1)
      md->cU1[i]     	= groups->grpnr[egcUser1][ag];
    if (md->cU2)
      md->cU2[i]     	= groups->grpnr[egcUser2][ag];

    if (ir->bQMMM) {
      if (groups->grpnr[egcQMMM] == 0 || 
	  groups->grpnr[egcQMMM][ag] < groups->grps[egcQMMM].nr-1) {
	md->bQM[i]      = TRUE;
      } else {
	md->bQM[i]      = FALSE;
      }
    }
    /* Initialize AdResS weighting functions to adressw */
    if (ir->bAdress){
       md->wf[i]           = 1.0;
        /* if no tf table groups specified, use default table */
       md->tf_table_index[i] = DEFAULT_TF_TABLE;
       if (ir->adress->n_tf_grps > 0){
            /* if tf table groups specified, tf is only applied to thoose energy groups*/
            md->tf_table_index[i] = NO_TF_TABLE;
            /* check wether atom is in one of the relevant energy groups and assign a table index */
            for (g=0; g<ir->adress->n_tf_grps; g++){
                if (md->cENER[i] == ir->adress->tf_table_index[g]){
                   md->tf_table_index[i] = g;
                }
            }
        }
    }
  }

  gmx_mtop_atomlookup_destroy(alook);

  md->start  = start;
  md->homenr = homenr;
  md->lambda = 0;
}

void update_mdatoms(t_mdatoms *md,real lambda)
{
  int    al,end;
  real   L1=1.0-lambda;
  
  end=md->nr;

  if (md->nMassPerturbed) {
    for(al=0; (al<end); al++) {
      if (md->bPerturbed[al]) {
	md->massT[al] = L1*md->massA[al]+ lambda*md->massB[al];
	if (md->invmass[al] > 1.1*ALMOST_ZERO)
	  md->invmass[al] = 1.0/md->massT[al];
      }
    }
    md->tmass = L1*md->tmassA + lambda*md->tmassB;
  } else {
    md->tmass = md->tmassA;
  }
  md->lambda = lambda;
}
