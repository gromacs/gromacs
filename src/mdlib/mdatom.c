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

t_mdatoms *atoms2md(FILE *fp,t_atoms *atoms,ivec nFreeze[],
		    int eI,real delta_t,real fric,real tau_t[],
		    bool bPert,bool bFree)
{
  int       i,g;
  real      fac;
  double    tm;
  t_mdatoms *md;
  
  snew(md,1);
  md->nr = atoms->nr;
  snew(md->massA,md->nr);
  snew(md->massB,md->nr);
  snew(md->massT,md->nr);
  snew(md->invmass,md->nr);
  snew(md->chargeA,md->nr);
  snew(md->chargeB,md->nr);
  snew(md->resnr,md->nr);
  snew(md->typeA,md->nr);
  snew(md->typeB,md->nr);
  snew(md->ptype,md->nr);
  snew(md->cTC,md->nr);
  snew(md->cENER,md->nr);
  snew(md->cACC,md->nr);
  snew(md->cFREEZE,md->nr);
  snew(md->cXTC,md->nr);
  snew(md->cVCM,md->nr);
  snew(md->cORF,md->nr);
  snew(md->bPerturbed,md->nr);

  snew(md->cU1,md->nr);
  snew(md->cU2,md->nr);
  
  /* QMMM additions */
  snew(md->cQMMM,md->nr);
  snew(md->nucnum,md->nr);
  /* This routine currently only works for ffG43a2 */
  /*  atomic_number(md->nr,atoms->atomtype,md->nucnum); */
  snew(md->bQM,md->nr);

  md->nPerturbed=0;
  md->bMassPerturbed=FALSE;
  md->bChargePerturbed=FALSE;
  tm=0.0;
  for(i=0; (i<md->nr); i++) {
    if (EI_ENERGY_MINIMIZATION(eI)) {
      md->massA[i]	= 1.0;
      md->massB[i]	= 1.0;
    } else if (eI == eiBD) {
      /* Make the mass proportional to the friction coefficient for BD.
       * This is necessary for the constraint algorithms.
       */
      if (fric) {
	md->massA[i]	= fric*delta_t;
	md->massB[i]	= fric*delta_t;
      } else {
	fac = delta_t/tau_t[atoms->atom[i].grpnr[egcTC]];
	md->massA[i]	= atoms->atom[i].m*fac;
	md->massB[i]	= atoms->atom[i].mB*fac;
      }
    } else {
      md->massA[i]	= atoms->atom[i].m;
      md->massB[i]	= atoms->atom[i].mB;
    }
    md->massT[i]	= md->massA[i];
    md->chargeA[i]	= atoms->atom[i].q;
    md->chargeB[i]	= atoms->atom[i].qB;
    md->resnr[i]	= atoms->atom[i].resnr;
    md->typeA[i]	= atoms->atom[i].type;
    md->typeB[i]	= atoms->atom[i].typeB;
    md->ptype[i]	= atoms->atom[i].ptype;
    md->cTC[i]		= atoms->atom[i].grpnr[egcTC];
    md->cENER[i]	= atoms->atom[i].grpnr[egcENER];
    md->cACC[i]		= atoms->atom[i].grpnr[egcACC];
    md->cFREEZE[i]	= atoms->atom[i].grpnr[egcFREEZE];
    md->cXTC[i]      	= atoms->atom[i].grpnr[egcXTC];
    md->cVCM[i]      	= atoms->atom[i].grpnr[egcVCM];
    md->cORF[i]      	= atoms->atom[i].grpnr[egcORFIT];
    if (md->massA[i] != 0.0) {
      tm               += md->massT[i];
      g = md->cFREEZE[i];
      if (nFreeze[g][XX] && nFreeze[g][YY] && nFreeze[g][ZZ])
	/* Set the mass of completely frozen particles to ALMOST_ZERO iso 0
	   to avoid div by zero in lincs or shake.
	   Note that constraints can still move a partially frozen particle. */
	md->invmass[i]	= ALMOST_ZERO;
      else if (md->massT[i] == 0)
	md->invmass[i]  = 0;
      else
	md->invmass[i]	= 1.0/md->massT[i];
    }
    if (bPert) {
      md->bPerturbed[i] = PERTURBED(atoms->atom[i]);
      if (md->bPerturbed[i]) {
	md->nPerturbed++;
	if (atoms->atom[i].mB != atoms->atom[i].m)
	  md->bMassPerturbed = TRUE;
	if (atoms->atom[i].qB != atoms->atom[i].q)
	  md->bChargePerturbed = TRUE;
      }
    }

    md->cU1[i]      	= atoms->atom[i].grpnr[egcUser1];
    md->cU2[i]      	= atoms->atom[i].grpnr[egcUser2];
    md->cQMMM[i]        = atoms->atom[i].grpnr[egcQMMM];
    if ((md->cQMMM[i])<(atoms->grps[egcQMMM].nr-1)){
      md->bQM[i]          = TRUE;
    } else {
      md->bQM[i] = FALSE;
    }
  }
  md->tmass  = tm;

  if (bFree) {  
    sfree(atoms->atom);
    atoms->atom=NULL;
    atoms->nr=0;
    atoms->nres=0;
    atoms->ngrpname=0;
  }
  
  if (bPert && fp)
    fprintf(fp,"There are %d atoms for free energy perturbation\n",
	    md->nPerturbed);
  
  return md;
}    

void md2atoms(t_mdatoms *md,t_atoms *atoms,bool bFree)
{
  int i;
  
  snew(atoms->atom,md->nr);
  for(i=0; (i<md->nr); i++) {
    atoms->atom[i].m                = md->massT[i];
    atoms->atom[i].resnr            = md->resnr[i];
    atoms->atom[i].type             = md->typeA[i];
    atoms->atom[i].ptype            = md->ptype[i];
    atoms->atom[i].grpnr[egcTC]     = md->cTC[i];
    atoms->atom[i].grpnr[egcENER]   = md->cENER[i];
    atoms->atom[i].grpnr[egcACC]    = md->cACC[i];
    atoms->atom[i].grpnr[egcFREEZE] = md->cFREEZE[i];
    atoms->atom[i].grpnr[egcVCM]    = md->cVCM[i];
    atoms->atom[i].grpnr[egcXTC]    = md->cXTC[i];
    atoms->atom[i].grpnr[egcORFIT]  = md->cORF[i];

    atoms->atom[i].grpnr[egcUser1]  = md->cU1[i];
    atoms->atom[i].grpnr[egcUser2]  = md->cU2[i];
    atoms->atom[i].grpnr[egcQMMM]   = md->cQMMM[i];

  }
  if (bFree) {
    sfree(md->massA);
    sfree(md->massB);
    sfree(md->massT);
    sfree(md->invmass);
    sfree(md->chargeA);
    sfree(md->chargeB);
    sfree(md->resnr);
    sfree(md->typeA);
    sfree(md->typeB);
    sfree(md->ptype);
    sfree(md->cTC);
    sfree(md->cENER);
    sfree(md->cACC);
    sfree(md->cFREEZE);
    sfree(md->cVCM);
    sfree(md->cXTC);
    sfree(md->cORF);
    
    sfree(md->cU1);
    sfree(md->cU2);
    sfree(md->cQMMM);
  }
}

void update_mdatoms(t_mdatoms *md,real lambda, bool bFirst)
{
  int    i,end;
  real   L1=1.0-lambda;
  
  end=md->nr;

  for(i=0; (i<end); i++) {
    if (md->bPerturbed[i] || bFirst) {
      md->massT[i]=L1*md->massA[i]+lambda*md->massB[i];
      if (md->invmass[i] > 1.1*ALMOST_ZERO)
	md->invmass[i]=1.0/md->massT[i];
    }
  }
}


