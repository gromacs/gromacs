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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_mdatom_c = "$Id$";

#include "typedefs.h"
#include "mdatoms.h"
#include "smalloc.h"
#include "main.h"

#define ALMOST_ZERO 1e-30

t_mdatoms *atoms2md(FILE *fp,t_atoms *atoms,ivec nFreeze[],
		    bool bBD,real delta_t,real fric,real tau_t[],
		    bool bPert,bool bFree)
{
  int       i,np,g;
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
  snew(md->chargeT,md->nr);
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
  snew(md->bPerturbed,md->nr);

  snew(md->cU1,md->nr);
  snew(md->cU2,md->nr);
  
  np=0;
  tm=0.0;
  for(i=0; (i<md->nr); i++) {
    if (bBD) {
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
    if (md->massA[i] != 0.0) {
      tm               += md->massT[i];
      g = md->cFREEZE[i];
      if (nFreeze[g][XX] && nFreeze[g][YY] && nFreeze[g][ZZ])
	/* Set the mass of completely frozen particles to ALMOST_ZERO iso 0
	   to avoid div by zero in lincs or shake.
	   Note that constraints can still move a partially frozen particle. */
	md->invmass[i]	= ALMOST_ZERO;
      else
	md->invmass[i]	= 1.0/md->massT[i];
    }
    if (bPert) {
      md->bPerturbed[i]   = PERTURBED(atoms->atom[i]);
      if (md->bPerturbed[i])
	np++;
    }

    md->cU1[i]      	= atoms->atom[i].grpnr[egcUser1];
    md->cU2[i]      	= atoms->atom[i].grpnr[egcUser2];
  }
  md->tmass  = tm;

  if (bFree) {  
    sfree(atoms->atom);
    atoms->atom=NULL;
  }
  
  if (fp)
    fprintf(fp,"There are %d atoms for free energy perturbation\n",np);
  
  return md;
}    

void md2atoms(t_mdatoms *md,t_atoms *atoms,bool bFree)
{
  int i;
  
  snew(atoms->atom,md->nr);
  for(i=0; (i<md->nr); i++) {
    atoms->atom[i].m                = md->massT[i];
    atoms->atom[i].q                = md->chargeT[i];
    atoms->atom[i].resnr            = md->resnr[i];
    atoms->atom[i].type             = md->typeA[i];
    atoms->atom[i].ptype            = md->ptype[i];
    atoms->atom[i].grpnr[egcTC]     = md->cTC[i];
    atoms->atom[i].grpnr[egcENER]   = md->cENER[i];
    atoms->atom[i].grpnr[egcACC]    = md->cACC[i];
    atoms->atom[i].grpnr[egcFREEZE] = md->cFREEZE[i];
    atoms->atom[i].grpnr[egcVCM]    = md->cVCM[i];
    atoms->atom[i].grpnr[egcXTC]    = md->cXTC[i];

    atoms->atom[i].grpnr[egcUser1]  = md->cU1[i];
    atoms->atom[i].grpnr[egcUser2]  = md->cU2[i];

  }
  if (bFree) {
    sfree(md->massA);
    sfree(md->massB);
    sfree(md->massT);
    sfree(md->invmass);
    sfree(md->chargeA);
    sfree(md->chargeB);
    sfree(md->chargeT);
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
    
    sfree(md->cU1);
    sfree(md->cU2);
  }
}

void init_mdatoms(t_mdatoms *md,real lambda,bool bFirst)
{
  static real lambda0;
  int    i,end;
  real   L1=1.0-lambda;
  
  if (bFirst)
    lambda0 = lambda;
  end=md->nr;
  
  /* Only do this loop the first time, or when lambda has changed.
   * One could also check whether there is any perturbed atom at all,
   * but if you don't have perturbed atoms, it does not make sense to modify lambda.
   * In principle this has to be parallellized, although it would mean extra 
   * communication. Basically only the charges are used on other nodes...
   */
  if (bFirst || (lambda0 != lambda)) {
    for(i=0; (i<end); i++) {
      if (md->bPerturbed[i] || bFirst) {
	md->massT[i]=L1*md->massA[i]+lambda*md->massB[i];
	if (md->invmass[i] > 1.1*ALMOST_ZERO)
	  md->invmass[i]=1.0/md->massT[i];
	md->chargeT[i]=L1*md->chargeA[i]+lambda*md->chargeB[i];
      }
    }
  }
  lambda0 = lambda;
}




