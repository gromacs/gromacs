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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_mdatom_c = "$Id$";

#include "typedefs.h"
#include "mdatoms.h"
#include "smalloc.h"
#include "main.h"

t_mdatoms *atoms2md(t_atoms *atoms,bool bPert)
{
  int       i,np;
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
  snew(md->bPerturbed,md->nr);

  snew(md->cU1,md->nr);
  snew(md->cU2,md->nr);
  snew(md->cU3,md->nr);
  
  np=0;
  tm=0.0;
  for(i=0; (i<md->nr); i++) {
    md->massA[i]	= atoms->atom[i].m;
    md->massB[i]	= atoms->atom[i].mB;
    md->massT[i]	= atoms->atom[i].m;
    if (md->massA[i] != 0.0) {
      tm               += md->massT[i];
      md->invmass[i]	= 1.0/md->massT[i];
    }
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
    if (bPert) {
      md->bPerturbed[i]   = PERTURBED(atoms->atom[i]);
      if (md->bPerturbed[i])
	np++;
    }

    md->cU1[i]      	= atoms->atom[i].grpnr[egcUser1];
    md->cU2[i]      	= atoms->atom[i].grpnr[egcUser2];
    md->cU3[i]      	= atoms->atom[i].grpnr[egcUser3];
  }
  md->tmass  = tm;
  
  /*sfree(atoms->atom);
  atoms->atom=NULL;*/
  
  fprintf(stdlog,"There are %d atoms for free energy perturbation\n",np);
  
  return md;
}    

void md2atoms(t_mdatoms *md,t_atoms *atoms)
{
  int i;
  
  snew(atoms->atom,md->nr);
  for(i=0; (i<md->nr); i++) {
    atoms->atom[i].m=md->massT[i];
    atoms->atom[i].q=md->chargeT[i];
    atoms->atom[i].resnr=md->resnr[i];
    atoms->atom[i].type=md->typeA[i];
    atoms->atom[i].ptype=md->ptype[i];
    atoms->atom[i].grpnr[egcTC]=md->cTC[i];
    atoms->atom[i].grpnr[egcENER]=md->cENER[i];
    atoms->atom[i].grpnr[egcACC]=md->cACC[i];
    atoms->atom[i].grpnr[egcFREEZE]=md->cFREEZE[i];

    atoms->atom[i].grpnr[egcUser1]=md->cU1[i];
    atoms->atom[i].grpnr[egcUser2]=md->cU2[i];
    atoms->atom[i].grpnr[egcUser3]=md->cU3[i];
    atoms->atom[i].grpnr[egcXTC]  =md->cXTC[i];

  }
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

  sfree(md->cU1);
  sfree(md->cU2);
  sfree(md->cU3);
  sfree(md->cXTC);
}

void init_mdatoms(t_mdatoms *md,real lambda,bool bFirst)
{
  int  i,end;
  real L1=1.0-lambda;
  
  end=md->nr;
  for(i=0; (i<end); i++) {
    if (md->bPerturbed[i] || bFirst) {
      md->massT[i]=L1*md->massA[i]+lambda*md->massB[i];
      if (md->massT[i] != 0.0)
	md->invmass[i]=1.0/md->massT[i];
      else
	md->invmass[i]=0.0;
      md->chargeT[i]=L1*md->chargeA[i]+lambda*md->chargeB[i];
    }
  }
}




