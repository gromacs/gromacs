/*
 *       @(#) init_sh.c 1.11 2/8/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.51
 * 
 * Copyright (c) 1990-1996,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * GROwing Monsters And Cloning Shrimps
 */
#include "init_sh.h"
#include "smalloc.h"
#include "assert.h"
#include "names.h"
	
void pr_ashell(FILE *log,int nas,t_atomshell as[])
{
  int i;
  
  fprintf(log,"Atom shells\n");
  for(i=0; (i<nas); i++)
    fprintf(log,"nucl1: %5d  shell: %5d katm: %10g\n",
	    as[i].nucl1,as[i].shell,1.0/as[i].k_1);
  fprintf(log,"\n");
}

void pr_bshell(FILE *log,int nbs,t_bondshell bs[])
{
  int i;
  
  fprintf(log,"Bond shells\n");
  for(i=0; (i<nbs); i++)
    fprintf(log,"nucl1: %5d  shell: %5d  nucl2: %5d kbnd: %10g\n",
	    bs[i].nucl1,bs[i].shell,bs[i].nucl2, 1.0/bs[i].k_1);
  fprintf(log,"\n");
}

void init_shells(FILE *log,int start,int homenr,
		 t_idef *idef,t_mdatoms *md,
		 int *nashell,t_atomshell **atom_shell,
		 int *nbshell,t_bondshell **bond_shell)
{
  t_atomshell *as;
  t_bondshell *bs;
  int         n[eptNR],nas,nbs;
  int         i,j,type,ftype,nra;
  int         pt1,pt2,a1,a2;
  bool        bS1,bS2;
  t_iatom     *ia;

  for(i=0; (i<eptNR); i++)
    n[i]=0;
  for(i=start; (i<start+homenr); i++) 
    n[md->ptype[i]]++;
  for(i=0; (i<eptNR); i++)
    if (n[i]!=0)
      fprintf(log,"There are: %d %s\n",n[i],ptype_str[i]);
  
  nas=n[eptShell];
  nbs=n[eptBond];
  snew(as,nas);
  snew(bs,nbs);
  
  /* Output variables (scalars and pointers) */
  *nashell=nas;
  *nbshell=nbs;
  *atom_shell=as;
  *bond_shell=bs;
  
  /* Initiate the shell structures */    
  for(i=0; (i<nas); i++) {
    as[i].shell=-1;
    as[i].nucl1=-1;
    as[i].k_1=0;
  }
  for(i=0; (i<nbs); i++) {
    bs[i].shell=-1;
    bs[i].nucl1=-1;
    bs[i].nucl2=-1;
    bs[i].k_1=0;
  }

  /* Now fill the structures */
  nas=nbs=0;
  ia=idef->il[F_BONDS].iatoms;
  for(i=0; (i<idef->il[F_BONDS].nr); ) {
    type=ia[0];
    ftype=idef->functype[type];
    nra=interaction_function[ftype].nratoms;
    
    /* Check whether we have a bond */
    a1=ia[1];
    a2=ia[2];
    pt1=md->ptype[a1];
    pt2=md->ptype[a2];
    bS1=((pt1 == eptShell) || (pt1 == eptBond));
    bS2=((pt2 == eptShell) || (pt2 == eptBond));
    
    /* Check whether one of the particles is a shell... */
    if ((pt1 == eptShell) && !bS2) {
      assert(nas < *nashell);
      as[nas].shell = a1;
      as[nas].nucl1 = a2;
      as[nas].k_1   = 1.0/idef->iparams[type].harmonic.krA;
      nas++;
    }
    else if ((pt2 == eptShell) && !bS1) {
      assert(nas < *nashell);
      as[nas].shell = a2;
      as[nas].nucl1 = a1;
      as[nas].k_1   = 1.0/idef->iparams[type].harmonic.krA;
      nas++;
    }
    else if ((pt1 == eptBond) && !bS2) {
      for(j=0; (j<nbs); j++)
	if (bs[j].shell == a1) {
	  bs[j].nucl2 = a2;
	  break;
	}
      if (j == nbs) {
	assert(nbs < *nbshell);
	bs[nbs].shell = a1;
	bs[nbs].nucl1 = a2;
	bs[nbs].k_1   = 1.0/idef->iparams[type].harmonic.krA;
	nbs++;
      }
    }
    else if ((pt2 == eptBond) && !bS1) {
      for(j=0; (j<nbs); j++)
	if (bs[j].shell == a2) {
	  bs[j].nucl2 = a1;
	  break;
	}
      if (j == nbs) {
	assert(nbs < *nbshell);
	bs[nbs].shell = a2;
	bs[nbs].nucl1 = a1;
	bs[nbs].k_1   = 1.0/idef->iparams[type].harmonic.krA;
	nbs++;
      }
    }
    ia += nra+1;
    i  += nra+1;
  }
  
  /* Verify whether it's all correct */
  for(i=0; (i<nas); i++) {
    if ((as[i].shell==-1) || (as[i].nucl1==-1) || (as[i].k_1==0))
      fatal_error(0,"atom_shell[%d] sucks: shell=%d, nucl1=%d, k_1=%g",
		  i,as[i].shell,as[i].nucl1,as[i].k_1);
  }
  for(i=0; (i<nbs); i++) {
    if ((bs[i].shell==-1) || (bs[i].nucl1==-1) ||
	(bs[i].nucl2==-1) || (bs[i].k_1==0))
      fatal_error(0,"bond_shell[%d] sucks: shell=%d, nucl1=%d,"
		  " nucl2=%d, k_1=%g",
		  i,bs[i].shell,bs[i].nucl1,bs[i].nucl2,bs[i].k_1);
  }
#ifdef DEBUG
  pr_ashell(log,nas,as);
  pr_bshell(log,nbs,bs);
  fflush(log);
#endif
}

