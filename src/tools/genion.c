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
static char *SRCID_genion_c = "$Id$";

#include "copyrite.h"
#include "string2.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "confio.h"
#include "assert.h"
#include "statutil.h"
#include "pbc.h"
#include "force.h"
#include "fatal.h"
#include "futil.h"
#include "maths.h"
#include "macros.h"
#include "physics.h"
#include "vec.h"
#include "tpxio.h"
#include "mdrun.h"
#include "calcpot.h"
#include "main.h"
#include "random.h"

static int *mk_windex(int w1,int nw)
{
  int *index;
  int i;

  snew(index,nw);
  for(i=0; (i<nw); i++)
    index[i]=w1+3*i;
  
  return index;
}

static void insert_ion(real q,int *nwater,bool bSet[],int index[],
		       real pot[],rvec x[],matrix box,
		       char *anm,char *resnm,
		       t_topology *top,real rmin,bool bRandom,int *seed)
{
  int  i,ii,ei,owater,wlast,m,nw;
  real extr_e,poti,rmin2;
  rvec xei,dx;
  bool bSub=FALSE;
  int  maxrand,natoms;
  
  ei=-1;
  nw = *nwater;
  maxrand = 100*nw;
  if (bRandom) {
    do {
      ei = nw*rando(seed);
      maxrand--;
    } while (bSet[ei] && (maxrand > 0));
  }
  else {
    for(i=0; (i<nw); i++) {
      if (!bSet[i]) {
	ii=index[i];
	poti=pot[ii];
	if (q > 0) {
	  if ((poti <= extr_e) || !bSub) {
	    extr_e = poti;
	    ei     = i;
	    bSub   = TRUE;
	  }
	}
	else {
	  if ((poti >= extr_e) || !bSub) {
	    extr_e = poti;
	    ei     = i;
	    bSub   = TRUE;
	  } 
	}
      }
    }
  }
  if (ei == -1)
    fatal_error(0,"No More replacable waters!");
  fprintf(stderr,"Replacing water %d (atom %d)\n",ei,index[ei]);
    
  /* Swap ei water at index ei with water at last index */
  owater = index[ei];
  wlast  = index[nw-1];
  for(m=0; (m<3); m++) {
    copy_rvec(x[owater+m],xei);
    copy_rvec(x[wlast+m],x[owater+m]);
    copy_rvec(xei,x[wlast+m]);
    pot[owater+m] = pot[wlast+m];
  }
  
  nw--;
  *nwater = nw;
    
  /* Mark all waters within rmin as unavailable for substitution */
  if (rmin > 0) {
    copy_rvec(x[wlast],xei);
    rmin2=rmin*rmin;
    for(i=0; (i<nw); i++) {
      if (!bSet[i]) {
	pbc_dx(xei,x[index[i]],dx);
	if (iprod(dx,dx) < rmin2)
	  bSet[i] = TRUE;
      }
    }
  }
  /* Replace the oxygen with the ion */
  natoms=top->atoms.nr;
  replace_atom(top,wlast,anm,resnm,q,1,1);
  delete_atom(top,wlast+2);
  delete_atom(top,wlast+1);
  /* Now copy all the remaining coordinates two places down the x array */
  for(i=wlast+3; (i<natoms); i++)
    copy_rvec(x[i],x[i-2]);
}

static void copy_atom(t_atoms *at1,int a1,t_atoms *at2,int a2,int r)
{
  int r1,r2;

  at2->atom[a2]=at1->atom[a1];
  snew(at2->atomname[a2],1);
  *at2->atomname[a2]=strdup(*at1->atomname[a1]);
  r1=at1->atom[a1].resnr;
  if (r == -1)
    r2=r1;
  else
    r2=r;
  at2->atom[a2].resnr=r2;
  if (at2->resname[r2] == NULL) {
    snew(at2->resname[r2],1);
    *at2->resname[r2]=strdup(*at1->resname[r1]);
  }
}  

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "genion calculates the electrostatic potential on all atoms, using",
    "normal GROMACS particle based methods (in contrast to other methods",
    "based on solving the Poisson-Boltzmann equation).",
    "If specified in the [TT]tpr[tt] file, a reaction field or",
    "generalized reaction field can be used, check out the manual of",
    "grompp for more information. The potential can be written as B-factors"
    "in a pdb file (for visualisation using e.g. rasmol)[PAR]",
    "When the potential has been calculated, ions can be generated",
    "at the positions of water molecules. This can be used to simulate",
    "a solution with a given ionic strength, or to neutralize the simulation",
    "system which is necessary when using e.g. the PPPM electrostatics",
    "method in the simulation. Note that PPPM can not be used in the",
    "genion program[PAR]",
    "To avoid",
    "clustering of ions it is advisable to set rmin (the radius around",
    "an ion in which no other ion will be placed) to a value high enough",
    "to allow solvent around the ion (> 0.6 nm)."
  };
  static char *bugs[] = {
    "Only monatomic ions can be used. For larger ions, e.g. sulfate we recommended to use genbox.",
    "The rmin option is currently out of order"
  };
  static int  p_num=0,n_num=0,nw=0,w1=1;
  static char *p_name="Na",*n_name="Cl";
  static real p_q=1,n_q=-1,rcut=0;
  static int  seed=1993;
  static bool bRandom=FALSE;
  static t_pargs pa[] = {
    { "-p",    FALSE, etINT,  &p_num, "Number of positive ions"       },
    { "-pn",   FALSE, etSTR,  &p_name,"Name of the positive ion"      },
    { "-pq",   FALSE, etREAL, &p_q,   "Charge of the positive ion"    },
    { "-n",    FALSE, etINT,  &n_num, "Number of negative ions"       },
    { "-nn",   FALSE, etSTR,  &n_name,"Name of the negative ion"      },
    { "-nq",   FALSE, etREAL, &n_q,   "Charge of the negative ion"    },
    { "-rmin", FALSE, etREAL, &rcut,  "Minimum distance between ions" },
    { "-w1",   FALSE, etINT,  &w1,    "First water atom to be cosidered" },
    { "-nw",   FALSE, etINT,  &nw,    "Number of water molecules" },
    { "-random",FALSE,etBOOL, &bRandom,"Use random placement of waters instead of based on potential. The rmin option should still work" },
    { "-seed", FALSE, etINT,  &seed,  "Seed for random number generator" }
  };
  t_topology  *top;
  t_parm      parm;
  t_commrec   cr;
  t_mdatoms   *mdatoms;
  t_nsborder  nsb;
  t_groups    grps;
  t_graph     *graph;
  t_forcerec  *fr;
  rvec        *x,*v;
  real        *pot;
  matrix      box;
  int         *index;
  bool        *bSet,bPDB;
  int         i;
  t_filenm fnm[] = {
    { efTPX, NULL,  NULL,     ffREAD  },
    { efSTO, "-o",  NULL,     ffWRITE },
    { efLOG, "-g",  "genion", ffWRITE },
    { efPDB, "-q",  "conf",   ffOPTWR },
    { efXTC, "-x",  "ctraj",  ffOPTRD },
    { efTRN, "-t",  "traj",   ffOPTRD },
    { efENX, "-e",  "ener",   ffOPTRD }
  };
#define NFILE asize(fnm)
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,TRUE,NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    asize(bugs),bugs);
  bPDB = ftp2bSet(efPDB,NFILE,fnm);
  if (bRandom && bPDB) {
    fprintf(stderr,"Not computing potential with random option!\n");
    bPDB = FALSE;
  }
    
  /* Check input for something sensible */
  w1--;
  if (nw <= 0) {
#define NOWAT "No water molecules in your system."
    if ((p_num != 0) || (n_num != 0))
      fatal_error(0,NOWAT" Can not generate ions");
    else if (!bPDB)
      fatal_error(0,NOWAT" And you don't want to know the potential either?");
    else
      fprintf(stderr,NOWAT"Will calculate potential anyway\n");
  }
  else 
    fprintf(stderr,"First water atom: %d  Number of water molecules: %d\n",
	    w1+1,nw);
  
  index  = mk_windex(w1,nw);
  snew(bSet,nw);
  snew(top,1);
  init_calcpot(NFILE,fnm,top,&x,&parm,&cr,
	       &graph,&mdatoms,&nsb,&grps,&fr,&pot,box);
  snew(v,top->atoms.nr);
  snew(top->atoms.pdbinfo,top->atoms.nr);
  
  /* Calculate potential once only */
  if (!bRandom) {
    calc_pot(stdlog,&nsb,&cr,&grps,&parm,top,x,fr,graph,mdatoms,pot);
    if (bPDB) {
      for(i=0; (i<top->atoms.nr); i++)
	top->atoms.pdbinfo[i].bfac = pot[i];
      write_sto_conf(ftp2fn(efPDB,NFILE,fnm),"Potential calculated by genion",
		     &top->atoms,x,v,box);
    }
  }
  /* Now loop over the ions that have to be placed */
  if ((p_num > 0) || (n_num > 0)) {
    do {
      if ((p_num > 0) && (p_num >= n_num))  {
	fprintf(stderr,"+");
	insert_ion(p_q,&nw,bSet,index,pot,x,box,
		   p_name,p_name,top,rcut,bRandom,&seed);
	p_num--;
      }
      else if (n_num > 0) {
	fprintf(stderr,"-");
	insert_ion(n_q,&nw,bSet,index,pot,x,box,
		   n_name,n_name,top,rcut,bRandom,&seed);
	n_num--;
      }
    } while (p_num+n_num > 0);
    fprintf(stderr,"\n");
    
    write_sto_conf(ftp2fn(efSTO,NFILE,fnm),*top->name,&top->atoms,x,v,box);
  }
  thanx(stdout);
  
  return 0;
}
