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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
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
#include "index.h"

static void insert_ion(int nsa,int *nwater,
		       bool bSet[],int repl[],atom_id index[],
		       real pot[],rvec x[],t_pbc *pbc,
		       int sign,real q,char *ionname,
		       t_mdatoms *mdatoms,
		       real rmin,bool bRandom,int *seed)
{
  int  i,ii,ei,owater,wlast,m,nw;
  real extr_e,poti,rmin2;
  rvec xei,dx;
  bool bSub=FALSE;
  int  maxrand;
  
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
    extr_e = 0;
    for(i=0; (i<nw); i++) {
      if (!bSet[i]) {
	ii=index[nsa*i];
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
    fatal_error(0,"No more replaceable solvent!");
  fprintf(stderr,"Replacing solvent molecule %d (atom %d) with %s\n",
	  ei,index[nsa*ei],ionname);
  
  /* Replace solvent molecule charges with ion charge */
  bSet[ei] = TRUE;
  repl[ei] = sign;
  mdatoms->chargeA[index[nsa*ei]] = q;
  for(i=1; i<nsa; i++)
    mdatoms->chargeA[index[nsa*ei+i]] = 0;

  /* Mark all solvent molecules within rmin as unavailable for substitution */
  if (rmin > 0) {
    rmin2=rmin*rmin;
    for(i=0; (i<nw); i++) {
      if (!bSet[i]) {
	pbc_dx(pbc,x[index[nsa*ei]],x[index[nsa*i]],dx);
	if (iprod(dx,dx) < rmin2)
	  bSet[i] = TRUE;
      }
    }
  }
}

static char *aname(char *mname)
{
  char *str;
  int  i;

  str = strdup(mname);
  i=strlen(str)-1;
  while (i>1 && (isdigit(str[i]) || (str[i]=='+') || (str[i]=='-'))) {
    str[i]='\0';
    i--;
  }

  return str;
}

void sort_ions(int nsa,int nw,int repl[],atom_id index[],
	       t_atoms *atoms,rvec x[],
	       char *p_name,char *n_name)
{
  int i,j,k,r,np,nn,starta,startr,npi,nni;
  rvec *xt;
  char **pptr=NULL,**nptr=NULL,**paptr=NULL,**naptr=NULL;

  snew(xt,atoms->nr);

  /* Put all the solvent in front and count the added ions */
  np=0;
  nn=0;
  j=index[0];
  for(i=0; i<nw; i++) {
    r = repl[i];
    if (r == 0)
      for(k=0; k<nsa; k++)
	copy_rvec(x[index[nsa*i+k]],xt[j++]);
    else if (r>0)
      np++;
    else if (r<0)
      nn++;
  }

  if (np+nn > 0) {
    /* Put the positive and negative ions at the end */
    starta = index[nsa*(nw - np - nn)];
    startr = atoms->atom[starta].resnr;

    if (np) {
      snew(pptr,1);
      pptr[0] = p_name;
      snew(paptr,1);
      paptr[0] = aname(p_name);
    }
    if (nn) {
      snew(nptr,1);
      nptr[0] = n_name;
      snew(naptr,1);
      naptr[0] = aname(n_name);
    }
    npi = 0;
    nni = 0;
    for(i=0; i<nw; i++) {
      r = repl[i];
      if (r > 0) {
	j = starta+npi;
	k = startr+npi;
	copy_rvec(x[index[nsa*i]],xt[j]);
	atoms->atomname[j] = paptr;
	atoms->atom[j].resnr = k ;
	atoms->resname[k] = pptr;
	npi++;
      } else if (r < 0) {
	j = starta+np+nni;
	k = startr+np+nni;
	copy_rvec(x[index[nsa*i]],xt[j]);
	atoms->atomname[j] = naptr;
	atoms->atom[j].resnr = k;
	atoms->resname[k] = nptr;
	nni++;
      }
    }
    for(i=index[nsa*nw-1]+1; i<atoms->nr; i++) {
      j = i-(nsa-1)*(np+nn);
      atoms->atomname[j] = atoms->atomname[i];
      atoms->atom[j] = atoms->atom[i];
      copy_rvec(x[i],xt[j]);
    }
    atoms->nr -= (nsa-1)*(np+nn);

    /* Copy the new positions back */
    for(i=index[0]; i<atoms->nr; i++)
      copy_rvec(xt[i],x[i]);
    sfree(xt);
  }
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "genion replaces solvent molecules by monoatomic ions at",
    "the position of the first atoms with the most favorable electrostatic",
    "potential or at random. The potential is calculated on all atoms, using",
    "normal GROMACS particle based methods (in contrast to other methods",
    "based on solving the Poisson-Boltzmann equation).",
    "The potential is recalculated after every ion insertion.",
    "If specified in the run input file, a reaction field, shift function",
    "or user function can be used. For the user function a table file",
    "can be specified with the option [TT]-table[tt].",
    "The group of solvent molecules should be continuous and all molecules",
    "should have the same number of atoms.",
    "The user should add the ion molecules to the topology file and include",
    "the file [TT]ions.itp[tt].",
    "Ion names for Gromos96 should include the charge.[PAR]",
    "With the option [TT]-pot[tt] the potential can be written as B-factors",
    "in a pdb file (for visualisation using e.g. rasmol).",
    "The unit of the potential is 1000 kJ/(mol e), the scaling be changed",
    "with the [TT]-scale[tt] option.[PAR]",
    "For larger ions, e.g. sulfate we recommended to use genbox."
  };
  static int  p_num=0,n_num=0;
  static char *p_name="Na",*n_name="Cl";
  static real p_q=1,n_q=-1,rmin=0.6,scale=0.001;
  static int  seed=1993;
  static bool bRandom=FALSE;
  static t_pargs pa[] = {
    { "-np",    FALSE, etINT,  {&p_num}, "Number of positive ions"       },
    { "-pname", FALSE, etSTR,  {&p_name},"Name of the positive ion"      },
    { "-pq",    FALSE, etREAL, {&p_q},   "Charge of the positive ion"    },
    { "-nn",    FALSE, etINT,  {&n_num}, "Number of negative ions"       },
    { "-nname", FALSE, etSTR,  {&n_name},"Name of the negative ion"      },
    { "-nq",    FALSE, etREAL, {&n_q},   "Charge of the negative ion"    },
    { "-rmin",  FALSE, etREAL, {&rmin},  "Minimum distance between ions" },
    { "-random",FALSE,etBOOL, {&bRandom},"Use random placement of ions instead of based on potential. The rmin option should still work" },
    { "-seed",  FALSE, etINT,  {&seed},  "Seed for random number generator" },
    { "-scale", FALSE, etREAL, {&scale}, "Scaling factor for the potential for -pot" }
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
  t_pbc       pbc;
  int         *repl;
  atom_id     *index;
  char        *grpname;
  bool        *bSet,bPDB;
  int         i,nw,nwa,nsa;
  t_filenm fnm[] = {
    { efTPX, NULL,  NULL,      ffREAD  },
    { efXVG, "-table","table", ffOPTRD },
    { efNDX, NULL,  NULL,      ffOPTRD },
    { efSTO, "-o",  NULL,      ffWRITE },
    { efLOG, "-g",  "genion",  ffWRITE },
    { efPDB, "-pot", "pot",    ffOPTWR }
  };
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_BE_NICE,NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    0,NULL);
  bPDB = ftp2bSet(efPDB,NFILE,fnm);
  if (bRandom && bPDB) {
    fprintf(stderr,"Not computing potential with random option!\n");
    bPDB = FALSE;
  }
    
  /* Check input for something sensible */
  if ((p_num<0) || (n_num<0))
    fatal_error(0,"Negative number of ions to add?");

  nsb.nodeid=0;

  snew(top,1);
  init_calcpot(ftp2fn(efLOG,NFILE,fnm),ftp2fn(efTPX,NFILE,fnm),
	       opt2fn("-table",NFILE,fnm),top,&parm,&cr,
	       &graph,&mdatoms,&nsb,&grps,&fr,&pot,box,&x);

  if ((p_num == 0) && (n_num == 0)) {
    if (!bPDB) {
      fprintf(stderr,"No ions to add and no potential to calculate.\n");
      exit(0);
    }
    nw  = 0;
    nsa = 0; /* too keep gcc happy */
  } else {
    printf("Select a continuous group of solvent molecules\n");
    get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&nwa,&index,&grpname);
    for(i=1; i<nwa; i++)
      if (index[i] != index[i-1]+1)
	fatal_error(0,"The solvent group is not continuous: index[%d]=%d, "
		    "index[%d]=%d",i,index[i-1]+1,i+1,index[i]+1);
    nsa = 1;
    while ((nsa<nwa) &&
	   (top->atoms.atom[index[nsa]].resnr ==
	    top->atoms.atom[index[nsa-1]].resnr))
      nsa++;
    if (nwa % nsa)
      fatal_error(0,"Your solvent group size (%d) is not a multiple of %d",
		  nwa,nsa);
    nw = nwa/nsa;
    fprintf(stderr,"Number of (%d-atomic) solvent molecules: %d\n",nsa,nw);
	if (p_num+n_num > nw)
      fatal_error(0,"Not enough solvent for adding ions");
  }
  snew(bSet,nw);
  snew(repl,nw);
  
  snew(v,top->atoms.nr);
  snew(top->atoms.pdbinfo,top->atoms.nr);

  set_pbc(&pbc,box);

  /* Now loop over the ions that have to be placed */
  do {
    if (!bRandom) {
      calc_pot(stdlog,&nsb,&cr,&grps,&parm,top,x,fr,mdatoms,pot,box,graph);
      if (bPDB || debug) {
	char buf[STRLEN];
	
	if (debug)
	  sprintf(buf,"%d_%s",p_num+n_num,ftp2fn(efPDB,NFILE,fnm));
	else
	  strcpy(buf,ftp2fn(efPDB,NFILE,fnm));
	for(i=0; (i<top->atoms.nr); i++)
	    top->atoms.pdbinfo[i].bfac = pot[i]*scale;
	write_sto_conf(buf,"Potential calculated by genion",
		       &top->atoms,x,v,box);
	bPDB = FALSE;
      }
    }
    if ((p_num > 0) && (p_num >= n_num))  {
      insert_ion(nsa,&nw,bSet,repl,index,pot,x,&pbc,
		 1,p_q,p_name,mdatoms,rmin,bRandom,&seed);
      p_num--;
    }
    else if (n_num > 0) {
      insert_ion(nsa,&nw,bSet,repl,index,pot,x,&pbc,
		 -1,n_q,n_name,mdatoms,rmin,bRandom,&seed);
      n_num--;
    }
  } while (p_num+n_num > 0);
  fprintf(stderr,"\n");

  if (nw)
    sort_ions(nsa,nw,repl,index,&top->atoms,x,p_name,n_name);
  
  sfree(top->atoms.pdbinfo);
  top->atoms.pdbinfo = NULL;
  write_sto_conf(ftp2fn(efSTO,NFILE,fnm),*top->name,&top->atoms,x,NULL,box);

  thanx(stderr);
  
  return 0;
}
