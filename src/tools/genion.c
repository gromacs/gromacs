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

static void insert_ion(int *nwater,bool bSet[],int repl[],int index[],
		       real pot[],rvec x[],
		       int sign,real q,char *ionname,
		       t_mdatoms *mdatoms,
		       real rmin,bool bRandom,int *seed)
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
    extr_e = 0;
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
    fatal_error(0,"No More replaceable waters!");
  fprintf(stderr,"Replacing water %d (atom %d) with %s\n",
	  ei,index[ei],ionname);
  
  /* Replace water charges with ion charge */
  bSet[ei] = TRUE;
  repl[ei] = sign;
  mdatoms->chargeA[index[ei]] = q;
  mdatoms->chargeA[index[ei]+1] = 0;
  mdatoms->chargeA[index[ei]+2] = 0;

  /* Mark all waters within rmin as unavailable for substitution */
  if (rmin > 0) {
    rmin2=rmin*rmin;
    for(i=0; (i<nw); i++) {
      if (!bSet[i]) {
	pbc_dx(x[index[ei]],x[index[i]],dx);
	if (iprod(dx,dx) < rmin2)
	  bSet[i] = TRUE;
      }
    }
  }
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

void sort_ions(int nw,int repl[],int index[],t_atoms *atoms,rvec x[],
	       char *p_name,char *n_name)
{
  int i,j,k,r,np,nn,starta,startr,npi,nni;
  rvec *xt;
  char **pptr,**nptr;

  snew(xt,atoms->nr);

  /* Put all the water in front and count the added ions */
  np=0;
  nn=0;
  j=0;
  for(i=0; i<nw; i++) {
    r = repl[i];
    if (r == 0)
      for(k=0; k<3; k++)
	copy_rvec(x[index[i]+k],xt[j++]);
    else if (r>0)
      np++;
    else if (r<0)
      nn++;
  }

  if (np+nn > 0) {
    /* Put the positive and negative ions at the end */
    starta = index[nw - np - nn];
    startr = atoms->atom[starta].resnr;
    if (np) {
      snew(pptr,1);
      pptr[0] = p_name;
    }
    if (nn) {
      snew(nptr,1);
      nptr[0] = n_name;
    }
    npi = 0;
    nni = 0;
    for(i=0; i<nw; i++) {
      r = repl[i];
      if (r > 0) {
	j = starta+npi;
	k = startr+npi;
	copy_rvec(x[index[i]],xt[j]);
	atoms->atomname[j] = pptr;
	atoms->atom[j].resnr = k ;
	atoms->resname[k] = pptr;
	npi++;
      } else if (r < 0) {
	j = starta+np+nni;
	k = startr+np+nni;
	copy_rvec(x[index[i]],xt[j]);
	atoms->atomname[j] = nptr;
	atoms->atom[j].resnr = k;
	atoms->resname[k] = nptr;
	nni++;
      }
    }
    atoms->nr -= 2*(np+nn);
  }
  /* Copy the new positions back */
  for(i=index[0]; i<atoms->nr; i++)
    copy_rvec(xt[i],x[i]);
  sfree(xt);
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "genion replaces water molecules by monoatomic ions. Ions can be placed",
    "at the water oxygen positions with the most favorable electrostatic",
    "potential or at random. The potential is calculated on all atoms, using",
    "normal GROMACS particle based methods (in contrast to other methods",
    "based on solving the Poisson-Boltzmann equation).",
    "The potential is recalculated after every ion insertion.",
    "If specified in the run input file, a reaction field or shift function",
    "can be used. The potential can be written as B-factors",
    "in a pdb file (for visualisation using e.g. rasmol)[PAR]"
    "For larger ions, e.g. sulfate we recommended to use genbox."
  };
  static int  p_num=0,n_num=0,nw=0,w1=1;
  static char *p_name="Na",*n_name="Cl";
  static real p_q=1,n_q=-1,rmin=0.6;
  static int  seed=1993;
  static bool bRandom=FALSE;
  static t_pargs pa[] = {
    { "-p",    FALSE, etINT,  &p_num, "Number of positive ions"       },
    { "-pn",   FALSE, etSTR,  &p_name,"Name of the positive ion"      },
    { "-pq",   FALSE, etREAL, &p_q,   "Charge of the positive ion"    },
    { "-n",    FALSE, etINT,  &n_num, "Number of negative ions"       },
    { "-nn",   FALSE, etSTR,  &n_name,"Name of the negative ion"      },
    { "-nq",   FALSE, etREAL, &n_q,   "Charge of the negative ion"    },
    { "-rmin", FALSE, etREAL, &rmin,  "Minimum distance between ions" },
    { "-w1",   FALSE, etINT,  &w1,    "First water atom to be cosidered (counting from 1)" },
    { "-nw",   FALSE, etINT,  &nw,    "Number of water molecules" },
    { "-random",FALSE,etBOOL, &bRandom,"Use random placement of ions instead of based on potential. The rmin option should still work" },
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
  int         *index,*repl;
  bool        *bSet,bPDB;
  int         i;
  t_filenm fnm[] = {
    { efTPX, NULL,  NULL,      ffREAD  },
    { efSTO, "-o",  NULL,      ffWRITE },
    { efLOG, "-g",  "genion",  ffWRITE },
    { efPDB, "-pot",  "pot",     ffOPTWR }
  };
#define NFILE asize(fnm)
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,TRUE,NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    0,NULL);
  bPDB = ftp2bSet(efPDB,NFILE,fnm);
  if (bRandom && bPDB) {
    fprintf(stderr,"Not computing potential with random option!\n");
    bPDB = FALSE;
  }
    
  /* Check input for something sensible */
  if (p_num+n_num<0)
    fatal_error(0,"No ions to add");
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
  snew(repl,nw);
  snew(top,1);
  init_calcpot(NFILE,fnm,top,&x,&parm,&cr,
	       &graph,&mdatoms,&nsb,&grps,&fr,&pot,box);
  
  snew(v,top->atoms.nr);
  snew(top->atoms.pdbinfo,top->atoms.nr);

  /* Now loop over the ions that have to be placed */
  do {
    if (!bRandom) {
      calc_pot(stdlog,&nsb,&cr,&grps,&parm,top,x,fr,graph,mdatoms,pot);
      if (bPDB || debug) {
	char buf[STRLEN];
	
	if (debug)
	  sprintf(buf,"%d_%s",p_num+n_num,ftp2fn(efPDB,NFILE,fnm));
	else
	  strcpy(buf,ftp2fn(efPDB,NFILE,fnm));
	for(i=0; (i<top->atoms.nr); i++)
	    top->atoms.pdbinfo[i].bfac = pot[i]*0.001;
	write_sto_conf(buf,"Potential calculated by genion",
		       &top->atoms,x,v,box);
	bPDB = FALSE;
      }
    }
    if ((p_num > 0) && (p_num >= n_num))  {
      insert_ion(&nw,bSet,repl,index,pot,x,
		 1,p_q,p_name,mdatoms,rmin,bRandom,&seed);
      p_num--;
    }
    else if (n_num > 0) {
      insert_ion(&nw,bSet,repl,index,pot,x,
		 -1,n_q,n_name,mdatoms,rmin,bRandom,&seed);
      n_num--;
    }
  } while (p_num+n_num > 0);
  fprintf(stderr,"\n");

  sort_ions(nw,repl,index,&top->atoms,x,p_name,n_name);
  
  sfree(top->atoms.pdbinfo);
  top->atoms.pdbinfo = NULL;
  write_sto_conf(ftp2fn(efSTO,NFILE,fnm),*top->name,&top->atoms,x,NULL,box);

  thanx(stdout);
  
  return 0;
}
