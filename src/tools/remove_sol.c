/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_remove_sol_c = "$Id$";

#include "copyrite.h"
#include "string2.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "confio.h"
#include "paramio.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include <physics.h>

/* center configuration to 0,,0,0 */
void move_to_origin(int natoms,rvec *x)
{
  int n,d;
  rvec xmin;
  
  /* determine xmin */
  for(d=0;(d<DIM);d++) {
    xmin[d]=x[0][d];
    for(n=0;(n<natoms);n++) {
      if (x[n][d]<xmin[d])
	xmin[d]=x[n][d];
    }
  }
  
  /* translate */
  for(n=0;(n<natoms);n++)
    for(d=0;(d<DIM);d++) 
      x[n][d]-=xmin[d];
}


void usage(char *prog,char *arg)
{
  fprintf(stderr,"Usage: %s %s \n",prog,common_args);
  exit(1);
}

void main(int argc, char *argv[])
{
  t_atoms *atoms,*atoms_o;
  int     natoms,natoms_o,resnr_o;
  rvec    *x,*v,*x_o,*v_o;
  char    title[STRLEN];
  matrix  box;
  bool    *bRemove;
  bool    bNewRes;
  real     MIN_X,MAX_X;
  char *env;
  
  int     n,i,d;

  t_filenm fnm[] = {
    { efGRO, "-f", NULL,      ffREAD  },
    { efGRO, "-o", "output",  ffWRITE }
  };
#define NFILE asize(fnm)

  /* parse it */
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,FALSE,NULL);

  /* read protein and determine size */
  get_coordnum(opt2fn("-f",NFILE,fnm),&natoms);
  for(i=0;(i<natoms);i++)
    bRemove=FALSE;
  snew(atoms,1);
  snew(atoms->resname,natoms);
  snew(atoms->atomname,natoms);
  snew(atoms->atom,natoms);
  snew(x,natoms);
  snew(v,natoms);
  snew(atoms_o,1);
  snew(atoms_o->resname,natoms);
  snew(atoms_o->atomname,natoms);
  snew(atoms_o->atom,natoms);
  snew(x_o,natoms);
  snew(v_o,natoms);
  natoms_o=0;
  read_whole_conf(opt2fn("-f",NFILE,fnm),title,atoms,x,v,box);

  /* perform geometric transformation */


  /* move molecule in such a way that the minimum values are at 0,0,0 */
  env = (char *)getenv("MIN_X");
  sscanf(env,"%f",&MIN_X);
  env = (char *)getenv("MAX_X");
  sscanf(env,"%f",&MAX_X);


  fprintf(stderr," Min X = %8.3f\n",MIN_X);
  fprintf(stderr," Max X = %8.3f\n",MAX_X);

  snew(bRemove,atoms->nres);
  for(n=0;(n<atoms->nr);n++)  {
    if ( n % 10  == 0 ) 
      fprintf(stderr,"\r %5d",n);
    if ( strcmp("SOL",*atoms->resname[atoms->atom[n].resnr])==0) {
      if (( x[n][XX] < MIN_X ) || (x[n][XX] > MAX_X ))
	bRemove[atoms->atom[n].resnr]=TRUE;
    }
  }
  fprintf(stderr,"asignment done\n");

  natoms_o=0;
  resnr_o=0;
  atoms_o->nr=0;
  atoms_o->nres=0;
  for(n=0;(n<atoms->nr);n++) {
    if ( n % 10 == 0 ) 
      fprintf(stderr,"\r %5d",n);
    if ( n==0 )
      bNewRes=TRUE;
    else {
      if ( atoms->atom[n].resnr != atoms->atom[n-1].resnr )
	bNewRes=TRUE;
      else
	bNewRes=FALSE;
    }
    if ( bRemove[atoms->atom[n].resnr]==FALSE) {
      if ( bNewRes == TRUE )
	atoms_o->nres++;
      atoms_o->nr++;
      for(d=0;(d<DIM);d++) {
	x_o[natoms_o][d]=x[n][d];
	v_o[natoms_o][d]=v[n][d];
      }
      /* copy resname */
      atoms_o->atomname[natoms_o]   = atoms->atomname[n];
      atoms_o->atom[natoms_o].resnr = atoms_o->nres-1;
      atoms_o->resname[atoms_o->atom[natoms_o].resnr]=atoms->resname[atoms->atom[n].resnr];
      natoms_o++;
	
    }
  }

  /* write away configuration */
  write_conf(opt2fn("-o",NFILE,fnm),title,atoms_o,x_o,v_o,box);
  

  thanx(stdout);
}
