/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
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
 * GROtesk MACabre and Sinister
 */
#include "typedefs.h"
#include "txtdump.h"
#include "statutil.h"
#include "smalloc.h"
#include "macros.h"
#include "confio.h"
#include "vec.h"
#include "princ.h"
#include "copyrite.h"

void main(int argc,char *argv[])
{
  int     i,natoms;
  atom_id *index;
  rvec    *x,*v,d,xcm;
  matrix  box,trans;
  char    title[256];
  t_atoms atoms;
  t_filenm fnm[] = {
    { efGRO, "-f", NULL, FALSE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,FALSE,NULL);

  atoms.nr=0;
  atoms.nres=0;
  get_coordnum(fnm[0].fn,&natoms);
  snew(atoms.atomname,natoms);
  snew(atoms.resname,natoms);
  snew(atoms.atom,natoms);
  snew(x,natoms);
  snew(v,natoms);
  read_whole_conf(ftp2fn(efGRO,NFILE,fnm),title,&atoms,x,v,box);
  snew(index,natoms);
  for(i=0; (i<natoms); i++) {
    index[i]=i;
    atoms.atom[i].m=1;
  }
    
  pr_rvecs(stdout,0,"x1",x,natoms);
  sub_xcm(x,natoms,index,atoms.atom,xcm,FALSE);
  /*pr_rvecs(stdout,0,"x2",x,natoms);*/
  principal_comp(natoms,index,atoms.atom,x,trans,d);
  pr_rvecs(stdout,0,"x3",x,natoms);
  rotate_atoms(natoms,index,x,trans);
  /*pr_rvecs(stdout,0,"x4",x,natoms);*/
  add_xcm(x,natoms,index,xcm);
  pr_rvecs(stdout,0,"x5",x,natoms);
}
