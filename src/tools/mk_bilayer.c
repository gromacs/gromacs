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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
#include <math.h>
#include "smalloc.h"
#include "typedefs.h"
#include "macros.h"
#include "string2.h"
#include "confio.h"
#include "vec.h"
#include "statutil.h"
#include "copyrite.h"
#include "3dview.h"

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "make_bilayer reads a monolayer and turns it into a bilayer."
    "The distance between the two can be set at the command line."
  };
  t_manual man = { asize(desc),desc,0,NULL,NULL,0,NULL};

  real    tr=0;
  char    title[STRLEN];
  int     i,n,natoms;
  rvec    *x,*v,min;
  mat4    rot1,tran0,tran1,tran2,total,tmp1,tmp2;
  vec4    *xrot;
  matrix  box;
  t_atoms atoms;
  t_filenm fnm[] = {
    { efGRO, "-f", NULL,      ffREAD  },
    { efGRO, "-o", "confout", ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,FALSE,&man);

  for(i=1; (i<argc); i++)
    if (strcmp(argv[i],"-tr") == 0) 
      tr=dscan(argc,argv,&i);
    else
      usage(argv[0],argv[i]);
  
  atoms.nr=0;
  atoms.nres=0;
  get_coordnum(fnm[0].fn,&natoms);
  snew(atoms.atomname,2*natoms);
  snew(atoms.resname,2*natoms);
  snew(atoms.atom,2*natoms);
  snew(x,2*natoms);
  snew(v,2*natoms);
  snew(xrot,natoms);

  read_whole_conf(fnm[0].fn,title,&atoms,x,v,box);

  copy_rvec(x[0],min);;

  for (i=0;i<natoms;i++)
  {
    if (x[i][XX] < min[XX]) min[XX] = x[i][XX];
    if (x[i][YY] < min[YY]) min[YY] = x[i][YY];
    if (x[i][ZZ] < min[ZZ]) min[ZZ] = x[i][ZZ];
  }

  translate(0,0,min[ZZ],tran1); 
  translate(-box[XX][XX],0,tr,tran2);
  rotate(YY, M_PI, rot1);

  mult_matrix(tmp1,tran2,rot1);
  mult_matrix(total,tmp1,tran1);

  for(i=0; i<natoms;i++)
    m4_op(total,x[i],xrot[i]);
  
  for(i=0; i<natoms; i++)
    for(n=0;n<DIM;n++)
      {
	x[natoms+i][n] = xrot[i][n];
	v[natoms+i][n] = 0.00;
      }

  for (i=natoms;i<2*natoms;i++)
  {  
    atoms.atomname[i] = atoms.atomname[i-natoms];  
    atoms.resname[i]  = atoms.resname[i-natoms];
    atoms.atom[i] = atoms.atom[i-natoms];
  }
  
  atoms.nr *= 2;
  atoms.nres *= 2;

  write_conf(fnm[1].fn,title,&atoms,x,v,box);

  thanx(stdout);

  return 0;
}


