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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_protonate_c = "$Id$";

#include <math.h>
#include "typedefs.h"
#include "macros.h"
#include "copyrite.h"
#include "smalloc.h"
#include "statutil.h"
#include "confio.h"
#include "genhydro.h"

int main (int argc,char *argv[])
{
  static char *desc[] = {
    "protonate protonates a protein molecule."
  };

  char        title[256];  
  t_atoms     *atoms;
  rvec        *x,*v;
  matrix      box;
  int         natoms;
  
  t_filenm fnm[] = {
    { efGRO, "-f", NULL,      ffREAD  },
    { efGRO, "-o", "confout", ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,0,NULL,
		    asize(desc),desc,0,NULL);

  snew(atoms,1);
  atoms->nr=0;
  atoms->nres=0;
  get_coordnum(opt2fn("-f",NFILE,fnm),&natoms);
  snew(atoms->atomname,natoms);
  snew(atoms->resname,natoms);
  snew(atoms->atom,natoms);
  snew(x,natoms);
  snew(v,natoms);
  read_whole_conf(opt2fn("-f",NFILE,fnm),title,atoms,x,v,box);
  sfree(v);
  
  fprintf(stderr,"Read title: %s\n",title);
  
  protonate(&atoms,&x);
  snew(v,atoms->nr);
  
  write_conf(opt2fn("-o",NFILE,fnm),title,atoms,x,v,box);
  
  thanx(stdout);
  
  return 0;
}


