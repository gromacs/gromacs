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
 * GRowing Old MAkes el Chrono Sweat
 */
#include <math.h>
#include <string.h>
#include "pdbio.h"
#include "confio.h"
#include "symtab.h"
#include "smalloc.h"
#include "macros.h"
#include "copyrite.h"
#include "statutil.h"

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "pdb2gro converts a protein data bank file to a gromos coordinate",
    "file. Angstroms are converted to nanometers and atoms and residues",
    "are renumbered. No box is generated (zero box)"
  };

  FILE      *in,*out;
  int       natom;
  t_pdbatom *pdba;
  t_atoms   atoms;
  t_symtab  symtab;
  rvec      *x,*v;
  matrix    box;

  static bool      bH=FALSE;
  t_pargs pa[] = {
    { "-nh", FALSE, etBOOL, &bH, "Filter out hydrogen atoms." }
  };
  t_filenm fnm[] = {
    { efPDB, "-f", NULL, ffREAD  },
    { efGRO, "-o", NULL, ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);
  
  in=ftp2FILE(efPDB,NFILE,fnm,"r");
  natom=read_pdbatoms(in,&pdba,box,bH);
  fclose(in);
  open_symtab(&symtab);
  pdb2atoms(natom,pdba,&atoms,&x,&symtab);
  snew(v,natom);
  out=ftp2FILE(efGRO,NFILE,fnm,"w");
  write_hconf(out,cool_quote(),&atoms,x,v,box);
  fclose(out);
  
  thanx(stdout);
  
  return 0;
}
