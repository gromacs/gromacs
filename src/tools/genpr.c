/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_genpr_c = "$Id$";

#include <math.h>
#include "sysstuff.h"
#include "statutil.h"
#include "string.h"
#include "copyrite.h"
#include "smalloc.h"
#include "typedefs.h"
#include "confio.h"
#include "futil.h"
#include "macros.h"
#include "rdgroup.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "genpr produces an include file for a topology containing",
    "a list of atom numbers and three force constants for the",
    "X, Y and Z direction. A single isotropic force constant may",
    "be given on the command line instead of three components.[PAR]",
    "WARNING: genpr only works for the first molecule.",
    "Position restraints are interactions within molecules, therefore",
    "they should be included within the correct [TT][ moleculetype ][tt]",
    "block in the topology. Since the atom numbers in every moleculetype",
    "in the topology start at 1 and the numbers in the input file for",
    "genpr number consecutively from 1, genpr will only produce a useful",
    "file for the first molecule."
  };
  static rvec    fc={1000.0,1000.0,1000.0};
  t_pargs pa[] = {
    { "-fc", FALSE, etRVEC, {&fc}, 
      "force constants (kJ mol-1 nm-2)" }
  };
  
  t_atoms atoms;
  int     i;
  FILE    *out;
  int          igrp;
  atom_id      *ind_grp;
  char         *gn_grp;
  char         title[STRLEN];
  matrix       box;
  
  t_filenm fnm[] = {
    { efSTX, "-f", NULL, ffREAD },
    { efNDX, "-n", NULL, ffOPTRD },
    { efITP, "-o", "posre", ffWRITE }
  };
#define NFILE asize(fnm)
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);
  
  if ( !ftp2bSet(efNDX,NFILE,fnm) ) {
    if ( !ftp2bSet(efSTX,NFILE,fnm) )
      fatal_error(0,"no index file and no structure file suplied");
    else {
      rvec *x,*v;
      
      get_stx_coordnum(ftp2fn(efSTX,NFILE,fnm),&(atoms.nr));
      init_t_atoms(&atoms,atoms.nr,TRUE);
      snew(x,atoms.nr);
      snew(v,atoms.nr);
      fprintf(stderr,"\nReading structure file\n");
      read_stx_conf(ftp2fn(efSTX,NFILE,fnm),title,&atoms,x,v,box);
      sfree(x);
      sfree(v);
    }
  }
  printf("Select group to position restrain\n");
  get_index(&atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&igrp,&ind_grp,&gn_grp);
  
  out=ftp2FILE(efITP,NFILE,fnm,"w");
  fprintf(out,"; position restraints for %s of %s\n\n",gn_grp,title);
  fprintf(out,"[ position_restraints ]\n");
  fprintf(out,";%3s %5s %9s %10s %10s\n","i","funct","fcx","fcy","fcz");
  for(i=0; i<igrp; i++) 
    fprintf(out,"%4u %4d %10g %10g %10g\n",
	    ind_grp[i]+1,1,fc[XX],fc[YY],fc[ZZ]);
  fclose(out);
  
  thanx(stdout);
  
  return 0;
}
