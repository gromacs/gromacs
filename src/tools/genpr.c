/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Glycine aRginine prOline Methionine Alanine Cystine Serine
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
    { "-fc", FALSE, etRVEC, {fc}, 
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
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
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
    fprintf(out,"%4d %4d %10g %10g %10g\n",
	    ind_grp[i]+1,1,fc[XX],fc[YY],fc[ZZ]);
  fclose(out);
  
  thanx(stderr);
  
  return 0;
}
