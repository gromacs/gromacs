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
#include "index.h"

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
    "file for the first molecule.[PAR]",
    "The -of option produces an index file that can be used for",
    "freezing atoms. In this case the input file must be a pdb file."
  };
  static rvec    fc={1000.0,1000.0,1000.0};
  static real    freeze_level;
  t_pargs pa[] = {
    { "-fc", FALSE, etRVEC, {fc}, 
      "force constants (kJ mol-1 nm-2)" },
    { "-freeze", FALSE, etREAL, {&freeze_level},
      "if the -of option or this one is given an index file will be written containing atom numbers of all atoms that have a B-factor less than the level given here" }
  };
  
  t_atoms atoms;
  int     i;
  FILE    *out;
  int          igrp;
  atom_id      *ind_grp;
  char         *gn_grp;
  char         title[STRLEN];
  matrix       box;
  bool         bFreeze;
  
  t_filenm fnm[] = {
    { efSTX, "-f",  NULL,    ffREAD },
    { efNDX, "-n",  NULL,    ffOPTRD },
    { efITP, "-o",  "posre", ffWRITE },
    { efNDX, "-of", "freeze",    ffOPTWR }
  };
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);
  
  bFreeze = opt2bSet("-of",NFILE,fnm) || opt2parg_bSet("-freeze",asize(pa),pa);
  
  if ( !opt2bSet("-n",NFILE,fnm) ) {
    if ( !ftp2bSet(efSTX,NFILE,fnm) )
      gmx_fatal(FARGS,"no index file and no structure file suplied");
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
  if (bFreeze) {
    if (atoms.pdbinfo == NULL) 
      gmx_fatal(FARGS,"No B-factors in input file %s, use a pdb file next time.",
		  ftp2fn(efSTX,NFILE,fnm));
    
    out=opt2FILE("-of",NFILE,fnm,"w");
    fprintf(out,"[ freeze ]\n");
    for(i=0; (i<atoms.nr); i++) {
      if (atoms.pdbinfo[i].bfac <= freeze_level)
	fprintf(out,"%d\n",i+1);
    }
    fclose(out);
  }
  else {
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
  }
  thanx(stderr);
  
  return 0;
}
