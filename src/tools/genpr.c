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
 * Gyas ROwers Mature At Cryogenic Speed
 */
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

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "genpr produces an include file for a topology containing",
    "a list of atom numbers and three force constants for the",
    "X, Y and Z direction. A single isotropic force constant may",
    "be given on the command line instead of three components.[PAR]",
    "This list is used as the position restraint list"
  };
  static int     a1=-1,a2=-1;
  static real    fc=1000,fx=1000,fy=1000,fz=1000;
  t_pargs pa[] = {
    { "-fc", FALSE, etREAL, &fc, 
      "force constant for isotropic restraining (kJ/mol nm^2)" },
    { "-fx", FALSE, etREAL, &fx, 
      "id. for X direction (id)" },
    { "-fy", FALSE, etREAL, &fy, 
      "id. for Y direction (id)" },
    { "-fz", FALSE, etREAL, &fz, 
      "id. for Z direction (id)" },
    { "-a1", FALSE, etINT,  &a1, "first atom (numbering from 1)" },
    { "-a2", FALSE, etINT,  &a2, "last atom" }
  };
  int     i;
  FILE    *out;
  t_filenm fnm[] = {
    { efITP, "-o", "posre", ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  if ((a1 == -1) || (a2 == -1)) 
    fatal_error(0,"a1 (%d) or a2 (%d) not set",a1,a2);
  if (opt2parg_bSet("-fc",asize(pa),pa))
    fx=fy=fz=fc;
      
  out=ftp2FILE(efITP,NFILE,fnm,"w");
  fprintf(out,"[ position_restraints ]\n");
  fprintf(out,";%7s%8s%8s\n","i","funct","fc");
  for(i=a1; (i<=a2); i++) 
    fprintf(out,"%8d%8d  %8.0f  %8.0f  %8.0f\n",
	    i,1,fx,fy,fz);
  fclose(out);
  
  thanx(stdout);
  
  return 0;
}
