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
  static rvec    fc={1000.0,1000.0,1000.0};
  t_pargs pa[] = {
    { "-fc", FALSE, etRVEC, &fc, 
      "force constant (kJ/mol nm^2)" },
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
      
  out=ftp2FILE(efITP,NFILE,fnm,"w");
  fprintf(out,"[ position_restraints ]\n");
  fprintf(out,";%7s%8s%8s\n","i","funct","fc");
  for(i=a1; (i<=a2); i++) 
    fprintf(out,"%8d%8d  %8.0f  %8.0f  %8.0f\n",
	    i,1,fc[XX],fc[YY],fc[ZZ]);
  fclose(out);
  
  thanx(stdout);
  
  return 0;
}
