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
 * Great Red Owns Many ACres of Sand 
 */
#include "statutil.h"
#include "macros.h"
#include "readinp.h"

int main(int argc,char *argv[])
{
  t_filenm fnm[] = {
    { efMDP, "-f",  NULL,    FALSE },
    { efMDP, "-po", "mdout", FALSE },
  };
#define NFILE asize(fnm)
  int       ninp;
  t_inpfile *inp;
  static char *boxes[] = { "cubic", "rectangular", NULL };
  
  /* Parse the command line */
  parse_common_args(&argc,argv,0,NFILE,fnm,FALSE,NULL);
  
  inp=read_inpfile(opt2fn("-f",NFILE,fnm),&ninp);
  write_inpfile(opt2fn("-po",NFILE,fnm),ninp,inp);

  printf("nsteps = %d\n",get_eint(ninp,inp,"nsteps",-1));
  printf("dt     = %g\n",get_ereal(ninp,inp,"dt",0.1));
  printf("cpp    = %s\n",get_estr(ninp,inp,"cpp","/lib/cpp"));
  printf("box    = %s\n",boxes[get_eenum(ninp,inp,"box",boxes)]);
    
  return 0;
}
