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
 * GROningen MAchine for Chemical Simulation
 */
#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "enerio.h"
#include "statutil.h"
#include "disre.h"
#include "names.h"
#include "copyrite.h"
#include "macros.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "edr2ene converts edr files to ene files which are in binary",
    "ene format. These files can be read on only one platform by g_energy.",
    "This is useful when you want to increase performance of analysis"
  };

  FILE      *out;
  t_energy  *ee;
  t_drblock *dr;
  XDR       xdr;
  int       nre,j;
  real      t;
  char      **enm;

  t_filenm fnm[] = {
    { efEDR, "-f", "enexdr", FALSE },
    { efENE, "-o", "ener",FALSE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,0,NULL,
		    asize(desc),desc,0,NULL);


  xdropen(&xdr,ftp2fn(efEDR,NFILE,fnm),"r");
  enm = NULL;
  edr_nms(&xdr,&nre,&enm);

  dr=get_drblock();
  snew(ee,nre);

  out=ftp2FILE(efENE,NFILE,fnm,"w");
  wr_ener_nms(out,nre,enm);


  j=0;
  while ( edr_io(&xdr,&t,&j,&nre,ee,dr)) {
    if ((j%10) == 0)
      fprintf(stderr,"\rFrame: %d",j);
    wr_ener(out,t,j,nre,ee,dr);
    j++;
  }
  fprintf(stderr,"\n");
    
  thanx(stdout);
  
  return 0;
}
