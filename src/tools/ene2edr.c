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
    "ene2edr converts ene files to edr files which are in portable",
    "xdr format. These files can be read on any platform by g_energy.",
    "This is useful when you perform your mdruns on a different",
    "CPU than you use for analysis."
  };

  FILE      *in;
  t_energy  *ee;
  t_drblock *dr;
  XDR       xdr;
  int       nre,j;
  real      t;
  char      **enm;

  t_filenm fnm[] = {
    { efENE, "-f", "ener", ffREAD },
    { efEDR, "-o", "enxdr",ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,0,NULL,asize(desc),desc,
		    0,NULL);

  in=ftp2FILE(efENE,NFILE,fnm,"r");
  
  dr=get_drblock();
  rd_ener_nms(in,&nre,&enm);
  snew(ee,nre);
  xdropen(&xdr,ftp2fn(efEDR,NFILE,fnm),"w");
  edr_nms(&xdr,&nre,&enm);
  
  while (rd_ener(in,&t,&j,ee,dr))  {
    if ((j%10) == 0)
      fprintf(stderr,"\rFrame: %d",j);
    edr_io(&xdr,&t,&j,&nre,ee,dr);
    j++;
  }
  fprintf(stderr,"\n");
    
  thanx(stdout);
  
  return 0;
}
