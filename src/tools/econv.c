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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#include "typedefs.h"
#include "smalloc.h"
#include "enerio.h"
#include "statutil.h"
#include "names.h"
#include "copyrite.h"
#include "macros.h"

char *EType[] = { 
  "Zero",
  "Angle4",  
  "Angle", 
  "BuckingHam",
  "Bonds",
  "Dis. Res", 
  "Impropers",
  "LJ",
  "LJ-14",
  "Coulomb (LR)",
  "Proper Dih.", 
  "Position Rest.",
  "Ryckaert-Bell.",
  "Shake",
  "Coulomb (SR)",
  "TA-Disres",
  "Potential",
  "Kinetic En.",
  "Temperature",
  "Total Energy"
  };
#define F_NRE asize(EType)
#define NBVP 27

int main(int argc,char *argv[])
{
  FILE     *in,*out;
  t_energy eold[F_NRE+NBVP];
  int      teller=0;
  real     t;
  t_filenm fnm[] = {
    { efTRJ, "-f", NULL, ffREAD  },
    { efENE, "-o", NULL, ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,NFILE,fnm,TRUE,NULL);
  
  in=ftp2FILE(efTRJ,NFILE,fnm,"r");
  out=ftp2FILE(efENE,NFILE,fnm,"w");
  wr_ener_nms(out,F_NRE,EType);
  while (next_e(in,&t,eold)) {
    if ((teller++ % 1) == 0)
      fprintf(stderr,"\rFrame: %d",teller-1);
    wr_ener(out,t,teller,F_NRE,eold,NULL);
  }
  fprintf(stderr,"\n");
  fclose(in);
  fclose(out);
  
  thanx(stdout);
  
  return 0;
}
