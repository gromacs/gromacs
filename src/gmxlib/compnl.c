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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_compnl_c = "$Id$";

#include "ns.h"
#include "smalloc.h"
#include "wnblist.h"
#include "futil.h"
#include "macros.h"
#include "statutil.h"
#include "copyrite.h"

int main(int argc,char *argv[])
{
  FILE *in,*out;
  int  i,j,nmiss,mod;
  int ***mat;
  t_filenm fnm[] = {
    { efLOG, "-f1", NULL, ffREAD },
    { efLOG, "-f2", NULL, ffREAD },
    { efOUT, "-o",  "compnl", ffWRITE }
  };
#define NFILE asize(fnm)
  static int natoms=648;
  static t_pargs pa[] = {
    { "-nat",      FALSE, etINT, &natoms, "Number of atoms" }
  };

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,TRUE,NFILE,fnm,asize(pa),pa,0,NULL,0,NULL);

  snew(mat,2);  
  snew(mat[0],natoms);
  snew(mat[1],natoms);
  for(i=0; (i<natoms); i++) {
    snew(mat[0][i],natoms);
    snew(mat[1][i],natoms);
  }
  out = ffopen(ftp2fn(efOUT,NFILE,fnm),"w");
  
  for(i=0; (i<2); i++) {
    in = ffopen(fnm[i].fn,"r");
    fprintf(stderr,"Reading %s\n",fnm[i].fn);
    fprintf(out,   "Reading %s\n",fnm[i].fn);
    read_nblist(in,out,mat[i],natoms);
    fclose(in);
  }
  mod=1;
  fprintf(stderr,"Comparing Interaction Matrices\n");
  fprintf(out,   "Comparing Interaction Matrices\n");
  nmiss = 0;
  for(i=0; (i<natoms); i+=mod) {
    for(j=0; (j<natoms); j+=mod) {
      if (mat[0][i][j] != mat[1][i][j]) {
	fprintf(out,"i: %5d, j: %5d, shift1: %5d, shift2: %5d\n",
		i,j,mat[0][i][j]-1,mat[1][i][j]-1);
	nmiss++;
      }
    }
  }
  fprintf(out,"There were %d mismatches\n",nmiss);
  fprintf(out,"Done.\n");
  fclose(out);
  fprintf(stderr,"There were %d mismatches\n",nmiss);
  fprintf(stderr,"Finished\n");
}



