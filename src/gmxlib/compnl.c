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
 * Great Red Oystrich Makes All Chemists Sane
 */
#include "ns_mgmt.h"
#include "smalloc.h"
#include "wnblist.h"
#include "futil.h"

int main(int argc,char *argv[])
{
  FILE *in;
  int  i,j,natoms,mod;
  int **mat1,**mat2;
  
  if (argc < 3) {
    fprintf(stderr,"Usage: %s ref.log test1.log ...\n",argv[0]);
    exit(1);
  }
  fprintf(stderr,"Natoms ? \n");
  scanf("%d",&natoms);
  snew(mat1,natoms);
  snew(mat2,natoms);
  for(i=0; (i<natoms); i++) {
    snew(mat1[i],natoms);
    snew(mat2[i],natoms);
  }
  in=ffopen(argv[1],"r");
  fprintf(stderr,"Reading %s\n",argv[1]);
  read_nblistshift(in,mat1,natoms);
  fclose(in);
  
  for(i=2; (i<argc); i++) {
    in=ffopen(argv[i],"r");
    fprintf(stderr,"Reading %s\n",argv[i]);
    read_nblistshift(in,mat2,natoms);
    fclose(in);
  }
  mod=1;
  fprintf(stderr,"Comparing Interaction Matrices\n");
  printf("Comparing Interaction Matrices\n");
  for(i=0; (i<natoms); i+=mod) {
    for(j=0; (j<natoms); j+=mod) {
      if (mat1[i][j] != mat2[i][j])
	printf("i: %5d, j: %5d, Grp-SEQ: %5d, Grp-PAR: %5d\n",
	       i,j,mat1[i][j]-1,mat2[i][j]-1);
    }
  }
  fprintf(stderr,"Finished\n");
}



