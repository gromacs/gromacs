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
#include "ns_mgmt.h"
#include "smalloc.h"
#include "wnblist.h"
#include "futil.h"


int main(int argc,char *argv[])
{
  FILE *in;
  int  i,j,natoms,mod;
  int **matrix;
  
  fprintf(stderr,"Natoms ? \n");
  scanf("%d",&natoms);
  snew(matrix,natoms);
  for(i=0; (i<natoms); i++)
    snew(matrix[i],natoms);
  for(i=1; (i<argc); i++) {
    in=ffopen(argv[i],"r");
    fprintf(stderr,"Reading %s\n",argv[i]);
#ifndef PRINT_SHIFT
    read_nblist(in,matrix);
#else
    read_nblistshift(in,matrix,natoms);
#endif
    fclose(in);
  }
  mod=1;
  for(i=0; (i<natoms); i+=mod) {
    printf("%3d-",i/mod);
    for(j=0; (j<=i); j+=mod)
      printf(" ");
    for(j=i+mod; (j<natoms); j+=mod) {
      printf("%1d",matrix[i][j]);
    }
    printf("-%3d\n",i/mod);
  }
  fprintf(stderr,"finshed\n");
  
  return 0;
}



