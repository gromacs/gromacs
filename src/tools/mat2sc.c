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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
#include <stdio.h>

int main(int argc, char *argv[]) 
{
  double y;
  int    i,j,k,ncol=1;
 
  if (argc < 2) {
    fprintf(stderr,"Usage: %s ncol\n",argv[0]);
    exit(1);
  }
    
  ncol=atoi(argv[1]);
  
  i=0;
  printf("# Hacked by gromacs\n# Don't edit\n\n");
  printf("set numeric\n");
  do {
    for(j=0; (j<ncol); j++) {
      if (scanf("%lf",&y) != 1)
	return 0;
      if (j < 26)
	printf("let %c%d = %g\n",'A'+j,i,y);
      else {
	printf("let %c%c%d = %g\n",'A'-1+j/26,'A'+j%26,i,y);
      }
    }
    i++;
  } while (1);
  
  return 0;
}


