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
 * Grunge ROck MAChoS
 */
#include "stdio.h"
#include "math.h"

#define eNRNB 23

void breakit(int nprocs)
{
  double n[eNRNB],a=0;
  int    i;
  
  for(i=0; (i<eNRNB); i++)
    n[i]=i*0.34729139483721;
    
  for(i=0; (i<eNRNB); i++)
    a+=n[i]/nprocs;
  
  printf("a=%e\n",a);
}

void works(double nprocs)
{
  double *n,a=0;
  int    i;
    
  n=(double *)malloc(sizeof(*n)*eNRNB);
  
  for(i=0; (i<eNRNB); i++)
    n[i]=i*0.34729139483721;
    
  for(i=0; (i<eNRNB); i++)
    n[i]/=nprocs;

  for(i=0; (i<eNRNB); i++)
    a+=n[i];
    
  printf("a=%e\n",a);
}

int main(int argc,char *argv[])
{
  breakit(4);
  works(4);
  
  return 0;
}
