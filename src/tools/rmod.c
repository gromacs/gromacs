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
 * GROtesk MACabre and Sinister
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef enum { FALSE, TRUE }  bool;

bool bRmod(double a,double b)
{
  int iq;
  double fac = 100000.0;
  
  iq = ((fac+0.5)*a)/(fac*b);
  
  if (fabs(a-b*iq) <= (a/fac))
    return TRUE;
  else
    return FALSE;
}

void main()
{
  int i,kk,nmod;
  float a,b,dt;
  
  dt   = 0.002;
  b    = 10.0;
  nmod = 50;
  for(i=0; (i<10000); i++) {
    a  = i*(nmod*dt);
    kk = bRmod(a,b);
    if (kk == 1)
      printf("%5d  %10f  %5d ***\n",i,a,kk);
    else
      printf("%5d  %10f  %5d\n",i,a,kk);
  }
}
