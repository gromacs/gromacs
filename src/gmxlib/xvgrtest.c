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
#include "xvgr.h"

int main(int argc,char *argv[])
{
  int i,j;
  int ij[5][5]= {
    1,2,2,1,3,
    1,2,2,1,2,
    1,2,1,1,2,
    1,2,2,1,2,
    1,2,2,1,3
  };
  
  xvgr_world(stdout,0,0,5,5);
  
  for(i=0; (i<5); i++)
    for(j=0; (j<5); j++)
      xvgr_box(stdout,
	       evWorld,
	       i,j,i+1,j+1,
	       elSolid,1,ij[i][j],
	       eppColor,ij[i][j],0);
}
