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
static char *SRCID_tstsqrtc_c = "$Id$";

#include "vec.h"

int main(int argc,char *argv[])
{
  real x,y,z,diff,av;
  int  i;

  printf("%12s  %12s  %12s  %12s  %12s\n","X","Y","Z","Dy","Dy/Z");
  for(i=1; (i<1000); i++) {
    x = i*1.0;
    y = invsqrt(x);
    z = 1.0/sqrt(x);
    diff = y-z;
    av   = 0.5*(y+z);
    printf("%12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n",x,y,z,diff,diff/z);
  }
}
