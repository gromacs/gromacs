/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
static char *SRCID_tstsqrtc_c = "$Id$";
#include "vec.h"

int main(int argc,char *argv[])
{
  real x,y,z,diff,av;
  int  i;

  printf("%12s  %12s  %12s  %12s  %12s\n","X","invsqrt(X)","1/sqrt(X)","Abs. Diff.","Rel. Diff.");
  for(i=1; (i<1000); i++) {
    x = i*1.0;
    y = invsqrt(x);
    z = 1.0/sqrt(x);
    diff = y-z;
    av   = 0.5*(y+z);
    printf("%12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n",x,y,z,diff,diff/z);
  }
  return 0;
}
