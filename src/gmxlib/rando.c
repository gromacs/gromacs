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
static char *SRCID_rando_c = "$Id$";
#include <time.h>
#include <unistd.h>
#include "typedefs.h"

int make_seed(void)
{
  return (int)(((long)time(NULL)+(long)getpid()) % (long)1000000);
}
	
real rando(int *ig)
     /* generate a random number. */
{
  int  irand;

  int  m    = 100000000;
  real rm   = 100000000.0;  /* same number as m, but real format */
  int  m1   = 10000;
  int  mult = 31415821;
  
  real r;
  int  irandh,irandl,multh,multl;

  irand = abs(*ig) % m;
  
  /* multiply irand by mult, but take into account that overflow
   * must be discarded, and do not generate an error.
   */
  irandh = irand / m1;
  irandl = irand % m1;
  multh  = mult / m1;
  multl  = mult % m1;
  irand  = ((irandh*multl+irandl*multh) % m1) * m1 + irandl*multl;
  irand  = (irand + 1) % m;

  /* convert irand to a real random number between 0 and 1. */
  r = (irand / 10);
  r = r * 10 / rm;
  if ((r <= 0) || (r > 1))
    r = 0.0;
  *ig = irand;
  
  return r;
}        



