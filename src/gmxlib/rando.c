/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <time.h>
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#include <process.h>
#endif
#include "sysstuff.h"
#include "typedefs.h"
#include "random.h"

int make_seed(void)
{
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
  return (int)_getpid();
#else
  return (int)getpid();
#endif
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



