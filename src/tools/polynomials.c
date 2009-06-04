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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "gmx_fatal.h"
#include "gstat.h"

real LegendreP(real x,unsigned long m)

{
  real polynomial=0,x2,x3;

  switch (m) {
  case eacP0:
    polynomial=1.0;
    break;
  case eacP1:
    polynomial=x;
    break;
  case eacP2:
    x2=x*x;
    polynomial=1.5*x2 - 0.5;
    break;
  case eacP3:
    x2=x*x;
    polynomial=(35*x2*x2 - 30*x2 + 3)/8;
    break;
  case eacP4:
    x2=x*x;
    x3=x2*x;
    polynomial=(63*x3*x2 - 70*x3 + 15*x)/8;
    break;
  default:
    gmx_fatal(FARGS,"Legendre polynomials of order %d are not supported, %s %d",
		m,__FILE__,__LINE__);
  }
  return (polynomial);
}
