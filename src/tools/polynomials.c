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
 * GRowing Old MAkes el Chrono Sweat
 */
#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "fatal.h"
#include "gstat.h"

real LegendreP(real x,unsigned long m)

{
  real polynomial,x2,x3;

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
    fatal_error(0,"Legendre polynomials of order %d are not supported, %s %d",
		m,__FILE__,__LINE__);
  }
  return (polynomial);
}
