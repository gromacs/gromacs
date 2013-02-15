/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include <math.h>
#include <ctype.h>
#include "gmx_lapack.h"


float
F77_FUNC(slanst,SLANST)(const char *norm,
	int *n,
	float *d,
	float *e)
{
  const char ch=toupper(*norm);
  float dtemp,max,val,scale,sum;
  int i,j;


  if(*n<=0)
    return 0.0;
  
  switch(ch) {
  case 'M':
    max = fabs(d[*n-1]);
      for(i=0;i<(*n-1);i++) {
	dtemp = fabs(d[i]);
	if(dtemp>max)
	  max = dtemp;
	dtemp = fabs(e[i]);
	if(dtemp>max)
	  max = dtemp;
      }
    val = max;
    break;
    
  case 'O':
  case '1':
  case 'I':

    if(*n==1)
      val = fabs(d[0]);
    else {
      max = fabs(d[0]) + fabs(e[0]);
      dtemp = fabs(e[*n-2]) + fabs(d[*n-1]);
      if(dtemp>max)
	max = dtemp;
      for(i=1;i<(*n-1);i++) {
	dtemp = fabs(d[i]) + fabs(e[i]) + fabs(e[i-1]);
	if(dtemp>max)
	  max = dtemp;
      }
      val = max;
    }
    break;

  case 'F':
  case 'E':
    scale = 0.0;
    sum   = 1.0;
    i = *n-1;
    j = 1;
    if(*n>1) {
      F77_FUNC(slassq,SLASSQ)(&i,e,&j,&scale,&sum);
      sum *= 2;
    }
    F77_FUNC(slassq,SLASSQ)(n,d,&j,&scale,&sum);
    val = scale * sqrt(sum);
    break;
    
  default:
    val = 0.0;
    break;
  }
  return val;
}
