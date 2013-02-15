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
#include<math.h>
#include<ctype.h>
#include "gmx_lapack.h"


double
F77_FUNC(dlange,DLANGE)(const char *norm,
	int *m,
	int *n,
	double *a,
	int *lda,
	double *work)
{
  const char ch=toupper(*norm);
  double dtemp,sum,max,val,scale;
  int i,j;

  switch(ch) {
  case 'M':
    max = 0.0;
    for(j=0;j<*n;j++)
      for(i=0;i<*m;i++) {
	dtemp = fabs(a[j*(*lda)+i]);
	if(dtemp>max)
	  max = dtemp;
      }
    val = max;
    break;

  case 'O':
  case '1':
    max = 0.0;
    for(j=0;j<*n;j++) {
      sum = 0.0;
      for(i=0;i<*m;i++) 
	sum += fabs(a[j*(*lda)+i]);
      if(sum>max)
	max = sum;
    }
    val = max;
    break;

  case 'I':
    for(i=0;i<*m;i++)
      work[i] = 0.0;
    for(j=0;j<*n;j++)
      for(i=0;i<*m;i++)
	work[i] += fabs(a[j*(*lda)+i]);
    max = 0;
    for(i=0;i<*m;i++)
      if(work[i]>max)
	max=work[i];
    val = max;
    break;

  case 'F':
  case 'E':
    scale = 0.0;
    sum   = 1.0;
    i = 1;
    for(j=0;j<*n;j++) 
      F77_FUNC(dlassq,DLASSQ)(m,&(a[j*(*lda)+0]),&i,&scale,&sum);
    val = scale*sqrt(sum);
    break;

  default:
    val = 0.0;
    break;
  }
  return val;
}
