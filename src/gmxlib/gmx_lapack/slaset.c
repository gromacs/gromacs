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
#include <ctype.h>
#include "gmx_lapack.h"


void
F77_FUNC(slaset,SLASET)(const char *uplo,
	int *m,
	int *n,
	float *alpha,
	float *beta,
	float *a,
	int *lda)
{
  int i,j,k;
  const char ch=toupper(*uplo);

  if(ch=='U') {
    for(j=1;j<*n;j++) {
      k = (j < *m) ? j : *m;
      for(i=0;i<k;i++)
	a[j*(*lda)+i] = *alpha;
    }
  } else if(ch=='L') {
    k = (*m < *n) ? *m : *n;
    for(j=0;j<k;j++) {
      for(i=j+1;i<*m;i++)
	a[j*(*lda)+i] = *alpha;
    }
  } else {
    for(j=0;j<*n;j++) {
      for(i=0;i<*m;i++)
	a[j*(*lda)+i] = *alpha;
    }    
  }

  k = (*m < *n) ? *m : *n;
  for(i=0;i<k;i++)
    a[i*(*lda)+i] = *beta;
}
