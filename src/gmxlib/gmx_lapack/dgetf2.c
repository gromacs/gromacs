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
#include <types/simple.h>

#include "gmx_blas.h"
#include "gmx_lapack.h"


void
F77_FUNC(dgetf2,DGETF2)(int *m,
	int *n,
	double *a,
	int *lda,
	int *ipiv,
	int *info)
{
  int j,jp,k,t1,t2,t3;
  double one,minusone;
  double tmp;

  one = 1.0;
  minusone = -1.0;

  if(*m<=0 || *n<=0)
    return;

  k = (*m < *n) ? *m : *n;
  for(j=1;j<=k;j++) {
    t1 = *m-j+1;
    t2 = 1;
    jp = j - 1 + F77_FUNC(idamax,IDAMAX)(&t1,&(a[(j-1)*(*lda)+(j-1)]),&t2);
    ipiv[j-1] = jp;
    if( fabs(a[(j-1)*(*lda)+(jp-1)])>GMX_DOUBLE_MIN ) {
      if(jp != j)
	F77_FUNC(dswap,DSWAP)(n,&(a[ j-1 ]),lda,&(a[ jp-1 ]),lda);
      
      if(j<*m) {
	t1 = *m-j;
	t2 = 1;
	tmp = 1.0/a[(j-1)*(*lda)+(j-1)];
	F77_FUNC(dscal,DSCAL)(&t1,&tmp,&(a[(j-1)*(*lda)+(j)]),&t2);
      }
    } else {
      *info = j;
    }

    if(j<k) {
      t1 = *m-j;
      t2 = *n-j;
      t3 = 1;
      F77_FUNC(dger,DGER)(&t1,&t2,&minusone,&(a[(j-1)*(*lda)+(j)]),&t3,
	    &(a[(j)*(*lda)+(j-1)]),lda, &(a[(j)*(*lda)+(j)]),lda);
    }
  }
  return;
}
