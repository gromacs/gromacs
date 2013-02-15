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
#include "gmx_lapack.h"
#include "lapack_limits.h"

void
F77_FUNC(dorgbr,DORGBR)(const char *vect,
	int *m,
	int *n,
	int *k,
	double *a,
	int *lda,
	double *tau,
	double *work,
	int *lwork,
	int *info)
{
  int wantq,iinfo,j,i,i1,wrksz;
  int mn = (*m < *n) ? *m : *n;

  wantq = (*vect=='Q' || *vect=='q');

  *info = 0;
  wrksz = mn*DORGBR_BLOCKSIZE;
  if(*lwork==-1) {
    work[0] = wrksz;
    return;
  }
  
  if(*m==0 || *n==0)
    return;

  if(wantq) {
    if(*m>=*k)
      F77_FUNC(dorgqr,DORGQR)(m,n,k,a,lda,tau,work,lwork,&iinfo);
    else {
      for(j=*m;j>=2;j--) {
	a[(j-1)*(*lda)+0] = 0.0;
	for(i=j+1;i<=*m;i++)
	  a[(j-1)*(*lda)+(i-1)] = a[(j-2)*(*lda)+(i-1)]; 
      }
      a[0] = 1.0;
      for(i=2;i<=*m;i++)
	a[i-1] = 0.0;
      if(*m>1) {
	i1 = *m-1;
	F77_FUNC(dorgqr,DORGQR)(&i1,&i1,&i1,&(a[*lda+1]),lda,tau,work,lwork,&iinfo);
      }
    }
  } else {
    if(*k<*n)
      F77_FUNC(dorglq,DORGLQ)(m,n,k,a,lda,tau,work,lwork,&iinfo);
    else {
      a[0] = 1.0;
      for(i=2;i<=*m;i++)
	a[i-1] = 0.0;
      for(j=2;j<=*n;j++) {
	for(i=j-1;i>=2;i--)
	  a[(j-1)*(*lda)+(i-1)] = a[(j-1)*(*lda)+(i-2)]; 
	a[(j-1)*(*lda)+0] = 0.0;
      }
      if(*n>1) {
	i1 = *n-1;
	F77_FUNC(dorglq,DORGLQ)(&i1,&i1,&i1,&(a[*lda+1]),lda,tau,work,lwork,&iinfo);
      }
    }
  }
  work[0] = wrksz;
  return;
}
 
