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


#include <types/simple.h>
#include "gmx_blas.h"

void
F77_FUNC(dsyr2k,DSYR2K)(const char *uplo, 
	const char *trans,
	int *n__,
	int *k__,
	double *alpha__,
	double *a,
	int *lda__,
	double *b,
	int *ldb__,
	double *beta__,
	double *c,
	int *ldc__)
{
  char ch1,ch2;
  int nrowa;
  int i,j,l;
  double temp1,temp2;

  
  int n = *n__;
  int k = *k__;
  int lda = *lda__;
  int ldb = *ldb__;
  int ldc = *ldc__;
  
  double alpha = *alpha__;
  double beta  = *beta__;
  
  ch1 = toupper(*uplo);
  ch2 = toupper(*trans);

  if(ch2 == 'N')
    nrowa = n;
  else
    nrowa = k;

  if(n==0 || ( ( fabs(alpha)<GMX_DOUBLE_MIN || k==0 ) && fabs(beta-1.0)<GMX_DOUBLE_EPS))
    return;

  if(fabs(alpha)<GMX_DOUBLE_MIN ) {
    if(ch1=='U') {
      if(fabs(beta)<GMX_DOUBLE_MIN) 
	for(j=1;j<=n;j++) 
	  for(i=1;i<=j;i++)
	    c[(j-1)*(ldc)+(i-1)] = 0.0;
      else
	for(j=1;j<=n;j++) 
	  for(i=1;i<=j;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
    } else {
      /* lower */
      if(fabs(beta)<GMX_DOUBLE_MIN) 
	for(j=1;j<=n;j++) 
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] = 0.0;
      else
	for(j=1;j<=n;j++) 
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
    }
    return;
  }

  if(ch2=='N') {
    if(ch1=='U') {
      for(j=1;j<=n;j++) {
	if(fabs(beta)<GMX_DOUBLE_MIN)
	  for(i=1;i<=j;i++)
	     c[(j-1)*(ldc)+(i-1)] = 0.0;
	else if(fabs(beta-1.0)>GMX_DOUBLE_EPS)
	  for(i=1;i<=j;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
	for(l=1;l<=k;l++) {
	  if( fabs(a[(l-1)*(lda)+(j-1)])>GMX_DOUBLE_MIN ||
	      fabs(b[(l-1)*(ldb)+(j-1)])>GMX_DOUBLE_MIN) {
	    temp1 = alpha * b[(l-1)*(ldb)+(j-1)];
	    temp2 = alpha * a[(l-1)*(lda)+(j-1)];
	    for(i=1;i<=j;i++)
	      c[(j-1)*(ldc)+(i-1)] += 
		a[(l-1)*(lda)+(i-1)] * temp1 + 
		b[(l-1)*(ldb)+(i-1)] * temp2;
	  }
	}
      }
    } else {
      /* lower */
      for(j=1;j<=n;j++) {
	if(fabs(beta)<GMX_DOUBLE_MIN)
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] = 0.0;
	else if(fabs(beta-1.0)>GMX_DOUBLE_EPS)
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
	for(l=1;l<=k;l++) {
	  if( fabs(a[(l-1)*(lda)+(j-1)])>GMX_DOUBLE_MIN ||
	      fabs(b[(l-1)*(ldb)+(j-1)])>GMX_DOUBLE_MIN) {
	    temp1 = alpha * b[(l-1)*(ldb)+(j-1)];
	    temp2 = alpha * a[(l-1)*(lda)+(j-1)];
	    for(i=j;i<=n;i++)
	      c[(j-1)*(ldc)+(i-1)] += 
		a[(l-1)*(lda)+(i-1)] * temp1 + 
		b[(l-1)*(ldb)+(i-1)] * temp2;
	  }
	}
      }
    }
  } else {
    /* transpose */
    if(ch1=='U') {
      for(j=1;j<=n;j++) 
	for(i=1;i<=j;i++) {
	  temp1 = 0.0;
	  temp2 = 0.0;
	  for (l=1;l<=k;l++) {
	     temp1 += a[(i-1)*(lda)+(l-1)] * b[(j-1)*(ldb)+(l-1)];
	     temp2 += b[(i-1)*(ldb)+(l-1)] * a[(j-1)*(lda)+(l-1)];
	  }
	  if(fabs(beta)<GMX_DOUBLE_MIN)
	    c[(j-1)*(ldc)+(i-1)] = alpha * (temp1 + temp2);
	  else
	    c[(j-1)*(ldc)+(i-1)] = beta * c[(j-1)*(ldc)+(i-1)] +
	      alpha * (temp1 + temp2);
	}
    } else {
      /* lower */
      for(j=1;j<=n;j++) 
	for(i=j;i<=n;i++) {
	  temp1 = 0.0;
	  temp2 = 0.0;
	  for (l=1;l<=k;l++) {
	     temp1 += a[(i-1)*(lda)+(l-1)] * b[(j-1)*(ldb)+(l-1)];
	     temp2 += b[(i-1)*(ldb)+(l-1)] * a[(j-1)*(lda)+(l-1)];
	  }
	  if(fabs(beta)<GMX_DOUBLE_MIN)
	    c[(j-1)*(ldc)+(i-1)] = alpha * (temp1 + temp2);
	  else
	    c[(j-1)*(ldc)+(i-1)] = beta * c[(j-1)*(ldc)+(i-1)] +
	      alpha * (temp1 + temp2);
	}
    }
  }
  return;
}
