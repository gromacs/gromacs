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
#include <math.h>

#include <types/simple.h>
#include "gmx_blas.h"

void
F77_FUNC(sgemm,SGEMM)(const char *transa,
       const char *transb,
       int *m__,
       int *n__,
       int *k__,
       float *alpha__,
       float *a,
       int *lda__,
       float *b,
       int *ldb__,
       float *beta__,
       float *c,
       int *ldc__)
{
  const char tra=toupper(*transa);
  const char trb=toupper(*transb);
  float temp;
  int i,j,l;
  int nrowa,ncola,nrowb;

  int m   = *m__;
  int n   = *n__;
  int k   = *k__;
  int lda = *lda__;
  int ldb = *ldb__;
  int ldc = *ldc__;
  
  float alpha = *alpha__;
  float beta  = *beta__;
  
  if(tra=='N') {
    nrowa = m;
    ncola = k;
  } else {
    nrowa = k;
    ncola = m;
  }

  if(trb=='N') 
    nrowb = k;
   else 
    nrowb = n;
  
  if(m==0 || n==0 || (( fabs(alpha)<GMX_FLOAT_MIN || k==0) && fabs(beta-1.0)<GMX_FLOAT_EPS))
    return;

  if(fabs(alpha)<GMX_FLOAT_MIN) {
    if(fabs(beta)<GMX_FLOAT_MIN) {
      for(j=0;j<n;j++)
	for(i=0;i<m;i++)
	  c[j*(ldc)+i] = 0.0;
    } else {
      /* nonzero beta */
      for(j=0;j<n;j++)
	for(i=0;i<m;i++)
	  c[j*(ldc)+i] *= beta;
    }
    return;
  }

  if(trb=='N') {
    if(tra=='N') {
      
      for(j=0;j<n;j++) {
	if(fabs(beta)<GMX_FLOAT_MIN) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] = 0.0;
	} else if(fabs(beta-1.0)>GMX_FLOAT_EPS) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] *= beta;
	} 
	for(l=0;l<k;l++) {
	  if( fabs(b[ j*(ldb) + l ])>GMX_FLOAT_MIN) {
	    temp = alpha * b[ j*(ldb) + l ];
	    for(i=0;i<m;i++)
	      c[j*(ldc)+i] += temp * a[l*(lda)+i]; 
	  }
	}
      }
    } else {
      /* transpose A, but not B */
      for(j=0;j<n;j++) {
	for(i=0;i<m;i++) {
	  temp = 0.0;
	  for(l=0;l<k;l++) 
	    temp += a[i*(lda)+l] * b[j*(ldb)+l];
	  if(fabs(beta)<GMX_FLOAT_MIN)
	    c[j*(ldc)+i] = alpha * temp;
	  else
	    c[j*(ldc)+i] = alpha * temp + beta * c[j*(ldc)+i];
	}
      }
    }
  } else {
    /* transpose B */
    if(tra=='N') {

      /* transpose B, but not A */

      for(j=0;j<n;j++) {
	if(fabs(beta)<GMX_FLOAT_MIN) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] = 0.0;
	} else if(fabs(beta-1.0)>GMX_FLOAT_EPS) {
	  for(i=0;i<m;i++)
	    c[j*(ldc)+i] *= beta;
	} 
	for(l=0;l<k;l++) {
	  if( fabs(b[ l*(ldb) + j ])>GMX_FLOAT_MIN) {
	    temp = alpha * b[ l*(ldb) + j ];
	    for(i=0;i<m;i++)
	      c[j*(ldc)+i] += temp * a[l*(lda)+i]; 
	  }
	}
      }
 
    } else {
      /* Transpose both A and B */
       for(j=0;j<n;j++) {
	for(i=0;i<m;i++) {
	  temp = 0.0;
	  for(l=0;l<k;l++) 
	    temp += a[i*(lda)+l] * b[l*(ldb)+j];
	  if(fabs(beta)<GMX_FLOAT_MIN)
	    c[j*(ldc)+i] = alpha * temp;
	  else
	    c[j*(ldc)+i] = alpha * temp + beta * c[j*(ldc)+i];
	}
       }
    }
  }
}
