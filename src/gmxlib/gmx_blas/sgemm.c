#include <ctype.h>
#include <math.h>

#include <types/simple.h>
#include "gmx_blas.h"

void
F77_FUNC(sgemm,SGEMM)(char *transa,
       char *transb,
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
  char tra=toupper(*transa);
  char trb=toupper(*transb);
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
