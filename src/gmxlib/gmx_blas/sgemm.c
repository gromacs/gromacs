#include <ctype.h>
#include "gmx_blas.h"

void
F77_FUNC(sgemm,SGEMM)(char *transa,
                      char *transb,
                      int *m,
                      int *n,
                      int *k,
                      float *alpha,
                      float *a,
                      int *lda,
                      float *b,
                      int *ldb,
                      float *beta,
                      float *c,
                      int *ldc)
{
  char tra=toupper(*transa);
  char trb=toupper(*transb);
  float temp;
  int i,j,l;
  int nrowa,ncola,nrowb;

  if(tra=='N') {
    nrowa = *m;
    ncola = *k;
  } else {
    nrowa = *k;
    ncola = *m;
  }

  if(trb=='N') 
    nrowb = *k;
   else 
    nrowb = *n;
  
  if(*m==0 || *n==0 || (( *alpha==0.0 || *k==0) && *beta==1.0))
    return;

  if(*alpha==0.0) {
    if(*beta==0.0) {
      for(j=0;j<*n;j++)
	for(i=0;i<*m;i++)
	  c[j*(*ldc)+i] = 0.0;
    } else {
      /* nonzero beta */
      for(j=0;j<*n;j++)
	for(i=0;i<*m;i++)
	  c[j*(*ldc)+i] *= *beta;
    }
    return;
  }

  if(trb=='N') {
    if(tra=='N') {
      
      for(j=0;j<*n;j++) {
	if(*beta==0.0) {
	  for(i=0;i<*m;i++)
	    c[j*(*ldc)+i] = 0.0;
	} else if(*beta != 1.0) {
	  for(i=0;i<*m;i++)
	    c[j*(*ldc)+i] *= *beta;
	} 
	for(l=0;l<*k;l++) {
	  if( b[ j*(*ldb) + l ] != 0.0) {
	    temp = *alpha * b[ j*(*ldb) + l ];
	    for(i=0;i<*m;i++)
	      c[j*(*ldc)+i] += temp * a[l*(*lda)+i]; 
	  }
	}
      }
    } else {
      /* transpose A, but not B */
      for(j=0;j<*n;j++) {
	for(i=0;i<*m;i++) {
	  temp = 0.0;
	  for(l=0;l<*k;l++) 
	    temp += a[i*(*lda)+l] * b[j*(*ldb)+l];
	  if(*beta==0.0)
	    c[j*(*ldc)+i] = *alpha * temp;
	  else
	    c[j*(*ldc)+i] = *alpha * temp + *beta * c[j*(*ldc)+i];
	}
      }
    }
  } else {
    /* transpose B */
    if(tra=='N') {

      /* transpose B, but not A */

      for(j=0;j<*n;j++) {
	if(*beta==0.0) {
	  for(i=0;i<*m;i++)
	    c[j*(*ldc)+i] = 0.0;
	} else if(*beta != 1.0) {
	  for(i=0;i<*m;i++)
	    c[j*(*ldc)+i] *= *beta;
	} 
	for(l=0;l<*k;l++) {
	  if( b[ l*(*ldb) + j ] != 0.0) {
	    temp = *alpha * b[ l*(*ldb) + j ];
	    for(i=0;i<*m;i++)
	      c[j*(*ldc)+i] += temp * a[l*(*lda)+i]; 
	  }
	}
      }
 
    } else {
      /* Transpose both A and B */
       for(j=0;j<*n;j++) {
	for(i=0;i<*m;i++) {
	  temp = 0.0;
	  for(l=0;l<*k;l++) 
	    temp += a[i*(*lda)+l] * b[l*(*ldb)+j];
	  if(*beta==0.0)
	    c[j*(*ldc)+i] = *alpha * temp;
	  else
	    c[j*(*ldc)+i] = *alpha * temp + *beta * c[j*(*ldc)+i];
	}
       }
    }
  }
}
