#include <ctype.h>
#include <math.h>

#include <types/simple.h>
#include "gmx_blas.h"

void
F77_FUNC(ssyr2k,SSYR2K)(char *uplo, 
                        char *trans,
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
  char ch1,ch2;
  int nrowa;
  int i,j,l;
  float temp1,temp2;

  int n = *n__;
  int k = *k__;
  int lda = *lda__;
  int ldb = *ldb__;
  int ldc = *ldc__;
  
  float alpha = *alpha__;
  float beta  = *beta__;
  
  ch1 = toupper(*uplo);
  ch2 = toupper(*trans);

  if(ch2 == 'N')
    nrowa = n;
  else
    nrowa = k;

  if(n==0 || ( ( fabs(alpha)<GMX_FLOAT_MIN || k==0 ) && fabs(beta-1.0)<GMX_FLOAT_EPS))
    return;

  if(fabs(alpha)<GMX_FLOAT_MIN) {
    if(ch1=='U') {
      if(fabs(beta)<GMX_FLOAT_MIN) 
	for(j=1;j<=n;j++) 
	  for(i=1;i<=j;i++)
	    c[(j-1)*(ldc)+(i-1)] = 0.0;
      else
	for(j=1;j<=n;j++) 
	  for(i=1;i<=j;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
    } else {
      /* lower */
      if(fabs(beta)<GMX_FLOAT_MIN) 
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
	if(fabs(beta)<GMX_FLOAT_MIN)
	  for(i=1;i<=j;i++)
	     c[(j-1)*(ldc)+(i-1)] = 0.0;
	else if(fabs(beta-1.0)>GMX_FLOAT_EPS)
	  for(i=1;i<=j;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
	for(l=1;l<=k;l++) {
	  if( fabs(a[(l-1)*(lda)+(j-1)])>GMX_FLOAT_MIN ||
	      fabs(b[(l-1)*(ldb)+(j-1)])>GMX_FLOAT_MIN) {
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
	if(fabs(beta)<GMX_FLOAT_MIN)
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] = 0.0;
	else if(fabs(beta-1.0)>GMX_FLOAT_EPS)
	  for(i=j;i<=n;i++)
	    c[(j-1)*(ldc)+(i-1)] *= beta;
	for(l=1;l<=k;l++) {
	  if( fabs(a[(l-1)*(lda)+(j-1)])>GMX_FLOAT_MIN ||
	      fabs(b[(l-1)*(ldb)+(j-1)])>GMX_FLOAT_MIN) {
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
	  if(fabs(beta)<GMX_FLOAT_MIN)
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
	  if(fabs(beta)<GMX_FLOAT_MIN)
	    c[(j-1)*(ldc)+(i-1)] = alpha * (temp1 + temp2);
	  else
	    c[(j-1)*(ldc)+(i-1)] = beta * c[(j-1)*(ldc)+(i-1)] +
	      alpha * (temp1 + temp2);
	}
    }
  }
  return;
}
