#include<ctype.h>
#include "gmx_blas.h"

void
F77_FUNC(dsyr2k,DSYR2K)(char *uplo, 
                        char *trans,
                        int *n,
                        int *k,
                        double *alpha,
                        double *a,
                        int *lda,
                        double *b,
                        int *ldb,
                        double *beta,
                        double *c,
                        int *ldc)
{
  char ch1,ch2;
  int nrowa;
  int i,j,l;
  double temp1,temp2;

  ch1 = toupper(*uplo);
  ch2 = toupper(*trans);

  if(ch2 == 'N')
    nrowa = *n;
  else
    nrowa = *k;

  if(*n==0 || ( ( alpha==0 || *k==0 ) && *beta==1))
    return;

  if(*alpha == 0 ) {
    if(ch1=='U') {
      if(*beta==0) 
	for(j=1;j<=*n;j++) 
	  for(i=1;i<=j;i++)
	    c[(j-1)*(*ldc)+(i-1)] = 0.0;
      else
	for(j=1;j<=*n;j++) 
	  for(i=1;i<=j;i++)
	    c[(j-1)*(*ldc)+(i-1)] *= *beta;
    } else {
      /* lower */
      if(*beta==0) 
	for(j=1;j<=*n;j++) 
	  for(i=j;i<=*n;i++)
	    c[(j-1)*(*ldc)+(i-1)] = 0.0;
      else
	for(j=1;j<=*n;j++) 
	  for(i=j;i<=*n;i++)
	    c[(j-1)*(*ldc)+(i-1)] *= *beta;
    }
    return;
  }

  if(ch2=='N') {
    if(ch1=='U') {
      for(j=1;j<=*n;j++) {
	if(*beta==0)
	  for(i=1;i<=j;i++)
	     c[(j-1)*(*ldc)+(i-1)] = 0.0;
	else if(*beta != 1.0)
	  for(i=1;i<=j;i++)
	    c[(j-1)*(*ldc)+(i-1)] *= *beta;
	for(l=1;l<=*k;l++) {
	  if( a[(l-1)*(*lda)+(j-1)] != 0 ||
	      b[(l-1)*(*ldb)+(j-1)] != 0) {
	    temp1 = *alpha * b[(l-1)*(*ldb)+(j-1)];
	    temp2 = *alpha * a[(l-1)*(*lda)+(j-1)];
	    for(i=1;i<=j;i++)
	      c[(j-1)*(*ldc)+(i-1)] += 
		a[(l-1)*(*lda)+(i-1)] * temp1 + 
		b[(l-1)*(*ldb)+(i-1)] * temp2;
	  }
	}
      }
    } else {
      /* lower */
      for(j=1;j<=*n;j++) {
	if(*beta==0)
	  for(i=j;i<=*n;i++)
	    c[(j-1)*(*ldc)+(i-1)] = 0.0;
	else if(*beta != 1.0)
	  for(i=j;i<=*n;i++)
	    c[(j-1)*(*ldc)+(i-1)] *= *beta;
	for(l=1;l<=*k;l++) {
	  if( a[(l-1)*(*lda)+(j-1)] != 0 ||
	      b[(l-1)*(*ldb)+(j-1)] != 0) {
	    temp1 = *alpha * b[(l-1)*(*ldb)+(j-1)];
	    temp2 = *alpha * a[(l-1)*(*lda)+(j-1)];
	    for(i=j;i<=*n;i++)
	      c[(j-1)*(*ldc)+(i-1)] += 
		a[(l-1)*(*lda)+(i-1)] * temp1 + 
		b[(l-1)*(*ldb)+(i-1)] * temp2;
	  }
	}
      }
    }
  } else {
    /* transpose */
    if(ch1=='U') {
      for(j=1;j<=*n;j++) 
	for(i=1;i<=j;i++) {
	  temp1 = 0.0;
	  temp2 = 0.0;
	  for (l=1;l<=*k;l++) {
	     temp1 += a[(i-1)*(*lda)+(l-1)] * b[(j-1)*(*ldb)+(l-1)];
	     temp2 += b[(i-1)*(*ldb)+(l-1)] * a[(j-1)*(*lda)+(l-1)];
	  }
	  if(*beta==0)
	    c[(j-1)*(*ldc)+(i-1)] = *alpha * (temp1 + temp2);
	  else
	    c[(j-1)*(*ldc)+(i-1)] = *beta * c[(j-1)*(*ldc)+(i-1)] +
	      *alpha * (temp1 + temp2);
	}
    } else {
      /* lower */
      for(j=1;j<=*n;j++) 
	for(i=j;i<=*n;i++) {
	  temp1 = 0.0;
	  temp2 = 0.0;
	  for (l=1;l<=*k;l++) {
	     temp1 += a[(i-1)*(*lda)+(l-1)] * b[(j-1)*(*ldb)+(l-1)];
	     temp2 += b[(i-1)*(*ldb)+(l-1)] * a[(j-1)*(*lda)+(l-1)];
	  }
	  if(*beta==0)
	    c[(j-1)*(*ldc)+(i-1)] = *alpha * (temp1 + temp2);
	  else
	    c[(j-1)*(*ldc)+(i-1)] = *beta * c[(j-1)*(*ldc)+(i-1)] +
	      *alpha * (temp1 + temp2);
	}
    }
  }
  return;
}
