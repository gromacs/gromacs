#include "gmx_blas.h"
#include "gmx_lapack.h"

void
F77_FUNC(dorgl2,DORGL2)(int *m,
	int *n,
	int *k,
	double *a,
	int *lda,
	double *tau,
	double *work,
	int *info)
{
  int i,j,l;
  int i1,i2,i3;
  double d1;

  if(*m<=0)
    return;

  *info = 0;

  if(*k < *m) {
    for(j=0;j<*n;j++) {
      for(l=*k;l<*m;l++)
	a[j*(*lda)+l] = 0.0;
      if(j>=*k && j<*m)
	a[j*(*lda)+j] = 1.0;
    }
  }

  for(i=*k-1;i>=0;i--) {
    if(i<(*n-1)) {
      if(i<(*m-1)) {
	a[i*(*lda)+i] = 1.0;
	i1 = *m - i - 1;
	i2 = *n - i ;
	i3 = 1;
	F77_FUNC(dlarf,DLARF)("R",&i1,&i2,&(a[i*(*lda)+i]),lda,&(tau[i]),
	       &(a[i*(*lda)+i+1]),lda,work);
      }
      i1 = *n - i - 1;
      d1 = -tau[i];
      F77_FUNC(dscal,DSCAL)(&i1,&d1,&(a[(i+1)*(*lda)+i]),lda);
    }
    a[i*(*lda)+i] = 1.0 - tau[i];

    for(l=0;l<(i-1);l++)
      a[l*(*lda)+i] = 0.0;
  }
  return;
}
