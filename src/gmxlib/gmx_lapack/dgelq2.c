#include "gmx_lapack.h"


void
F77_FUNC(dgelq2,DGELQ2)(int *m,
	int *n,
	double *a,
	int *lda,
	double *tau,
	double *work,
	int *info)
{
  int k = (*m < *n) ? *m : *n;
  int i,i1,i2;
  double aii;

  for(i=0;i<k;i++) {
    i1 = *n-i;
    i2 = ( i < (*n-1)) ? i : (*n-1);
    F77_FUNC(dlarfg,DLARFG)(&i1,&(a[i*(*lda)+i]),&(a[i2*(*lda)+i]),lda,&(tau[i]));
    
    if(i<(*m-1)) {
      aii = a[i*(*lda)+i];
      a[i*(*lda)+i] = 1.0;
      i1 = *m - i - 1;
      i2 = *n - i;
      F77_FUNC(dlarf,DLARF)("R",&i1,&i2,&(a[i*(*lda)+i]),lda,&(tau[i]),&(a[i*(*lda)+i+1]),lda,work);
      a[i*(*lda)+i] = aii;
    }
  }
  *info = 0;
  return;
}
