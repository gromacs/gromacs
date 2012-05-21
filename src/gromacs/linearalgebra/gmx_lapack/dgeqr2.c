#include "../gmx_lapack.h"


void
F77_FUNC(dgeqr2,DGEQR2)(int *m,
	int *n,
	double *a,
	int *lda,
	double *tau,
	double *work,
	int *info)
{
  int k = (*m < *n) ? *m : *n;
  int i,i1,i2,i3;
  double aii;

  *info = 0;
  
  for(i=0;i<k;i++) {
    i1 = *m - i;
    i2 = ( (i+1) < (*m-1) ) ? (i+1) : (*m-1);
    i3 = 1;
    F77_FUNC(dlarfg,DLARFG)(&i1,&(a[i*(*lda)+i]),&(a[i*(*lda)+i2]),&i3,&(tau[i]));
    if(i<(*n-1)) {
      aii = a[i*(*lda)+i];
      a[i*(*lda)+i] = 1.0;
      i2 = *n - i - 1;
      F77_FUNC(dlarf,DLARF)("L",&i1,&i2,&(a[i*(*lda)+i]),&i3,&(tau[i]),
	     &(a[(i+1)*(*lda)+i]),lda,work);
      a[i*(*lda)+i] = aii;
    }
  }
  return;
}
