#include "gmx_blas.h"
#include "gmx_lapack.h"

void
F77_FUNC(sorg2r,SORG2R)(int *m,
	int *n,
	int *k,
	float *a,
	int *lda,
	float *tau,
	float *work,
	int *info)
{
  int i,j,l;
  int i1,i2,i3;
  float d1;

  if(*n<=0)
    return;

  *info = 0;

  for(j=*k;j<*n;j++) {
    for(l=0;l<*m;l++)
      a[j*(*lda)+l] = 0.0;
    a[j*(*lda)+j] = 1.0;
  }

  for(i=*k-1;i>=0;i--) {
    if(i<(*n-1)) {
      a[i*(*lda)+i] = 1.0;
      i1 = *m - i;
      i2 = *n - i - 1;
      i3 = 1;
      F77_FUNC(slarf,SLARF)("L",&i1,&i2,&(a[i*(*lda)+i]),&i3,&(tau[i]),
	     &(a[(i+1)*(*lda)+i]),lda,work);
    }
    if(i<(*m-1)) {
      i1 = *m - i - 1;
      d1 = -tau[i];
      i2 = 1;
      F77_FUNC(sscal,SSCAL)(&i1,&d1,&(a[(i)*(*lda)+i+1]),&i2);
      a[i*(*lda)+i] = 1.0 - tau[i];

      for(l=0;l<(i-1);l++)
	a[i*(*lda)+l] = 0.0;
    }
  }
  return;
}
