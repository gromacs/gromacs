#include "gmx_lapack.h"

/* LAPACK */
void
F77_FUNC(dlaswp,DLASWP)(int *n,
	double *a,
	int *lda,
	int *k1,
	int *k2,
	int *ipiv,
	int *incx)
{
  int ix0,i1,i2,inc,n32;
  int ix,i,j,ip,k;
  double temp;

  if(*incx>0) {
    ix0 = *k1 - 1;
    i1 = *k1 - 1;
    i2 = *k2;
    inc = 1;
  } else if(*incx<0) {
    ix0 = *incx * (1- *k2);
    i1 = *k2 - 1;
    i2 = *k1;
    inc = -1;
  } else
    return;

  n32 = *n / 32;
  
  n32 *= 32;


  if(n32!=0) {
    for(j=0;j<n32;j+=32) {
      ix = ix0;
      for(i=i1;i<i2;i+=inc,ix+=*incx) {
	ip = ipiv[ix] - 1;
	if(ip != i) {
	  for(k=j;k<j+32;k++) {
	    temp = a[(k)*(*lda)+i];
	    a[(k)*(*lda)+i] = a[(k)*(*lda)+ip];
	    a[(k)*(*lda)+ip] = temp;
	  }
	}
      }
    }
  }
  if(n32!=*n) {
    ix = ix0;
    for(i=i1;i<i2;i+=inc,ix+=*incx) {
      ip = ipiv[ix] - 1;
      if(ip != i) {
	for(k=n32;k<*n;k++) {
	    temp = a[(k)*(*lda)+i];
	    a[(k)*(*lda)+i] = a[(k)*(*lda)+ip];
	    a[(k)*(*lda)+ip] = temp;
	}
      }
    }
  }
  return;
}
