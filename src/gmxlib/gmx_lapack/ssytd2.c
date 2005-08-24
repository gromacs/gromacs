#include <ctype.h>
#include <math.h>

#include <types/simple.h>

#include "gmx_blas.h"
#include "gmx_lapack.h"

void
F77_FUNC(ssytd2,SSYTD2)(char *    uplo,
	int *     n,
	float *  a,
	int *     lda,
	float *  d,
	float *  e,
	float *  tau,
	int *     info)
{
  float minusone,zero;
  float taui,alpha,tmp;
  int ti1,ti2,ti3,i;
  char ch=toupper(*uplo);

  zero = 0.0;
  minusone = -1.0;

  if(*n<=0)
    return;

  if(ch=='U') {
    for(i=*n-1;i>=1;i--) {

      ti1 = 1;
      F77_FUNC(slarfg,SLARFG)(&i,&(a[i*(*lda)+(i-1)]),&(a[i*(*lda)+0]),&ti1,&taui);
      e[i-1] = a[i*(*lda) + (i-1)];
      if(fabs(taui)>GMX_FLOAT_MIN) {
	a[i*(*lda)+(i-1)] = 1.0;
      
	ti1 = 1;
	F77_FUNC(ssymv,SSYMV)("U",&i,&taui,a,lda,&(a[i*(*lda)+0]),&ti1,&zero,tau,&ti1);

	tmp = F77_FUNC(sdot,SDOT)(&i,tau,&ti1,&(a[i*(*lda)+0]),&ti1);

	alpha = -0.5*taui*tmp;

	F77_FUNC(saxpy,SAXPY)(&i,&alpha,&(a[i*(*lda)+0]),&ti1,tau,&ti1);

	F77_FUNC(ssyr2,SSYR2)("U",&i,&minusone,&(a[i*(*lda)+0]),&ti1,tau,&ti1,a,lda);

	a[i*(*lda)+(i-1)] = e[i-1]; 

      }
      d[i] = a[i*(*lda)+i];
      tau[i-1] = taui;
    }
    d[0] = a[0];
    
  } else {
    /* lower */

    for(i=1;i<*n;i++) {

      ti1 = *n - i;
      ti2 = ( *n < i+2) ? *n : i+2;
      ti3 = 1;
      F77_FUNC(slarfg,SLARFG)(&ti1,&(a[(i-1)*(*lda)+(i)]),&(a[(i-1)*(*lda)+ti2-1]),&ti3,&taui);

      e[i-1] = a[(i-1)*(*lda) + (i)];

      if(fabs(taui)>GMX_FLOAT_MIN) {
	a[(i-1)*(*lda)+(i)] = 1.0;
      
	ti1 = *n - i;
	ti2 = 1;
	F77_FUNC(ssymv,SSYMV)(uplo,&ti1,&taui,&(a[i*(*lda)+i]),lda,&(a[(i-1)*(*lda)+i]),
	       &ti2,&zero,&(tau[i-1]),&ti2);
	
	tmp = F77_FUNC(sdot,SDOT)(&ti1,&(tau[i-1]),&ti2,&(a[(i-1)*(*lda)+i]),&ti2);

	alpha = -0.5*taui*tmp;

	F77_FUNC(saxpy,SAXPY)(&ti1,&alpha,&(a[(i-1)*(*lda)+i]),&ti2,&(tau[i-1]),&ti2);

	F77_FUNC(ssyr2,SSYR2)(uplo,&ti1,&minusone,&(a[(i-1)*(*lda)+i]),&ti2,&(tau[i-1]),&ti2,
	       &(a[(i)*(*lda)+i]),lda);

	a[(i-1)*(*lda)+(i)] = e[i-1]; 

      }
      d[i-1] = a[(i-1)*(*lda)+i-1];
      tau[i-1] = taui;
    }
    d[*n-1] = a[(*n-1)*(*lda)+(*n-1)];
 
  }
  return;
}
