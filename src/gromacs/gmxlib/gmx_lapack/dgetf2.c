#include <math.h>
#include <types/simple.h>

#include "gmx_blas.h"
#include "gmx_lapack.h"


void
F77_FUNC(dgetf2,DGETF2)(int *m,
	int *n,
	double *a,
	int *lda,
	int *ipiv,
	int *info)
{
  int j,jp,k,t1,t2,t3;
  double one,minusone;
  double tmp;

  one = 1.0;
  minusone = -1.0;

  if(*m<=0 || *n<=0)
    return;

  k = (*m < *n) ? *m : *n;
  for(j=1;j<=k;j++) {
    t1 = *m-j+1;
    t2 = 1;
    jp = j - 1 + F77_FUNC(idamax,IDAMAX)(&t1,&(a[(j-1)*(*lda)+(j-1)]),&t2);
    ipiv[j-1] = jp;
    if( fabs(a[(j-1)*(*lda)+(jp-1)])>GMX_DOUBLE_MIN ) {
      if(jp != j)
	F77_FUNC(dswap,DSWAP)(n,&(a[ j-1 ]),lda,&(a[ jp-1 ]),lda);
      
      if(j<*m) {
	t1 = *m-j;
	t2 = 1;
	tmp = 1.0/a[(j-1)*(*lda)+(j-1)];
	F77_FUNC(dscal,DSCAL)(&t1,&tmp,&(a[(j-1)*(*lda)+(j)]),&t2);
      }
    } else {
      *info = j;
    }

    if(j<k) {
      t1 = *m-j;
      t2 = *n-j;
      t3 = 1;
      F77_FUNC(dger,DGER)(&t1,&t2,&minusone,&(a[(j-1)*(*lda)+(j)]),&t3,
	    &(a[(j)*(*lda)+(j-1)]),lda, &(a[(j)*(*lda)+(j)]),lda);
    }
  }
  return;
}
