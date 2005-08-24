#include <math.h>
#include <ctype.h>

#include "gmx_blas.h"
#include "gmx_lapack.h"

#include <types/simple.h>

void
F77_FUNC(dlarf,DLARF)(char *side,
       int *m,
       int *n,
       double *v,
       int *incv,
       double *tau,
       double *c,
       int *ldc,
       double *work)
{
  char ch=toupper(*side);
  double one = 1.0;
  double zero = 0.0;
  double minustau = -(*tau);
  int i1 = 1;


  if(ch=='L') {
    if(fabs(*tau)>GMX_DOUBLE_MIN) {
      F77_FUNC(dgemv,DGEMV)("T",m,n,&one,c,ldc,v,incv,&zero,work,&i1);
      F77_FUNC(dger,DGER)(m,n,&minustau,v,incv,work,&i1,c,ldc);
    }
  } else {
    if(fabs(*tau)>GMX_DOUBLE_MIN) {
      F77_FUNC(dgemv,DGEMV)("N",m,n,&one,c,ldc,v,incv,&zero,work,&i1);
      F77_FUNC(dger,DGER)(m,n,&minustau,work,&i1,v,incv,c,ldc);
    }
  }
  return;
}
