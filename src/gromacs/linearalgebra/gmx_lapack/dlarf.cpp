#include <cctype>
#include <cmath>

#include "../gmx_blas.h"
#include "../gmx_lapack.h"

#include "gromacs/utility/real.h"

void
F77_FUNC(dlarf,DLARF)(const char *side,
       int *m,
       int *n,
       double *v,
       int *incv,
       double *tau,
       double *c,
       int *ldc,
       double *work)
{
  const char ch=std::toupper(*side);
  double one = 1.0;
  double zero = 0.0;
  double minustau = -(*tau);
  int i1 = 1;


  if(ch=='L') {
    if(std::abs(*tau)>GMX_DOUBLE_MIN) {
      F77_FUNC(dgemv,DGEMV)("T",m,n,&one,c,ldc,v,incv,&zero,work,&i1);
      F77_FUNC(dger,DGER)(m,n,&minustau,v,incv,work,&i1,c,ldc);
    }
  } else {
    if(std::abs(*tau)>GMX_DOUBLE_MIN) {
      F77_FUNC(dgemv,DGEMV)("N",m,n,&one,c,ldc,v,incv,&zero,work,&i1);
      F77_FUNC(dger,DGER)(m,n,&minustau,work,&i1,v,incv,c,ldc);
    }
  }
  return;
}
