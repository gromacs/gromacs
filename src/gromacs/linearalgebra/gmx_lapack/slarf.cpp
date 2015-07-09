#include <cctype>
#include <cmath>

#include "../gmx_blas.h"
#include "../gmx_lapack.h"

#include "gromacs/utility/real.h"

void
F77_FUNC(slarf,SLARF)(const char *side,
       int *m,
       int *n,
       float *v,
       int *incv,
       float *tau,
       float *c,
       int *ldc,
       float *work)
{
  const char ch=std::toupper(*side);
  float one = 1.0;
  float zero = 0.0;
  float minustau = -(*tau);
  int i1 = 1;


  if(ch=='L') {
    if(std::abs(*tau)>GMX_FLOAT_MIN) {
      F77_FUNC(sgemv,SGEMV)("T",m,n,&one,c,ldc,v,incv,&zero,work,&i1);
      F77_FUNC(sger,SGER)(m,n,&minustau,v,incv,work,&i1,c,ldc);
    }
  } else {
    if(std::abs(*tau)>GMX_FLOAT_MIN) {
      F77_FUNC(sgemv,SGEMV)("N",m,n,&one,c,ldc,v,incv,&zero,work,&i1);
      F77_FUNC(sger,SGER)(m,n,&minustau,work,&i1,v,incv,c,ldc);
    }
  }
  return;
}
