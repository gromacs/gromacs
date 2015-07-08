#include <cmath>
#include "../gmx_blas.h"

int
F77_FUNC(idamax,IDAMAX)(int *n__,
                        double *dx,
                        int *incx__)
{
  int i,ix,idxmax;
  double dmax,tmp;

  int n    = *n__;
  int incx = *incx__;
  
  if(n<1 || incx<=0)
    return -1;

  if(n==1)
    return 1;

  dmax = std::abs(dx[0]);
  idxmax = 1;

  if(incx==1) {
    for(i=1;i<n;i++) {
      tmp = std::abs(dx[i]);
      if(tmp>dmax) {
	dmax = tmp;
	idxmax = i+1;
      }
    }
  } else {
    /* Non-unit increments */
    ix = incx; /* this is really 0 + an increment */
    for(i=1;i<n;i++,ix+=incx) {
      tmp = std::abs(dx[ix]);
      if(tmp>dmax) {
	dmax = tmp;
	idxmax = ix+1;
      }
    }    
  }
  return idxmax;
}
