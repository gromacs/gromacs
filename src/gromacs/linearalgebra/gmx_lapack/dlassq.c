#include <math.h>
#include "gromacs/utility/real.h"
#include "../gmx_lapack.h"

void
F77_FUNC(dlassq,DLASSQ)(int *n,
                        double *x,
                        int *incx,
                        double *scale,
                        double *sumsq)
{
  int ix;
  double absxi,t;

  if(*n>0) {
    for(ix=0;ix<=(*n-1)*(*incx);ix+=*incx) {
      if(fabs(x[ix])>GMX_DOUBLE_MIN) {
	absxi = fabs(x[ix]);
	if(*scale<absxi) {
	  t = *scale/absxi;
	  t = t*t;
	  *sumsq = 1.0 + (*sumsq)*t;
	  *scale = absxi;
	} else {
	  t = absxi/(*scale);
	  *sumsq += t*t;
	}
      }
    }
  }
  return;
}
