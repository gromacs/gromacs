#include <math.h>
#include <types/simple.h>
#include "gmx_lapack.h"

void
F77_FUNC(slassq,SLASSQ)(int *n,
                        float *x,
                        int *incx,
                        float *scale,
                        float *sumsq)
{
  int ix;
  float absxi,t;

  if(*n>0) {
    for(ix=0;ix<=(*n-1)*(*incx);ix+=*incx) {
      if(fabs(x[ix])>GMX_FLOAT_MIN) {
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
