#include "gmx_blas.h"


void
F77_FUNC(daxpy,DAXPY)(int   *   n_arg,
                      double *   da_arg,
                      double *   dx,
                      int *      incx_arg,
                      double *   dy,
                      int *      incy_arg)
{
  int i,ix,iy;
  int n=*n_arg;
  double da=*da_arg;
  int incx = *incx_arg;
  int incy = *incy_arg;

  if (n<=0)
    return;

  if(incx!=1 || incy!=1) {
    ix = 0;
    iy = 0;
    if(incx<0)
      ix = (1-n)*incx;
    if(incy<0)
      iy = (1-n)*incy;
    
    for(i=0;i<n;i++,ix+=incx,iy+=incy) 
      dy[iy] += da*dx[ix];

    return;

  } else {

    /* unroll */
    
    for(i=0;i<(n-4);i+=4) {
      dy[i]   += da*dx[i];
      dy[i+1] += da*dx[i+1];
      dy[i+2] += da*dx[i+2];
      dy[i+3] += da*dx[i+3];
    }
    /* continue with current value of i */
    for(;i<n;i++)
      dy[i]   += da*dx[i];
  }
}
