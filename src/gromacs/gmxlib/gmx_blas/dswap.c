#include "gmx_blas.h"

void
F77_FUNC(dswap,DSWAP)(int *n__,
                      double *dx,
                      int *incx__,
                      double *dy,
                      int *incy__)
{
  int i,ix,iy;
  double d1,d2,d3;

  int n = *n__;
  int incx = *incx__;
  int incy = *incy__;
  
  if(n<=0)
    return;

  if(incx==1 && incy==1) {
    for(i=0;i<(n-3);i+=3) {
      d1      = dx[i];
      d2      = dx[i+1];
      d3      = dx[i+2];
      dx[i]   = dy[i];
      dx[i+1] = dy[i+1];
      dx[i+2] = dy[i+2];
      dy[i]   = d1;
      dy[i+1] = d2;
      dy[i+2] = d3;
    }
    /* continue with last i value */
    for(;i<n;i++) {
      d1      = dx[i];
      dx[i]   = dy[i];
      dy[i]   = d1;
    }

  } else {
    ix = 0;
    iy = 0;
    if(incx<0)
      ix = incx * (1 - n);
    if(incy<0)
      iy = incy * (1 - n);

    for(i=0;i<n;i++,ix+=incx,iy+=incy) {
      d1     = dx[ix];
      dx[ix] = dy[iy];
      dy[iy] = d1;
    }
  }
  return;
}
 
