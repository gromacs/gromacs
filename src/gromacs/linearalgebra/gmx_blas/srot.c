#include "../gmx_blas.h"

void
F77_FUNC(srot,SROT)(int *n__,
                    float *dx,
                    int *incx__,
                    float *dy,
                    int *incy__,
                    float *c__,
                    float *s__)
{
  int i,ix,iy;
  float dtemp;

  int n = *n__;
  int incx = *incx__;
  int incy = *incy__;
  float c = *c__;
  float s = *s__;
  
  if(incx!=1 || incy!=1) {
    ix = 0;
    iy = 0;
    if(incx<0)
      ix = (1-n)*(incx);
    if(incy<0)
      iy = (1-n)*(incy);
    
    for(i=0;i<n;i++,ix+=incx,iy+=incy) {
      dtemp  = (c) * dx[ix] + (s) * dy[iy];
      dy[iy] = (c) * dy[iy] - (s) * dx[ix];
      dx[ix] = dtemp;
    }

    return;

  } else {

    /* unit increments */   
    for(i=0;i<n;i++) {
      dtemp = (c) * dx[i] + (s) * dy[i];
      dy[i] = (c) * dy[i] - (s) * dx[i];
      dx[i] = dtemp;      
    }

  }
}
