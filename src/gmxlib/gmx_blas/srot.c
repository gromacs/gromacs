#include "gmx_blas.h"

void
F77_FUNC(srot,SROT)(int *n,
                    float *dx,
                    int *incx,
                    float *dy,
                    int *incy,
                    float *c,
                    float *s)
{
  int i,ix,iy;
  float dtemp;

  if(*incx!=1 || *incy!=1) {
    ix = 0;
    iy = 0;
    if(incx<0)
      ix = (1-*n)*(*incx);
    if(incy<0)
      iy = (1-*n)*(*incy);
    
    for(i=0;i<*n;i++,ix+=*incx,iy+=*incy) {
      dtemp  = (*c) * dx[ix] + (*s) * dy[iy];
      dy[iy] = (*c) * dy[iy] - (*s) * dx[ix];
      dx[ix] = dtemp;
    }

    return;

  } else {

    /* unit increments */   
    for(i=0;i<*n;i++) {
      dtemp = (*c) * dx[i] + (*s) * dy[i];
      dy[i] = (*c) * dy[i] - (*s) * dx[i];
      dx[i] = dtemp;      
    }

  }
}
