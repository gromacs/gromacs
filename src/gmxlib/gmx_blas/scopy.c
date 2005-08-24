#include "gmx_blas.h"

void
F77_FUNC(scopy,SCOPY)(int *n__,
                      float *dx,
                      int *incx__,
                      float *dy,
                      int *incy__)
{
    int i,ix,iy;

    int n= *n__;
    int incx = *incx__;
    int incy = *incy__;
    
    if(incx!=1 || incy!=1) 
    {
        ix = 0;
        iy = 0;
        if(incx<0)
            ix = (1-n)*(incx);
        if(incy<0)
            iy = (1-n)*(incy);
        
        for(i=0;i<n;i++,ix+=incx,iy+=incy) 
            dy[iy] = dx[ix];
        
        return;
        
    } else {
        
        /* unroll */
        
        for(i=0;i<(n-8);i+=8) {
            dy[i]   = dx[i];
            dy[i+1] = dx[i+1];
            dy[i+2] = dx[i+2];
            dy[i+3] = dx[i+3];
            dy[i+4] = dx[i+4];
            dy[i+5] = dx[i+5];
            dy[i+6] = dx[i+6];
            dy[i+7] = dx[i+7];
        }
        /* continue with current value of i */
        for(;i<n;i++)
            dy[i] = dx[i];
    }
}
