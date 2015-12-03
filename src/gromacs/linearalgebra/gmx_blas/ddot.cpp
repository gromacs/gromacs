#include "../gmx_blas.h"

double
F77_FUNC(ddot,DDOT)(int *n_arg,
                    double *dx,
                    int *incx_arg,
                    double *dy,
                    int *incy_arg)
{
    int i,ix,iy,m;
    int n=*n_arg;
    int incx = *incx_arg;
    int incy = *incy_arg;
    double t1;
    
    if(n<=0)
        return 0.0;
    
    t1 = 0.0;
    
    if(incx!=1 || incy!=1) {
        ix = 0;
        iy = 0;
        if(incx<0)
            ix = (1-n)*incx;
        if(incy<0)
            iy = (1-n)*incy;
        
        for(i=0;i<n;i++,ix+=incx,iy+=incy) 
            t1 += dx[ix] * dy[iy];
        
        return t1;
        
    } else {
        
        m = n%5;
        
        for(i=0;i<m;i++)
            t1 += dx[i] * dy[i];
        
        /* unroll */
        for(i=m;i<n;i+=5) 
            t1  =  t1 + dx[i] * dy[i]   
                +    dx[i+1] * dy[i+1] 
                +    dx[i+2] * dy[i+2] 
                +    dx[i+3] * dy[i+3]   
                +    dx[i+4] * dy[i+4];   
        
        return t1;
    }
}

 
