#include <math.h>

#include <types/simple.h>
#include "gmx_blas.h"

void
F77_FUNC(sger,SGER)(int *m__,
                    int *n__,
                    float *alpha__,
                    float *x,
                    int *incx__,
                    float *y,
                    int *incy__,
                    float *a,
                    int *lda__)
{
    int ix,kx,jy;
    int i,j;
    float temp;
    
    int m = *m__;
    int n = *n__;
    int incx = *incx__;
    int incy = *incy__;
    int lda = *lda__;
    float alpha = *alpha__;
    
    if(m<=0 || n<=0 || fabs(alpha)<GMX_FLOAT_MIN)
        return;
    
    if(incy>0)
        jy = 0;
    else
        jy = incy * (1 - n);
    
    if(incx==1) {
        for(j=0;j<n;j++,jy+=incy)
            if(fabs(y[jy])>GMX_FLOAT_MIN) {
                temp = alpha * y[jy];
                for(i=0;i<m;i++)
                    a[j*(lda)+i] += temp*x[i];
            }
    } else {
        /* non-unit incx */
        if(incx>0) 
            kx = 0;
        else
            kx = incx * (1 - m);
        
        for(j=0;j<n;j++,jy+=incy) {
            if(fabs(y[jy])>GMX_FLOAT_MIN) {
                temp = alpha * y[jy];
                ix = kx;
                for(i=0;i<m;i++,ix+=incx)
                    a[j*(lda)+i] += temp*x[ix];
            }
        }
    }
        return;
}
