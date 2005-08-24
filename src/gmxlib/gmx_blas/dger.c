#include <math.h>

#include <types/simple.h>

#include "gmx_blas.h"

void
F77_FUNC(dger,DGER)(int *m__,
                    int *n__,
                    double *alpha__,
                    double *x,
                    int *incx__,
                    double *y,
                    int *incy__,
                    double *a,
                    int *lda__)
{
    int ix,kx,jy;
    int i,j;
    double temp;
    
    
    int m = *m__;
    int n = *n__;
    int incx = *incx__;
    int incy = *incy__;
    int lda = *lda__;
    double alpha = *alpha__;
    
    if(m<=0 || n<=0 || fabs(alpha)<GMX_DOUBLE_MIN)
        return;
    
    if(incy>0)
        jy = 0;
    else
        jy = incy * (1 - n);
    
    if(incx==1) {
        for(j=0;j<n;j++,jy+=incy)
            if(fabs(y[jy])>GMX_DOUBLE_MIN) {
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
            if(fabs(y[jy])>GMX_DOUBLE_MIN) {
                temp = alpha * y[jy];
                ix = kx;
                for(i=0;i<m;i++,ix+=incx)
                    a[j*(lda)+i] += temp*x[ix];
            }
        }
    }
        return;
}
