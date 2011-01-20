#include "gmx_blas.h"

void 
F77_FUNC(sscal,SSCAL)(int  *n__,
                      float *fact__,
                      float *dx,
                      int   *incx__)
{
    int nincx,i;

    int n = *n__;
    float fact = *fact__;
    int incx = *incx__;
    
    if(n<=0 || incx<=0)
        return;
    
    if(incx==1) {
        /* Unrool factor 5 */
        for(i=0;i<(n-5);i+=5) {
            dx[i]   *= fact;
            dx[i+1] *= fact;
            dx[i+2] *= fact;
            dx[i+3] *= fact;
            dx[i+4] *= fact;
        }    
        /* continue with current value of i */
        for(;i<n;i++)
            dx[i]   *= fact;
        
        return;
    } else {
        /* inc != 1 */
        nincx = n * (incx);
        for (i=0;i<nincx;i+=incx)
            dx[i] *= fact;
        
        return;
    } 
    
}
