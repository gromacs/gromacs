#include <math.h>


#include "gromacs/utility/real.h"
#include "../gmx_blas.h"

float
F77_FUNC(snrm2,SNRM2)(int  *     n__,
                      float *    x,
                      int    *    incx__)
{
    int ix,max_ix;
    float ssq,scale,absxi,t;
    
    int n = *n__;
    int incx = *incx__;
    
    if(n<1 || incx<1)
        return 0;
    else if (n==1) {
        t = x[0];
        if(t>=0)
            return t;
        else 
            return -t;
    }
    
    scale = 0.0;
    ssq   = 1.0;
    
    max_ix = 1+(n-1)*(incx);
    for(ix=1;ix<=max_ix;ix+=incx) {
        t = x[ix-1];
        if(fabs(t)>GMX_FLOAT_MIN) {
            absxi = (t>=0) ? t : (-t);
            if(scale<absxi) {
                t = scale/absxi;
                t = t*t;
                ssq = ssq*t + 1.0;
                scale = absxi;
            } else {
                t = absxi/scale;
                ssq += t*t;
            }
        }
    }
    return scale*sqrt(ssq);
    
}


 
