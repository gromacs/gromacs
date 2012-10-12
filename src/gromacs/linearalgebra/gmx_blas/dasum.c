#include <math.h>
#include "../gmx_blas.h"

double
F77_FUNC(dasum,DASUM)(int *n__, 
                      double *dx, 
                      int *incx__)
{
    int i__1, i__2;
    
    int i__, m, mp1;
    double dtemp;
    int nincx;
    
    int n = *n__;
    int incx = *incx__;
    
    --dx;
    
    dtemp = 0.;
    if (n <= 0 || incx <= 0) {
        return 0.0;
    }
    if (incx != 1) {
        nincx = n * incx;
        i__1 = nincx;
        i__2 = incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
            dtemp += fabs(dx[i__]);
        }
        return dtemp;
    }
    
    m = n % 6;
    if (m != 0) {
        i__2 = m;
        for (i__ = 1; i__ <= i__2; ++i__) {
            dtemp += fabs(dx[i__]);
        }
        if (n < 6) {
            return dtemp;
        }
    }
    mp1 = m + 1;
    i__2 = n;
    for (i__ = mp1; i__ <= i__2; i__ += 6) {
        dtemp = dtemp + fabs(dx[i__]) + fabs(dx[i__ + 1]) + 
        fabs(dx[i__ + 2]) + fabs(dx[i__+ 3]) + fabs(dx[i__ + 4]) +
        fabs(dx[i__ + 5]);
    }
    return dtemp;
}


