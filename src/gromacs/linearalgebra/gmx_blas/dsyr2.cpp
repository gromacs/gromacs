#include <cctype>
#include <cmath>

#include "gromacs/utility/real.h"

#include "../gmx_blas.h"

void
F77_FUNC(dsyr2,DSYR2)(const char *    uplo,
                      int *     n__,
                      double *  alpha__,
                      double *  x,
                      int *     incx__,
                      double *  y,
                      int *     incy__,
                      double *  a,
                      int *     lda__)
{
    int kx,ky,ix,iy,jx,jy,j,i;
    double temp1,temp2;
    const char ch=std::toupper(*uplo);
    
    int n = *n__;
    int lda = *lda__;
    int incx = *incx__;
    int incy = *incy__;
    float alpha = *alpha__;
    
    
    if(n<=0 || std::abs(alpha)<GMX_DOUBLE_MIN || incx==0 || incy==0 ||
       (ch != 'U' && ch != 'L'))
        return;
    
    jx = jy = kx = ky = 0;
    
    /* init start points for non-unit increments */
    if(incx!=1 || incy!=1) {
        if(incx>0)
            kx = 1;
        else
            kx = 1 - (n - 1)*(incx);
        if(incy>0)
            ky = 1;
        else
            ky = 1 - (n - 1)*(incy);
        
        jx = kx;
        jy = ky;
    }
    
    if(ch == 'U') {
        /* Data in upper part of A */
        if(incx==1 && incy==1) {
            /* Unit increments for both x and y */
            for(j=1;j<=n;j++) {
                if( std::abs(x[j-1])>GMX_DOUBLE_MIN  || std::abs(y[j-1])>GMX_DOUBLE_MIN ) {
                    temp1 = alpha * y[j-1];
                    temp2 = alpha * x[j-1];
                    for(i=1;i<=j;i++)
                        a[(j-1)*(lda)+(i-1)] += x[i-1]*temp1 + y[i-1]*temp2;
                }
            }
        } else {
            
            /* non-unit increments */
            for(j=1;j<=n;j++) {
                
                if( std::abs(x[jx-1])>GMX_DOUBLE_MIN || std::abs(y[jy-1])>GMX_DOUBLE_MIN ) {
                    temp1 = alpha * y[jy-1];
                    temp2 = alpha * x[jx-1];
                    ix = kx;
                    iy = ky;
                    for(i=1;i<=j;i++) {
                        a[(j-1)*(lda)+(i-1)] += x[ix-1]*temp1 + y[iy-1]*temp2;
                        ix += incx;
                        iy += incy;
                    }
                }
                jx += incx;
                jy += incy;
            }
        }
    } else {
        /* Data in lower part of A */
        if(incx==1 && incy==1) {
            /* Unit increments for both x and y */
            for(j=1;j<=n;j++) {
                if( std::abs(x[j-1])>GMX_DOUBLE_MIN  || std::abs(y[j-1])>GMX_DOUBLE_MIN ) {
                    temp1 = alpha * y[j-1];
                    temp2 = alpha * x[j-1];
                    for(i=j;i<=n;i++)
                        a[(j-1)*(lda)+(i-1)] += x[i-1]*temp1 + y[i-1]*temp2;
                }
            }
        } else {
            
            /* non-unit increments */
            for(j=1;j<=n;j++) {
                
                if( std::abs(x[jx-1])>GMX_DOUBLE_MIN || std::abs(y[jy-1])>GMX_DOUBLE_MIN ) {
                    temp1 = alpha * y[jy-1];
                    temp2 = alpha * x[jx-1];
                    ix = jx;
                    iy = jy;
                    for(i=j;i<=n;i++) {
                        a[(j-1)*(lda)+(i-1)] += x[ix-1]*temp1 + y[iy-1]*temp2;
                        ix += incx;
                        iy += incy;
                    }
                }
                jx += incx;
                jy += incy;
            }
        }
    }
    
    return;
}
