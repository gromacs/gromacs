#include <ctype.h>
#include "gmx_blas.h"

void
F77_FUNC(ssyr2,SSYR2)(char *    uplo,
                      int *     n,
                      float *  alpha,
                      float *  x,
                      int *     incx,
                      float *  y,
                      int *     incy,
                      float *  a,
                      int *     lda)
{
  int kx,ky,ix,iy,jx,jy,j,i;
  float temp1,temp2;
  char ch=toupper(*uplo);


  if(*n<=0 || *alpha==0 || *incx==0 || *incy==0 ||
     (ch != 'U' && ch != 'L'))
    return;
  
  jx = jy = kx = ky = 0;

  /* init start points for non-unit increments */
  if(*incx!=1 || *incy!=1) {
    if(*incx>0)
      kx = 1;
    else
      kx = 1 - (*n - 1)*(*incx);
    if(*incy>0)
      ky = 1;
    else
      ky = 1 - (*n - 1)*(*incy);
    
    jx = kx;
    jy = ky;
  }

  if(ch == 'U') {
    /* Data in upper part of A */
    if(*incx==1 && *incy==1) {
      /* Unit increments for both x and y */
      for(j=1;j<=*n;j++) {
	if( x[j-1] != 0  || y[j-1] != 0 ) {
	  temp1 = *alpha * y[j-1];
	  temp2 = *alpha * x[j-1];
	  for(i=1;i<=j;i++)
	    a[(j-1)*(*lda)+(i-1)] += x[i-1]*temp1 + y[i-1]*temp2;
	}
      }
    } else {
	
      /* non-unit increments */
      for(j=1;j<=*n;j++) {
	
	if( x[jx-1] != 0 || y[jy-1] != 0 ) {
	  temp1 = *alpha * y[jy-1];
	  temp2 = *alpha * x[jx-1];
	  ix = kx;
	  iy = ky;
	  for(i=1;i<=j;i++) {
	    a[(j-1)*(*lda)+(i-1)] += x[ix-1]*temp1 + y[iy-1]*temp2;
	    ix += *incx;
	    iy += *incy;
	  }
	}
	jx += *incx;
	jy += *incy;
      }
    }
  } else {
    /* Data in lower part of A */
    if(*incx==1 && *incy==1) {
      /* Unit increments for both x and y */
      for(j=1;j<=*n;j++) {
	if( x[j-1] != 0  || y[j-1] != 0 ) {
	  temp1 = *alpha * y[j-1];
	  temp2 = *alpha * x[j-1];
	  for(i=j;i<=*n;i++)
	    a[(j-1)*(*lda)+(i-1)] += x[i-1]*temp1 + y[i-1]*temp2;
	}
      }
    } else {
	
      /* non-unit increments */
      for(j=1;j<=*n;j++) {
	
	if( x[jx-1] != 0 || y[jy-1] != 0 ) {
	  temp1 = *alpha * y[jy-1];
	  temp2 = *alpha * x[jx-1];
	  ix = jx;
	  iy = jy;
	  for(i=j;i<=*n;i++) {
	    a[(j-1)*(*lda)+(i-1)] += x[ix-1]*temp1 + y[iy-1]*temp2;
	    ix += *incx;
	    iy += *incy;
	  }
	}
	jx += *incx;
	jy += *incy;
      }
    }
  }

  return;
}
