#include <math.h>
#include <ctype.h>

#include <types/simple.h>

#include "gmx_blas.h"

void
F77_FUNC(dgemv,DGEMV)(char *trans, 
       int *m__,
       int *n__,
       double *alpha__,
       double *a,
       int *lda__,
       double *x,
       int *incx__,
       double *beta__,
       double *y,
       int *incy__)
{
  char ch=toupper(*trans);
  int lenx,leny,kx,ky;
  int i,j,jx,jy,ix,iy;
  double temp;

  int m = *m__;
  int n = *n__;
  double alpha = *alpha__;
  double beta = *beta__;
  int incx = *incx__;
  int incy = *incy__;
  int lda = *lda__;
  
  if(n<=0 || m<=0 || (fabs(alpha)<GMX_DOUBLE_MIN && fabs(beta-1.0)<GMX_DOUBLE_EPS))
    return;

  if(ch=='N') {
    lenx = n;
    leny = m;
  } else {
    lenx = m;
    leny = n;
  }
  
   if(incx>0)
    kx = 1;
  else
    kx = 1 - (lenx -1)*(incx);

  if(incy>0)
    ky = 1;
  else
    ky = 1 - (leny -1)*(incy);
 
  if(fabs(beta-1.0)>GMX_DOUBLE_EPS) {
    if(incy==1) {
      if(fabs(beta)<GMX_DOUBLE_MIN)
	for(i=0;i<leny;i++)
	  y[i] = 0.0;
      else
	for(i=0;i<leny;i++)
	  y[i] *= beta;
    } else {
      /* non-unit incr. */
      iy = ky;
      if(fabs(beta)<GMX_DOUBLE_MIN) 
	for(i=0;i<leny;i++,iy+=incy)
	  y[iy] = 0.0;
      else
	for(i=0;i<leny;i++,iy+=incy)
	  y[iy] *= beta;
    }
  }
  
  if(fabs(alpha)<GMX_DOUBLE_MIN)
    return;
  
  if(ch=='N') {
    jx = kx;
    if(incy==1) {
      for(j=1;j<=n;j++,jx+=incx) 
	if(fabs(x[jx-1])>GMX_DOUBLE_MIN) {
	  temp = alpha * x[jx-1];
	  for(i=1;i<=m;i++)
	    y[i-1] += temp * a[(j-1)*(lda)+(i-1)];
	}
    } else {
      /* non-unit y incr. */
      for(j=1;j<=n;j++,jx+=incx) 
	if(fabs(x[jx-1])>GMX_DOUBLE_MIN) {
	  temp = alpha * x[jx-1];
	  iy = ky;
	  for(i=1;i<=m;i++,iy+=incy)
	    y[iy-1] += temp * a[(j-1)*(lda)+(i-1)];
	}
    }
  } else {
    /* transpose */
    jy = ky;
    if(incx==1) {
      for(j=1;j<=n;j++,jy+=incy) {
	temp = 0.0;
	for(i=1;i<=m;i++)
	  temp += a[(j-1)*(lda)+(i-1)] * x[i-1];
	y[jy-1] += alpha * temp;
      }
    } else {
      /* non-unit y incr. */
      for(j=1;j<=n;j++,jy+=incy) {
	temp = 0.0;
	ix = kx;
	for(i=1;i<=m;i++,ix+=incx)
	  temp += a[(j-1)*(lda)+(i-1)] * x[ix-1];
	y[jy-1] += alpha * temp;
      }
    }
  }
} 
   
