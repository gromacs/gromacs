#include <ctype.h>
#include "gmx_blas.h"

void
F77_FUNC(dgemv,DGEMV)(char *trans, 
                      int *m,
                      int *n,
                      double *alpha,
                      double *a,
                      int *lda,
                      double *x,
                      int *incx,
                      double *beta,
                      double *y,
                      int *incy)
{
  char ch=toupper(*trans);
  int lenx,leny,kx,ky;
  int i,j,jx,jy,ix,iy;
  double temp;

  if(*n<=0 || *m<=0 || (*alpha==0 && *beta==1.0))
    return;

  if(ch=='N') {
    lenx = *n;
    leny = *m;
  } else {
    lenx = *m;
    leny = *n;
  }
  
   if(*incx>0)
    kx = 1;
  else
    kx = 1 - (lenx -1)*(*incx);

  if(*incy>0)
    ky = 1;
  else
    ky = 1 - (leny -1)*(*incy);
 
  if(*beta != 1.0) {
    if(*incy==1) {
      if(*beta==0.0)
	for(i=0;i<leny;i++)
	  y[i] = 0.0;
      else
	for(i=0;i<leny;i++)
	  y[i] *= *beta;
    } else {
      /* non-unit incr. */
      iy = ky;
      if(*beta==0.0) 
	for(i=0;i<leny;i++,iy+=*incy)
	  y[iy] = 0.0;
      else
	for(i=0;i<leny;i++,iy+=*incy)
	  y[iy] *= *beta;
    }
  }
  
  if(*alpha==0.0)
    return;
  
  if(ch=='N') {
    jx = kx;
    if(*incy==1) {
      for(j=1;j<=*n;j++,jx+=*incx) 
	if(x[jx-1] != 0.0) {
	  temp = *alpha * x[jx-1];
	  for(i=1;i<=*m;i++)
	    y[i-1] += temp * a[(j-1)*(*lda)+(i-1)];
	}
    } else {
      /* non-unit y incr. */
      for(j=1;j<=*n;j++,jx+=*incx) 
	if(x[jx-1] != 0.0) {
	  temp = *alpha * x[jx-1];
	  iy = ky;
	  for(i=1;i<=*m;i++,iy+=*incy)
	    y[iy-1] += temp * a[(j-1)*(*lda)+(i-1)];
	}
    }
  } else {
    /* transpose */
    jy = ky;
    if(*incx==1) {
      for(j=1;j<=*n;j++,jy+=*incy) {
	temp = 0.0;
	for(i=1;i<=*m;i++)
	  temp += a[(j-1)*(*lda)+(i-1)] * x[i-1];
	y[jy-1] += *alpha * temp;
      }
    } else {
      /* non-unit y incr. */
      for(j=1;j<=*n;j++,jy+=*incy) {
	temp = 0.0;
	ix = kx;
	for(i=1;i<=*m;i++,ix+=*incx)
	  temp += a[(j-1)*(*lda)+(i-1)] * x[ix-1];
	y[jy-1] += *alpha * temp;
      }
    }
  }
} 
   
