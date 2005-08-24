#include <math.h>
#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "lapack_limits.h"

#include <types/simple.h>

void
F77_FUNC(dlasq1,DLASQ1)(int *n,
	double *d,
	double *e,
	double *work,
	int *info)
{
  double sigmx = 0.0;
  int i,j,k,iinfo;
  double minval,safemin;
  double dtemp,scale;
  double eps;

  eps = GMX_DOUBLE_EPS;
  minval = GMX_DOUBLE_MIN;
  safemin = minval*(1.0+GMX_DOUBLE_EPS);
  *info = 0;

  if(*n<0) {
    *info = -2;
    return;
  }
  
  for(i=0;i<*n-1;i++) {
    d[i] = fabs(d[i]);
    dtemp = fabs(e[i]);
    if(dtemp>sigmx)
      sigmx=dtemp;
  }
  d[*n-1] = fabs(d[*n-1]);
  
  if(fabs(sigmx)<GMX_DOUBLE_MIN) {
    F77_FUNC(dlasrt,DLASRT)("D",n,d,&iinfo);
    return;
  }

  for(i=0;i<*n;i++) {
    if(d[i]>sigmx)
      sigmx=d[i];
  }

  /* Copy d and e into work (z format) and scale.
   * Squaring input data makes scaling by a power of the
   * radix pointless.
   */
  scale = sqrt(eps/safemin);
  i = 1;
  j = 2;
  F77_FUNC(dcopy,DCOPY)(n,d,&i,work,&j);
  k = *n-1;
  F77_FUNC(dcopy,DCOPY)(&k,e,&i,work+1,&j);
  i = 0;
  j = 2*(*n)-1;
  k = 1;
  F77_FUNC(dlascl,DLASCL)("G",&i,&i,&sigmx,&scale,&j,&k,work,&j,&iinfo);


  /* Compute q and e elements */
  for(i=0;i<2*(*n)-1;i++)
    work[i] = work[i]*work[i];

  work[2*(*n)-1] = 0.0;

  F77_FUNC(dlasq2,DLASQ2)(n,work,info);

  j = 0;
  k = 1;
  if(*info==0) {
    for(i=0;i<*n;i++)
      d[i]=sqrt(work[i]);
    F77_FUNC(dlascl,DLASCL)("G",&j,&j,&scale,&sigmx,n,&k,d,n,&iinfo);
  }
  return;
}
