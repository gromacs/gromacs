#include <math.h>
#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include "lapack_limits.h"

#include "gromacs/utility/real.h"

void
F77_FUNC(slasq1,SLASQ1)(int *n,
	float *d,
	float *e,
	float *work,
	int *info)
{
  float sigmx = 0.0;
  int i,j,k,iinfo;
  float minval,safemin;
  float dtemp,scale;
  float eps;

  eps = GMX_FLOAT_EPS;
  minval = GMX_FLOAT_MIN;
  safemin = minval*(1.0+GMX_FLOAT_EPS);
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
  
  if(fabs(sigmx)<GMX_FLOAT_MIN) {
    F77_FUNC(slasrt,SLASRT)("D",n,d,&iinfo);
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
  F77_FUNC(scopy,SCOPY)(n,d,&i,work,&j);
  k = *n-1;
  F77_FUNC(scopy,SCOPY)(&k,e,&i,work+1,&j);
  i = 0;
  j = 2*(*n)-1;
  k = 1;
  F77_FUNC(slascl,SLASCL)("G",&i,&i,&sigmx,&scale,&j,&k,work,&j,&iinfo);


  /* Compute q and e elements */
  for(i=0;i<2*(*n)-1;i++)
    work[i] = work[i]*work[i];

  work[2*(*n)-1] = 0.0;

  F77_FUNC(slasq2,SLASQ2)(n,work,info);

  j = 0;
  k = 1;
  if(*info==0) {
    for(i=0;i<*n;i++)
      d[i]=sqrt(work[i]);
    F77_FUNC(slascl,SLASCL)("G",&j,&j,&scale,&sigmx,n,&k,d,n,&iinfo);
  }
  return;
}
