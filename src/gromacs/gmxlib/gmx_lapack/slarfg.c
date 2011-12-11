#include <math.h>
#include <types/simple.h>

#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "lapack_limits.h"


void
F77_FUNC(slarfg,SLARFG)(int   *n,
                        float *alpha,
                        float *x,
                        int    *incx,
                        float *tau)
{
  float xnorm,t;
  int    ti1,knt,j;
  float minval,safmin,rsafmn,beta;

  if(*n<=1) {
    *tau = 0;
    return;
  }

  ti1 = *n-1;

  xnorm = F77_FUNC(snrm2,SNRM2)(&ti1,x,incx);

  if(fabs(xnorm)<GMX_FLOAT_MIN) {
    *tau = 0.0;
  } else {

    t = F77_FUNC(slapy2,SLAPY2)(alpha,&xnorm);

    if(*alpha<0)
      beta = t;
    else
      beta = -t;

    minval = GMX_FLOAT_MIN;
    
    safmin = minval*(1.0+GMX_FLOAT_EPS) / GMX_FLOAT_EPS;

        
    if(fabs(beta)<safmin) {

      knt = 0;
      rsafmn = 1.0 / safmin;
      
      while(fabs(beta)<safmin) {
	knt++;
	ti1 = *n-1;
	F77_FUNC(sscal,SSCAL)(&ti1,&rsafmn,x,incx);
	beta *= rsafmn;
	*alpha *= rsafmn;
      }
      
      /* safmin <= beta <= 1 now */
      ti1 = *n-1;
      xnorm = F77_FUNC(snrm2,SNRM2)(&ti1,x,incx);
      t = F77_FUNC(slapy2,SLAPY2)(alpha,&xnorm);
      
      if(*alpha<0)
	beta = t;
      else
	beta = -t;
      
      *tau = (beta-*alpha)/beta;

      ti1= *n-1;
      t = 1.0/(*alpha-beta);
      F77_FUNC(sscal,SSCAL)(&ti1,&t,x,incx);
   
      *alpha = beta;
      for(j=0;j<knt;j++)
	*alpha *= safmin;
    } else {
      *tau = (beta-*alpha)/beta;
      ti1= *n-1;
      t = 1.0/(*alpha-beta);
      F77_FUNC(sscal,SSCAL)(&ti1,&t,x,incx);
      *alpha = beta;
    }
  }
   
  return;
}
