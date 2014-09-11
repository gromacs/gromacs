#include <math.h>
#include "gromacs/utility/real.h"

#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include "lapack_limits.h"


void
F77_FUNC(dlarfg,DLARFG)(int   *n,
                        double *alpha,
                        double *x,
                        int    *incx,
                        double *tau)
{
  double xnorm,t;
  int    ti1,knt,j;
  double minval,safmin,rsafmn,beta;

  if(*n<=1) {
    *tau = 0;
    return;
  }

  ti1 = *n-1;

  xnorm = F77_FUNC(dnrm2,DNRM2)(&ti1,x,incx);

  if(fabs(xnorm)<GMX_DOUBLE_MIN) {
    *tau = 0.0;
  } else {

    t = F77_FUNC(dlapy2,DLAPY2)(alpha,&xnorm);

    if(*alpha<0)
      beta = t;
    else
      beta = -t;

    minval = GMX_DOUBLE_MIN;
    
    safmin = minval*(1.0+GMX_DOUBLE_EPS) / GMX_DOUBLE_EPS;

        
    if(fabs(beta)<safmin) {

      knt = 0;
      rsafmn = 1.0 / safmin;
      
      while(fabs(beta)<safmin) {
	knt++;
	ti1 = *n-1;
	F77_FUNC(dscal,DSCAL)(&ti1,&rsafmn,x,incx);
	beta *= rsafmn;
	*alpha *= rsafmn;
      }
      
      /* safmin <= beta <= 1 now */
      ti1 = *n-1;
      xnorm = F77_FUNC(dnrm2,DNRM2)(&ti1,x,incx);
      t = F77_FUNC(dlapy2,DLAPY2)(alpha,&xnorm);
      
      if(*alpha<0)
	beta = t;
      else
	beta = -t;
      
      *tau = (beta-*alpha)/beta;

      ti1= *n-1;
      t = 1.0/(*alpha-beta);
      F77_FUNC(dscal,DSCAL)(&ti1,&t,x,incx);
   
      *alpha = beta;
      for(j=0;j<knt;j++)
	*alpha *= safmin;
    } else {
      *tau = (beta-*alpha)/beta;
      ti1= *n-1;
      t = 1.0/(*alpha-beta);
      F77_FUNC(dscal,DSCAL)(&ti1,&t,x,incx);
      *alpha = beta;
    }
  }
   
  return;
}
