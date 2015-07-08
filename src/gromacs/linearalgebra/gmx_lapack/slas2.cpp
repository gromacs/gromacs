#include <cmath>
#include "gromacs/utility/real.h"

#include "../gmx_lapack.h"

void
F77_FUNC(slas2,SLAS2)(float *f,
       float *g,
       float *h,
       float *ssmin,
       float *ssmax)
{
  float fa = std::abs(*f);
  float ga = std::abs(*g);
  float ha = std::abs(*h);
  float fhmin,fhmax,tmax,tmin,tmp1,tmp2;
  float as,at,au,c;

  fhmin = (fa<ha) ? fa : ha;
  fhmax = (fa>ha) ? fa : ha;
  
  if(std::abs(fhmin)<GMX_FLOAT_MIN) {
    *ssmin = 0.0;
    if(std::abs(fhmax)<GMX_FLOAT_MIN) 
      *ssmax = ga;
    else {
      tmax = (fhmax>ga) ? fhmax : ga;
      tmin = (fhmax<ga) ? fhmax : ga;
      tmp1 = tmin / tmax;
      tmp1 = tmp1 * tmp1;
      *ssmax = tmax* std::sqrt(1.0 + tmp1);
    }
  } else {
    if(ga<fhmax) {
      as = 1.0 + fhmin / fhmax;
      at = (fhmax-fhmin) / fhmax;
      au = (ga/fhmax);
      au = au * au;
      c = 2.0 / (  std::sqrt(as*as+au) + std::sqrt(at*at+au) );
      *ssmin = fhmin * c;
      *ssmax = fhmax / c;
    } else {
      au = fhmax / ga;
      if(std::abs(au)<GMX_FLOAT_MIN) {
	*ssmin = (fhmin*fhmax)/ga;
	*ssmax = ga;
      } else {
	as = 1.0 + fhmin / fhmax;
	at = (fhmax-fhmin)/fhmax;
	tmp1 = as*au;
	tmp2 = at*au;
	c = 1.0 / (  std::sqrt(1.0+tmp1*tmp1) + std::sqrt(1.0+tmp2*tmp2));
	*ssmin = (fhmin*c)*au;
	*ssmin = *ssmin + *ssmin;
	*ssmax = ga / (c+c);
      }
    }
  }
  return;
}
