#include <cmath>
#include "../gmx_lapack.h"
#include "lapack_limits.h"

#include "gromacs/utility/real.h"

void
F77_FUNC(dlartg,DLARTG)(double *f,
	double *g,
	double *cs,
	double *sn,
	double *r)
{
  double minval,safemin, safemin2, safemx2, eps;
  double f1,g1,f1a,g1a,scale;
  int i,n,count;

  eps = GMX_DOUBLE_EPS;
  minval = GMX_DOUBLE_MIN;
  safemin = minval*(1.0+eps);
  n = static_cast<int>(0.5*std::log( safemin/eps ) / std::log(2.0));
  safemin2 = std::pow(2.0,static_cast<double>(n));

  safemx2 = 1.0 / safemin2;

  if(std::abs(*g)<GMX_DOUBLE_MIN) {
    *cs = 1.0;
    *sn = 0.0;
    *r = *f;
  } else if (std::abs(*f)<GMX_DOUBLE_MIN) {
    *cs = 0.0;
    *sn = 1.0;
    *r = *g;
  } else {
    f1 = *f;
    g1 = *g;
    f1a = std::abs(f1);
    g1a = std::abs(g1);
    scale = (f1a > g1a) ? f1a : g1a;
    if(scale >= safemx2) {
      count = 0;
      while(scale >= safemx2) {
	count++;
	f1 *= safemin2;
	g1 *= safemin2;
	f1a = std::abs(f1);
	g1a = std::abs(g1);
	scale = (f1a > g1a) ? f1a : g1a;
      }
      *r =  std::sqrt(f1*f1 + g1*g1);
      *cs = f1 / *r;
      *sn = g1 / *r;
      for(i=0;i<count;i++)
	*r *= safemx2;
    } else if (scale<=safemin2) {
      count = 0;
      while(scale <= safemin2) {
	count++;
	f1 *= safemx2;
	g1 *= safemx2;
	f1a = std::abs(f1);
	g1a = std::abs(g1);
	scale = (f1a > g1a) ? f1a : g1a;
      }
      *r =  std::sqrt(f1*f1 + g1*g1);
      *cs = f1 / *r;
      *sn = g1 / *r;
      for(i=0;i<count;i++)
	*r *= safemin2;
    } else {
      *r =  std::sqrt(f1*f1 + g1*g1);
      *cs = f1 / *r;
      *sn = g1 / *r;
    }
    if(std::abs(*f)>std::abs(*g) && *cs<0.0) {
      *cs *= -1.0;
      *sn *= -1.0;
      *r  *= -1.0;
    }
  }
  return;
}
      
