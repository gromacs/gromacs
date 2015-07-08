#include <cmath>
#include "../gmx_lapack.h"

#include "gromacs/utility/real.h"

double
F77_FUNC(dlapy2,DLAPY2)(double * x, double * y)
{
  double xabs,yabs;
  double w,z;

  xabs = std::abs(*x);
  yabs = std::abs(*y);
  
  if(xabs>yabs) {
    w = xabs;
    z = yabs;
  } else {
    w = yabs;
    z = xabs;
  }

  if( std::abs(z)<GMX_DOUBLE_MIN) 
    return w;
  else {
    z = z/w;
    return w* std::sqrt(1.0+z*z);
  }
}
  
