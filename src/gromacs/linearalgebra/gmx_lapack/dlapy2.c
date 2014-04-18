#include <math.h>
#include "../gmx_lapack.h"

#include "gromacs/utility/real.h"

double
F77_FUNC(dlapy2,DLAPY2)(double * x, double * y)
{
  double xabs,yabs;
  double w,z;

  xabs = fabs(*x);
  yabs = fabs(*y);
  
  if(xabs>yabs) {
    w = xabs;
    z = yabs;
  } else {
    w = yabs;
    z = xabs;
  }

  if( fabs(z)<GMX_DOUBLE_MIN) 
    return w;
  else {
    z = z/w;
    return w*sqrt(1.0+z*z);
  }
}
  
