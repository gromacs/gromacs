#include <math.h>
#include "gmx_lapack.h"

#include <types/simple.h>

float
F77_FUNC(slapy2,SLAPY2)(float * x, float * y)
{
  float xabs,yabs;
  float w,z;

  xabs = fabs(*x);
  yabs = fabs(*y);
  
  if(xabs>yabs) {
    w = xabs;
    z = yabs;
  } else {
    w = yabs;
    z = xabs;
  }

  if( fabs(z)<GMX_FLOAT_MIN) 
    return w;
  else {
    z = z/w;
    return w*sqrt(1.0+z*z);
  }
}
  
