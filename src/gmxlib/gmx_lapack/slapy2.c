#include <math.h>
#include "gmx_lapack.h"

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

  if(z==0) 
    return w;
  else {
    z = z/w;
    return w*sqrt(1.0+z*z);
  }
}
  
