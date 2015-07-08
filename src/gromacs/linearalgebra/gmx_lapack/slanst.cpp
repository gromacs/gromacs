#include <cctype>
#include <cmath>

#include "../gmx_lapack.h"


float
F77_FUNC(slanst,SLANST)(const char *norm,
	int *n,
	float *d,
	float *e)
{
  const char ch=std::toupper(*norm);
  float dtemp,max,val,scale,sum;
  int i,j;


  if(*n<=0)
    return 0.0;
  
  switch(ch) {
  case 'M':
    max = std::abs(d[*n-1]);
      for(i=0;i<(*n-1);i++) {
	dtemp = std::abs(d[i]);
	if(dtemp>max)
	  max = dtemp;
	dtemp = std::abs(e[i]);
	if(dtemp>max)
	  max = dtemp;
      }
    val = max;
    break;
    
  case 'O':
  case '1':
  case 'I':

    if(*n==1)
      val = std::abs(d[0]);
    else {
      max = std::abs(d[0]) + std::abs(e[0]);
      dtemp = std::abs(e[*n-2]) + std::abs(d[*n-1]);
      if(dtemp>max)
	max = dtemp;
      for(i=1;i<(*n-1);i++) {
	dtemp = std::abs(d[i]) + std::abs(e[i]) + std::abs(e[i-1]);
	if(dtemp>max)
	  max = dtemp;
      }
      val = max;
    }
    break;

  case 'F':
  case 'E':
    scale = 0.0;
    sum   = 1.0;
    i = *n-1;
    j = 1;
    if(*n>1) {
      F77_FUNC(slassq,SLASSQ)(&i,e,&j,&scale,&sum);
      sum *= 2;
    }
    F77_FUNC(slassq,SLASSQ)(n,d,&j,&scale,&sum);
    val = scale *  std::sqrt(sum);
    break;
    
  default:
    val = 0.0;
    break;
  }
  return val;
}
