#include <math.h>
#include <ctype.h>
#include "gmx_lapack.h"


double
F77_FUNC(dlanst,DLANST)(char *norm,
	int *n,
	double *d,
	double *e)
{
  char ch=toupper(*norm);
  double dtemp,max,val,scale,sum;
  int i,j;


  if(*n<=0)
    return 0.0;
  
  switch(ch) {
  case 'M':
    max = fabs(d[*n-1]);
      for(i=0;i<(*n-1);i++) {
	dtemp = fabs(d[i]);
	if(dtemp>max)
	  max = dtemp;
	dtemp = fabs(e[i]);
	if(dtemp>max)
	  max = dtemp;
      }
    val = max;
    break;
    
  case 'O':
  case '1':
  case 'I':

    if(*n==1)
      val = fabs(d[0]);
    else {
      max = fabs(d[0]) + fabs(e[0]);
      dtemp = fabs(e[*n-2]) + fabs(d[*n-1]);
      if(dtemp>max)
	max = dtemp;
      for(i=1;i<(*n-1);i++) {
	dtemp = fabs(d[i]) + fabs(e[i]) + fabs(e[i-1]);
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
      F77_FUNC(dlassq,DLASSQ)(&i,e,&j,&scale,&sum);
      sum *= 2;
    }
    F77_FUNC(dlassq,DLASSQ)(n,d,&j,&scale,&sum);
    val = scale * sqrt(sum);
    break;
    
  default:
    val = 0.0;
    break;
  }
  return val;
}
