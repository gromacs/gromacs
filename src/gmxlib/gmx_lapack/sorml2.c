#include <ctype.h>
#include "gmx_lapack.h"

void
F77_FUNC(sorml2,SORML2)(char *side,
                        char *trans,
                        int *m,
                        int *n,
                        int *k,
                        float *a,
                        int *lda,
                        float *tau,
                        float *c,
                        int *ldc,
                        float *work,
                        int *info)
{
  char xside=toupper(*side);
  char xtrans=toupper(*trans);
  int i,i1,i2,i3,ni,mi,ic,jc;
  float aii;

  if(*m<=0 || *n<=0 || *k<=0)
    return;

  ic = jc = 0;

  if((xside=='L' && xtrans=='N') || (xside!='L' && xtrans!='N')) {
    i1 = 0;
    i2 = *k;
    i3 = 1;
  } else {
    i1 = *k-1;
    i2 = -1;
    i3 = -1;
  }
  
  if(xside=='L') {
    ni = *n;
    jc = 0;
  } else {
    mi = *m;
    ic = 0;
  }

  for(i=i1;i!=i2;i+=i3) {
    if(xside=='L') {
      mi = *m - i;
      ic = i;
    } else {
      ni = *n - i;
      jc = i;
    }
    aii = a[i*(*lda)+i];
    a[i*(*lda)+i] = 1.0;
    F77_FUNC(slarf,SLARF)(side,&mi,&ni,&(a[i*(*lda)+i]),lda,tau+i,
	   &(c[jc*(*ldc)+ic]),ldc,work);
    a[i*(*lda)+i] = aii;
  }
  return;
}
	     
