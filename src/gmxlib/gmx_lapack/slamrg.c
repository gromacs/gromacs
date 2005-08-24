#include "gmx_lapack.h"


/* LAPACK */


void
F77_FUNC(slamrg,SLAMRG)(int *n1,
                        int *n2,
                        float *a,
                        int *dtrd1,
                        int *dtrd2,
                        int *index)
{
  int n1sv = *n1;
  int n2sv = *n2;
  int i,ind1,ind2;

  if(*dtrd1>0)
    ind1 = 0;
  else
    ind1 = *n1-1;

  if(*dtrd2>0)
    ind2 = *n1;
  else
    ind2 = *n1+*n2-1;

  i = 0;
  
  while(n1sv>0 && n2sv>0) {
    if(a[ind1]<=a[ind2]) {
      index[i] = ind1 + 1;
      i++;
      ind1 += *dtrd1;
      n1sv--;
    } else {
      index[i] = ind2 + 1;
      i++;
      ind2 += *dtrd2;
      n2sv--;
    }
  }

  if(n1sv==0) {
    for(n1sv=1;n1sv<=n2sv;n1sv++) {
      index[i] = ind2 + 1;
      i++;
      ind2 += *dtrd2;
    } 
  } else {
    for(n2sv=1;n2sv<=n1sv;n2sv++) {
      index[i] = ind1 + 1;
      i++;
      ind1 += *dtrd1;
    } 
  }
  return;
}
