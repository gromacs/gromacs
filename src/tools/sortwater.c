#include <stdlib.h>
#include "typedefs.h"
#include "random.h"
#include "smalloc.h"
#include "vec.h"
#include "sortwater.h"

void randwater(int astart,int nwater,int nwatom,rvec x[],rvec v[],int *seed)
{
  int  i,j,wi,wj,*tab;
  rvec buf;
  
  snew(tab,nwater);
  for(i=0; (i<nwater); i++)
    tab[i]=i;
  for(j=0; (j<7*nwater); j++) {
    wi = (int) (nwater*rando(seed)) % nwater;
    do {
      wj = (int) (nwater*rando(seed)) % nwater;
    } while (wi == wj);
    wi = astart+wi*nwatom;
    wj = astart+wj*nwatom;
    /* Swap coords for wi and wj */
    for(i=0; (i<nwatom); i++) {
      copy_rvec(x[wi+i],buf);
      copy_rvec(x[wj+i],x[wi+i]);
      copy_rvec(buf,x[wj+i]);
      if (v) {
	copy_rvec(v[wi+i],buf);
	copy_rvec(v[wj+i],v[wi+i]);
	copy_rvec(buf,v[wj+i]);
      }
    }
  }
  sfree(tab);
}

static int rvcomp(const void *a,const void *b)
{
  rvec *aa,*bb;
  
  aa = (rvec *)a;
  bb = (rvec *)b;
  if (aa[XX] < bb[XX])
    return -1;
  else if (aa[XX] > bb[XX])
    return 1;
  else 
    return 0;
}

void sortwater(int nwater,int nwatom,rvec x[])
{
  qsort(x,nwater,nwatom*sizeof(x[0]),rvcomp);
}
