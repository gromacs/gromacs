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
  for(j=0; (j<23*nwater); j++) {
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

static rvec *xptr;
static int  nwat;

static int rvcomp(const void *a,const void *b)
{
  int aa,bb;
  
  aa = nwat*(*((int *)a));
  bb = nwat*(*((int *)b));
  if (xptr[aa][XX] < xptr[bb][XX])
    return -1;
  else if (xptr[aa][XX] > xptr[bb][XX])
    return 1;
  else 
    return 0;
}

void sortwater(int astart,int nwater,int nwatom,rvec x[],rvec v[])
{
  int  i,j,i0,size,rvi;
  int  *rvindex;
  rvec *tmp;
  
  /* Sort indices to rvecs */
  snew(rvindex,nwater);
  for(i=0; (i<nwater); i++) 
    rvindex[i] = i;
  xptr = x+astart;
  nwat = nwatom;
  
  qsort(rvindex,nwater,sizeof(rvindex[0]),rvcomp);
  if (debug)
    for(i=0; (i<nwater); i++) {
      rvi = rvindex[i]*nwatom;
      fprintf(debug,"rvindex[%5d] = %5d (x = %8.3f  %8.3f  %8.3f)\n",
	      i,rvi,x[astart+rvi][XX],x[astart+rvi][YY],x[astart+rvi][ZZ]);
    }
  snew(tmp,nwater*nwatom);
  
  for(i=0; (i<nwater); i++) {
    i0 = astart+nwatom*rvindex[i];
    for(j=0; (j<nwatom); j++) 
      copy_rvec(x[i0+j],tmp[nwatom*i+j]);
  }
  for(i=0; (i<nwater*nwatom); i++)
    copy_rvec(tmp[i],x[astart+i]);
    
  for(i=0; (i<nwater); i++) {
    i0 = astart+nwatom*rvindex[i];
    for(j=0; (j<nwatom); j++) 
      copy_rvec(v[i0+j],tmp[nwatom*i+j]);
  }
  for(i=0; (i<nwater*nwatom); i++)
    copy_rvec(tmp[i],v[astart+i]);
    
  sfree(tmp);
  sfree(rvindex);
}

