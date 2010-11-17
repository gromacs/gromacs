#include <math.h>
#include "gmx_lapack.h"

void
F77_FUNC(dlasdt,DLASDT)(int *n,
	int *lvl,
	int *nd,
	int *inode,
	int *ndiml,
	int *ndimr,
	int *msub)
{
  int maxn = (*n > 1) ? *n : 1;
  double temp;
  int i,il,ir,llst,nlvl,ncrnt;

  temp = log( ((double) maxn) / ((double)(*msub+1))) / log(2.0);
  
  *lvl = 1 + (int) temp;

  i = *n / 2;
  inode[0] = i + 1;
  ndiml[0] = i;
  ndimr[0] = *n - i - 1;
  il = -1;
  ir = 0;
  llst = 1;

  for(nlvl=1;nlvl<*lvl;nlvl++) {
    for(i=0;i<llst;i++) {
      il += 2;
      ir += 2;
      ncrnt = llst + i - 1;
      ndiml[il] = ndiml[ncrnt] / 2;
      ndimr[il] = ndiml[ncrnt] - ndiml[il] - 1;
      inode[il] = inode[ncrnt] - ndimr[il] - 1;
      ndiml[ir] = ndimr[ncrnt] / 2;
      ndimr[ir] = ndimr[ncrnt] - ndiml[ir] - 1;
      inode[ir] = inode[ncrnt] + ndiml[ir] + 1;
    }
    llst *= 2;
  }
  *nd = llst*2 - 1;
  return;
}
