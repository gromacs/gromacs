#include "vec.h"
#include "typedefs.h"
#include "smalloc.h"
#include "fatal.h"
#include "calcgrid.h"

#define facNR 6
int factor[facNR] = {2,3,5,7,11,13};
int decomp[facNR];
int ng,ng_max,*list,nlist,nlist_alloc;

static void make_list(int start_fac)
{
  int i;
  
  if (ng < ng_max) {
    if (nlist >= nlist_alloc) {
      nlist_alloc += 100;
      srenew(list,nlist_alloc);
    }
    list[nlist] = ng;
    nlist++;

    for(i=start_fac; i<facNR; i++) {
      if (i<4 || (factor[4]+factor[5]==0)) {
	ng*=factor[i];
	decomp[i]++;
	make_list(i);
	ng/=factor[i];
	decomp[i]--;
      }
    }
  }
}

static int list_comp(const void *a,const void *b)
{
  return (list[*((int *)a)] - list[*((int *)b)]);
}

void calc_grid(matrix box,real gr_sp,int *nx,int *ny,int *nz,int nprocs)
{
  int  d,n[DIM];
  int  i,j=0,x,y,z,nmin[DIM],*list_index;
  rvec box_size,spac;
  real spac_av;
  bool bDone;
  
  if (gr_sp <= 0)
    fatal_error(0,"invalid fourier grid spacing: %g",gr_sp);
  
  for(d=0; d<DIM; d++)
    box_size[d] = box[d][d];

  n[XX] = *nx;
  n[YY] = *ny;
  n[ZZ] = *nz;

  ng = 1;
  ng_max = 1;
  for(d=0; d<DIM; d++) {
    nmin[d] = (int)(box_size[d]/gr_sp + 0.999);
    if (2*nmin[d] > ng_max)
      ng_max = 2*nmin[d];
  }
  nlist=0;
  nlist_alloc=0;
  list=NULL;
  make_list(0);

  if ((*nx<=0) || (*ny<=0) || (*nz<=0))
    fprintf(stderr,"Calculating fourier grid dimensions for%s%s%s\n",
	    *nx > 0 ? "":" X",*ny > 0 ? "":" Y",*nz > 0 ? "":" Z");

  snew(list_index,nlist);
  for(i=0; i<nlist; i++)
    list_index[i] = i;
  qsort(list_index,nlist,sizeof(list_index[0]),list_comp);
  if (debug)
    for(i=0; i<nlist; i++)
      fprintf(debug,"grid: %d\n",list[list_index[i]]);
  
  if (((*nx>0) && (*nx != nprocs*(*nx/nprocs))) ||
      ((*ny>0) && (*ny != nprocs*(*ny/nprocs))))
    fatal_error(0,"the x or y grid spacing (nx %d, ny %d) is not divisible by the number of processors (%d)",*nx,*ny,nprocs);
  
  for(d=0; d<DIM; d++) {
    if (n[d] <= 0) {
      bDone = FALSE;
      for(i=0; (i<nlist) && !bDone; i++) {
	j = list_index[i];
	bDone = ((list[j] >= nmin[d]) && 
		 ((d == ZZ) || (list[j] == nprocs*(list[j]/nprocs))));
      }
      if (!bDone)
	fatal_error(0,"could not find a grid spacing with nx and ny divisible by the number of processors (%d)",nprocs);
      n[d] = list[j];
    }
  }
  
  *nx = n[XX];
  *ny = n[YY];
  *nz = n[ZZ];
  fprintf(stderr,"Using a fourier grid of %dx%dx%d, spacing %.3f %.3f %.3f\n",
	  *nx,*ny,*nz,
	  box_size[XX]/(*nx),box_size[YY]/(*ny),box_size[ZZ]/(*nz));

  sfree(list_index);
}


