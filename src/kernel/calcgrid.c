#include "vec.h"
#include "typedefs.h"
#include "smalloc.h"
#include "fatal.h"
#include "calcgrid.h"

#define facNR 6
int factor[facNR] = {2,3,5,7,11,13};
int decomp[facNR];
int ng,ng_max,*list,n_list,n_list_alloc;

static void make_list(int start_fac)
{
  int i;
  
  if (ng < ng_max) {
    if (n_list >= n_list_alloc) {
      n_list_alloc += 100;
      srenew(list,n_list_alloc);
    }
    list[n_list] = ng;
    n_list++;

    for(i=start_fac; i<facNR; i++) {
      /* allow any power of 2, 3, 5 and 7, but only one of 11 or 13 */
      if (i<4 || (decomp[4]+decomp[5]==0)) {
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
  return (*((int *)a) - *((int *)b));
}

real calc_grid(matrix box,real gr_sp,int *nx,int *ny,int *nz,int nprocs)
{
  int  d,n[DIM];
  int  i,nmin[DIM];
  rvec box_size,spacing;
  real max_spacing;
  
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
  n_list=0;
  n_list_alloc=0;
  list=NULL;
  for(i=0; i<facNR; i++)
    decomp[i]=0;
  make_list(0);

  if ((*nx<=0) || (*ny<=0) || (*nz<=0))
    fprintf(stderr,"Calculating fourier grid dimensions for%s%s%s\n",
	    *nx > 0 ? "":" X",*ny > 0 ? "":" Y",*nz > 0 ? "":" Z");

  qsort(list,n_list,sizeof(list[0]),list_comp);
  if (debug)
    for(i=0; i<n_list; i++)
      fprintf(debug,"grid: %d\n",list[i]);
  
  if (((*nx>0) && (*nx != nprocs*(*nx/nprocs))) ||
      ((*ny>0) && (*ny != nprocs*(*ny/nprocs))))
    fatal_error(0,"the x or y grid spacing (nx %d, ny %d) is not divisible by the number of processors (%d)",*nx,*ny,nprocs);
  
  for(d=0; d<DIM; d++) {
    for(i=0; (i<n_list) && (n[d]<=0); i++)
      if ((list[i] >= nmin[d]) && 
	  ((d == ZZ) || (list[i] == nprocs*(list[i]/nprocs))))
	n[d] = list[i];
    if (n[d] <= 0)
      fatal_error(0 ,"could not find a grid spacing with nx and ny divisible by the number of processors (%d)",nprocs);
  }
  
  max_spacing = 0;
  for(d=0; d<DIM; d++) {
    spacing[d] = box_size[d]/n[d];
    if (spacing[d] > max_spacing)
      max_spacing = spacing[d];
  }
  *nx = n[XX];
  *ny = n[YY];
  *nz = n[ZZ];
  fprintf(stderr,"Using a fourier grid of %dx%dx%d, spacing %.3f %.3f %.3f\n",
	  *nx,*ny,*nz,spacing[XX],spacing[YY],spacing[ZZ]);

  return max_spacing;
}


