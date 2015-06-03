/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013-2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#include "pme.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/math/vectypes.h"
#include "check.h"

#include <cuda.h>

#define PME_ORDER_MAX 12
typedef real *splinevec[DIM];

extern gpu_flags spread_gpu_flags;
extern gpu_flags spread_bunching_gpu_flags;

#include "thread_mpi/mutex.h"

/* This has to be a macro to enable full compiler optimization with xlC (and probably others too) */
#define DO_BSPLINE(order)                            \
     for (int ithx = 0; ithx < order; ithx++)	     \
     {  \
        int index_x = (i0[i]+ithx)*pny*pnz;                     \
        real valx    = coefficient[i]*thx[i*order+ithx];                          \
                                                     \
        {							      \
            real valxy    = valx*thy[i*order+ithy];                   \
            int index_xy = index_x+(j0[i]+ithy)*pnz;            \
                                                     \
            {					     \
                int index_xyz        = index_xy+(k0[i]+ithz);   \
                /*grid[index_xyz] += valxy*thz[i*order+ithz];*/      \
		atomicAdd(&grid[index_xyz], valxy*thz[i*order+ithz]); \
	    }						 \
	}						 \
    }


#define SPREAD_COEFFICIENTS_KERNEL(order) \
  {  \
    int i = blockIdx.z * blockDim.z + threadIdx.z;  \
    if (i < n) {  \
      DO_BSPLINE(order);			     \
    }  \
  }

__global__ void spread_coefficients_kernel_4(int n,
					     real *grid,
					     int *i0, int *j0, int *k0,
					     int pny, int pnz,
					     real *coefficient,
					     real *thx, real *thy, real *thz)
{
  int ithz = threadIdx.x;
  int ithy = threadIdx.y;
  SPREAD_COEFFICIENTS_KERNEL(4);
}

__global__ void spread_coefficients_kernel_5(int n,
					     real *grid,
					     int *i0, int *j0, int *k0,
					     int pny, int pnz,
					     real *coefficient,
					     real *thx, real *thy, real *thz)
{
  int ithz = threadIdx.x;
  int ithy = threadIdx.y;
  SPREAD_COEFFICIENTS_KERNEL(5);
}

__global__ void spread_coefficients_kernel_n(int order,
					     int n,
					     real *grid,
					     int *i0, int *j0, int *k0,
					     int pny, int pnz,
					     real *coefficient,
					     real *thx, real *thy, real *thz)
{
  for (int ithy = 0; (ithy < order); ithy++) {
    for (int ithz = 0; (ithz < order); ithz++) {
      SPREAD_COEFFICIENTS_KERNEL(order);
    }
  }
}

#include "th_a.cuh"

static tMPI::mutex print_mutex;

void try_bunching(int pnx, int pny, int pnz,
		  ivec *atc_idx, int *spline_ind, int spline_n,
		  real *atc_coefficient, int atc_n_foo) {
  int n = atc_n_foo; // spline_n;

  const int bnx = 32, bny = 32, bnz = 1;
  int bn = 0;
  int mark[bnx][bny];
  int marko[bnx][bny];
  for (int x = 0; x < bnx; ++x) {
    for (int y = 0; y < bny; ++y) {
      mark[x][y] = 0;
      marko[x][y] = 0;
    }
  }
  int maxmark = 0;
  int order = 4, bunch = 1; //(int) (order * 1.25 + .5);

  print_mutex.lock();
  fprintf(stderr, "%d %dx%dx%d\n", n, pnx, pny, pnz);
  for (int i = 0; i < atc_n_foo; ++i) {
    //for (int ii = 0; ii < n; ++ii) {
    //int i = spline_ind[ii];
    real coefficient_i = atc_coefficient[i];
    if (coefficient_i == 0) {
      continue;
    }
    int *idxptr = atc_idx[i];
    int x = idxptr[XX], y = idxptr[YY], z = idxptr[ZZ];

    if (!(x < bnx && y < bny && z < bnz)) {
	continue;
    }

    if (++mark[x / bunch][y / bunch] > maxmark) {
      maxmark = mark[x / bunch][y / bunch];
    }
    for (int dxb = 0; dxb < order / bunch; ++dxb) {
      for (int dyb = 0; dyb < order / bunch; ++dyb) {
	if (x / bunch + dxb < bnx && y / bunch + dyb < bny) {
	  ++marko[x / bunch + dxb][y / bunch + dyb];
	}
      }
    }

    ++bn;

    //fprintf(stderr, "(%d,%d,%d) %f\n", x, y, z, coefficient_i);
  }
    fprintf(stderr, "bn %d maxmark %d\n", bn, maxmark);
    int avgsum = 0, avgcount = 0;
    for (int x = 0; x < bnx / bunch; ++x) {
      for (int y = 0; y < bny / bunch; ++y) {
	fprintf(stderr, " %d", mark[x][y]);
	avgsum += mark[x][y]; ++avgcount;
      }
      fprintf(stderr, " mark avg %d\n", avgsum / avgcount);
    }
  print_mutex.unlock();
}

void spread_coefficients_bsplines_thread_gpu_2
(int pnx, int pny, int pnz, int offx, int offy, int offz,
 real *grid, int order, ivec *atc_idx, int *spline_ind, int spline_n,
 real *atc_coefficient, splinevec *spline_theta, int atc_n_foo,
 int thread)
{
  //fprintf(stderr, "Hello spread! %d %d\n", thread, spline_n);

    int ndatatot = pnx*pny*pnz;
    int size_grid = ndatatot * sizeof(real);

    real *grid_check;
    if (check_vs_cpu(spread_gpu_flags)) {
      grid_check = th_a(TH_ID_GRID, thread, size_grid, TH_LOC_HOST);
      memcpy(grid_check, grid, ndatatot * sizeof(real));
    }

    for (int i = 0; i < ndatatot; i++)
    {
      // FIX clear grid on device instead
        grid[i] = 0;
    }

    real *grid_d = th_a(TH_ID_GRID, thread, size_grid, TH_LOC_CUDA);
    cudaMemcpy(grid_d, grid, size_grid, cudaMemcpyHostToDevice);

    int size_real = spline_n * sizeof(real);
    int size_int = spline_n * sizeof(int);
    int *i0 = th_i(TH_ID_I0, thread, size_int, TH_LOC_HOST);
    int *j0 = th_i(TH_ID_J0, thread, size_int, TH_LOC_HOST);
    int *k0 = th_i(TH_ID_K0, thread, size_int, TH_LOC_HOST);
    real *coefficient = th_a(TH_ID_COEFFICIENT, thread, size_real, TH_LOC_HOST);
    real *thx = th_a(TH_ID_THX, thread, size_real * order, TH_LOC_HOST);
    real *thy = th_a(TH_ID_THY, thread, size_real * order, TH_LOC_HOST);
    real *thz = th_a(TH_ID_THZ, thread, size_real * order, TH_LOC_HOST);

    int *i0_d = th_i(TH_ID_I0, thread, size_int, TH_LOC_CUDA);
    int *j0_d = th_i(TH_ID_J0, thread, size_int, TH_LOC_CUDA);
    int *k0_d = th_i(TH_ID_K0, thread, size_int, TH_LOC_CUDA);
    real *coefficient_d = th_a(TH_ID_COEFFICIENT, thread, size_real, TH_LOC_CUDA);
    real *thx_d = th_a(TH_ID_THX, thread, size_real * order, TH_LOC_CUDA);
    real *thy_d = th_a(TH_ID_THY, thread, size_real * order, TH_LOC_CUDA);
    real *thz_d = th_a(TH_ID_THZ, thread, size_real * order, TH_LOC_CUDA);

    int oo = 0;

    if (check_vs_cpu(spread_bunching_gpu_flags)) {
      try_bunching(pnx, pny, pnz,
		   atc_idx, spline_ind, spline_n, atc_coefficient, atc_n_foo);
    }

    for (int ii = 0; ii < spline_n; ii++)
    {
        int i           = spline_ind[ii];
        real coefficient_i = atc_coefficient[i];
	if (coefficient_i == 0) {
	   continue;
	}

	coefficient[oo] = coefficient_i;

	int *idxptr = atc_idx[i];
	int iiorder = ii*order;
	int ooorder = oo*order;

	i0[oo]   = idxptr[XX] - offx;
	j0[oo]   = idxptr[YY] - offy;
	k0[oo]   = idxptr[ZZ] - offz;

	for (int o = 0; o < order; ++o) {
	  thx[ooorder + o] = (*spline_theta)[XX][iiorder + o];
	  thy[ooorder + o] = (*spline_theta)[YY][iiorder + o];
	  thz[ooorder + o] = (*spline_theta)[ZZ][iiorder + o];
	}
	++oo;
    }

    int n = oo;

    //fprintf(stderr, "World! %d %d/%d\n", thread, n, spline_n);

    cudaMemcpy(i0_d, i0, size_int, cudaMemcpyHostToDevice);
    cudaMemcpy(j0_d, j0, size_int, cudaMemcpyHostToDevice);
    cudaMemcpy(k0_d, k0, size_int, cudaMemcpyHostToDevice);
    cudaMemcpy(coefficient_d, coefficient, size_real, cudaMemcpyHostToDevice);
    cudaMemcpy(thx_d, thx, size_real * order, cudaMemcpyHostToDevice);
    cudaMemcpy(thy_d, thy, size_real * order, cudaMemcpyHostToDevice);
    cudaMemcpy(thz_d, thz, size_real * order, cudaMemcpyHostToDevice);

  int block_size = 32;
  int n_blocks = (n + block_size - 1) / block_size;
  dim3 dimGrid(1, 1, n_blocks);
  dim3 dimBlockOrder(order, order, block_size);
  dim3 dimBlockOne(1, 1, block_size);
    switch (order)
    {
    case 4: spread_coefficients_kernel_4<<<dimGrid, dimBlockOrder>>>
	(n, grid_d, i0_d, j0_d, k0_d, pny, pnz,
	 coefficient_d, thx_d, thy_d, thz_d); break;
    case 5: spread_coefficients_kernel_5<<<dimGrid, dimBlockOrder>>>
	(n, grid_d, i0_d, j0_d, k0_d, pny, pnz,
	 coefficient_d, thx_d, thy_d, thz_d); break;
    default: spread_coefficients_kernel_n<<<dimGrid, dimBlockOne>>>
	(order,
	 n, grid_d, i0_d, j0_d, k0_d, pny, pnz,
	 coefficient_d, thx_d, thy_d, thz_d); break;
    }

  if (check_vs_cpu(spread_gpu_flags)) {
    print_mutex.lock();
    fprintf(stderr, "Check %d  (%d x %d x %d)\n",
	    thread, pnx, pny, pnz);
    for (int i = 0; i < ndatatot; ++i) {
      real diff = grid_check[i];
      real cpu_v = grid_check[i];
      cudaMemcpy(&grid_check[i], &grid_d[i], sizeof(real), cudaMemcpyDeviceToHost);
      diff -= grid_check[i];
      real gpu_v = grid_check[i];
      if (diff != 0) {
	real absdiff = fabs(diff) / fabs(cpu_v);
	if (absdiff > .000001) {
	  fprintf(stderr, "%dppm", (int) (absdiff * 1e6));
	  if (absdiff > .0001) {
	    fprintf(stderr, " value %f ", cpu_v);
	  }
	} else {
	  fprintf(stderr, "~");
	}
	//fprintf(stderr, "(%f - %f)", cpu_v, gpu_v);
      } else {
	if (gpu_v == 0) {
	  fprintf(stderr, "0");
	} else {
	  fprintf(stderr, "=");
	}
      }
      if ((i + 1) % pnz == 0) {
	fprintf(stderr, "\n");
      }
    }
    print_mutex.unlock();
  }
  cudaMemcpy(grid, grid_d, size_grid, cudaMemcpyDeviceToHost);
}
