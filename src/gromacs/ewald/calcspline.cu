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
#include "check.h"

#include "gromacs/utility/basedefinitions.h"
//#include "gromacs/utility/real.h"
typedef float real;
#include "gromacs/math/vectypes.h"

#include <cuda.h>
#include <vector>

#define PME_ORDER_MAX 12
typedef real *splinevec[DIM];
extern gpu_flags calcspline_gpu_flags;

/* Macro to force loop unrolling by fixing order.
 * This gives a significant performance gain.
 */
#define CALC_SPLINE(order, max_order)		   \
    {                                              \
        real dr, div;                               \
        real data[max_order];                  \
        real ddata[max_order];                 \
                                               \
        for (int j = 0; j < DIM; j++)                     \
        {                                          \
	    dr  = fractx[i*DIM + j];		   \
                                               \
            /* dr is relative offset from lower cell limit */ \
            data[order-1] = 0;                     \
            data[1]       = dr;                          \
            data[0]       = 1 - dr;                      \
                                               \
            _Pragma("unroll")  \
            for (int k = 3; k < order; k++)               \
            {                                      \
                div       = 1.0/(k - 1.0);               \
                data[k-1] = div*dr*data[k-2];      \
                _Pragma("unroll")  \
                for (int l = 1; l < (k-1); l++)           \
                {                                  \
                    data[k-l-1] = div*((dr+l)*data[k-l-2]+(k-l-dr)* \
                                       data[k-l-1]);                \
                }                                  \
                data[0] = div*(1-dr)*data[0];      \
            }                                      \
            /* differentiate */                    \
            ddata[0] = -data[0];                   \
            _Pragma("unroll")  \
            for (int k = 1; k < order; k++)               \
            {                                      \
                ddata[k] = data[k-1] - data[k];    \
            }                                      \
                                               \
            div           = 1.0/(order - 1);                 \
            data[order-1] = div*dr*data[order-2];  \
            _Pragma("unroll")  \
	    for (int l = 1; l < (order-1); l++)  \
            {                                      \
                data[order-l-1] = div*((dr+l)*data[order-l-2]+    \
                                       (order-l-dr)*data[order-l-1]); \
            }                                      \
            data[0] = div*(1 - dr)*data[0];        \
                                               \
            _Pragma("unroll")  \
            for (int k = 0; k < order; k++)                 \
            {                                      \
	        theta[j*order*nr + i*order+k]  = data[k];  \
                dtheta[j*order*nr + i*order+k] = ddata[k];	\
            }                                      \
        }                                          \
    }


#define MAKE_BSPLINES_KERNEL(order, max_order) \
  {  \
    int i = blockIdx.x * blockDim.x + threadIdx.x;  \
    if (i < nr && (bDoSplines || coefficient[i] != 0.0)) {  \
      CALC_SPLINE(order, max_order);  \
    }  \
  }


__global__ void make_bsplines_kernel_4(real *theta, real *dtheta,
				       real *fractx, int nr, real *coefficient,
				       bool bDoSplines) {
  MAKE_BSPLINES_KERNEL(4, 4);
}


__global__ void make_bsplines_kernel_5(real *theta, real *dtheta,
				       real *fractx, int nr, real *coefficient,
				       bool bDoSplines) {
  MAKE_BSPLINES_KERNEL(5, 5);
}


__global__ void make_bsplines_kernel_n(int order,
				       real *theta, real *dtheta,
				       real *fractx, int nr, real *coefficient,
				       bool bDoSplines) {
  MAKE_BSPLINES_KERNEL(order, PME_ORDER_MAX);
}

/* ######## 12 reg no cmem or lmem, vs 15 reg, 4 cmem, 16 lmem */

__global__ void make_bsplines_kernel_42(real *theta, real *dtheta,
					real *fractx, int nr, real *coefficient,
					bool bDoSplines) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < nr && (bDoSplines || coefficient[i] != 0.0)) {
    {
      real data[4];
      real ddata[4];

      for (int j = 0; j < DIM; j++)
        {
	  real dr  = fractx[i*DIM + j];
	  real dr2 = dr * dr;

	  /* dr is relative offset from lower cell limit */
	  data[3] = 0;
	  data[2] = .5*dr2;
	  data[1] = .5*(1-2*dr2+2*dr);
	  data[0] = .5*(1-dr)*(1-dr);

	  /* differentiate */
	  ddata[0] = -data[0];
	  ddata[1] = data[0] - data[1];
	  ddata[2] = data[1] - data[2];
	  ddata[3] = data[2] - data[3];

	  real div     = 1.0/3;
	  data[3] = div*dr*data[2];
	  data[2] = div*((dr+1)*data[1]+(3-dr)*data[2]);
	  data[1] = div*((dr+2)*data[0]+(2-dr)*data[1]);
	  data[0] = div*(1-dr)*data[0];

	  int t0 = j*nr+i;
	  theta[t0+0] = data[0];
	  theta[t0+1] = data[1];
	  theta[t0+2] = data[2];
	  theta[t0+3] = data[3];
	  dtheta[t0+0] = ddata[0];
	  dtheta[t0+1] = ddata[1];
	  dtheta[t0+2] = ddata[2];
	  dtheta[t0+3] = ddata[3];
        }
    }
  }
}

/* ######## */

#include "th_a.cuh"

void make_bsplines_gpu(splinevec theta, splinevec dtheta, int order,
		       rvec fractx[], int nr, int ind[], real coefficient[],
		       gmx_bool bDoSplines, int thread) {
  int size = nr * sizeof(real);
  int size_dim = DIM * size;
  int size_order = order * size;
  int size_order_dim = order * size_dim;

  real *theta_d = th_a(TH_ID_THETA, thread, size_order_dim, TH_LOC_CUDA);
  //cudaMalloc((void **) &theta_d, size_order_dim);
  real *dtheta_d = th_a(TH_ID_DTHETA, thread, size_order_dim, TH_LOC_CUDA);
  //cudaMalloc((void **) &dtheta_d, size_order_dim);

  real *fractx_d = th_a(TH_ID_FRACTX, thread, size_dim, TH_LOC_CUDA);
  //cudaMalloc((void **) &fractx_d, size_dim);
  real *fractx_h = th_a(TH_ID_FRACTX, thread, size_dim, TH_LOC_HOST);
  //real *fractx_h = (real *) malloc(size_dim);
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < DIM; ++j) {
      fractx_h[i * DIM + j] = fractx[ind[i]][j];
    }
  }
  cudaMemcpy(fractx_d, fractx_h, size_dim, cudaMemcpyHostToDevice);

  real *coefficient_d = th_a(TH_ID_COEFFICIENT, thread, size, TH_LOC_CUDA);
  //cudaMalloc((void **) &coefficient_d, size);
  real *coefficient_h = th_a(TH_ID_COEFFICIENT, thread, size, TH_LOC_HOST);
  //real *coefficient_h = (real *) malloc(size);
  for (int i = 0; i < nr; ++i) {
    coefficient_h[i] = coefficient[ind[i]];
  }
  cudaMemcpy(coefficient_d, coefficient_h, size, cudaMemcpyHostToDevice);

  int block_size = 32;
  int n_blocks = (nr + block_size - 1) / block_size;
  switch (order)
  {
  case 4: make_bsplines_kernel_4<<<n_blocks, block_size>>>
      (theta_d, dtheta_d, fractx_d, nr, coefficient_d, bDoSplines); break;
  case 5: make_bsplines_kernel_5<<<n_blocks, block_size>>>
      (theta_d, dtheta_d, fractx_d, nr, coefficient_d, bDoSplines); break;
  default: make_bsplines_kernel_n<<<n_blocks, block_size>>>
      (order,
       theta_d, dtheta_d, fractx_d, nr, coefficient_d, bDoSplines); break;
  }

  if (check_vs_cpu(calcspline_gpu_flags)) {
    for (int j = 0; j < DIM; ++j) {
      for (int i = 0; i < order; ++i) {
	real diff = theta[j][i];
	cudaMemcpy(&theta[j][i], theta_d + j * nr * order + i, sizeof(real), cudaMemcpyDeviceToHost);
	diff -= theta[j][i];
	fprintf(stderr, " %f", diff);
      }
      fprintf(stderr, "\n");
    }
  }


  for (int j = 0; j < DIM; ++j) {
    cudaMemcpy(theta[j], theta_d + j * nr * order, size_order, cudaMemcpyDeviceToHost);
  }
  for (int j = 0; j < DIM; ++j) {
    cudaMemcpy(dtheta[j], dtheta_d + j * nr * order, size_order, cudaMemcpyDeviceToHost);
  }

  /*
  cudaFree(theta_d);
  cudaFree(dtheta_d);
  cudaFree(fractx_d);
  cudaFree(coefficient_d);
  free(fractx_h);
  free(coefficient_h);
  */
}

void free_bsplines_gpu(int thread) {
  make_bsplines_gpu(NULL, NULL, 0, NULL, 0, NULL, NULL, false, thread);
}
