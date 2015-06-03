#include "pme.h"
#include "pme-internal.h"

#include "gromacs/fft/fft.h"
#include "gromacs/fft/parallel_3dfft.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"
#include "gromacs/math/vectypes.h"

#include "check.h"

#include <cuda.h>
#include <cufft.h>

extern gpu_flags fft_gpu_flags;

struct gmx_parallel_3dfft_gpu {
  real *real_data;
  t_complex *complex_data;

  /* unused */
  MPI_Comm                  comm[2];
  gmx_bool                  bReproducible;
  int                       nthreads;

  ivec                      complex_order;
  ivec                      local_ndata;
  ivec                      local_offset;
  ivec                      local_size;

  int n[3];
  cufftHandle plan;
  cufftReal *rdata;
  cufftComplex *cdata;
};


struct test_3dfft_t {
  gmx_parallel_3dfft_t pfft_setup;
  gmx_parallel_3dfft_gpu_t pfft_setup_gpu;
};
test_3dfft_t *test_3dfft_init(int n0, int n1, int n2) {
  test_3dfft_t *t = new test_3dfft_t();
  ivec ndata;
  ndata[0] = n0;
  ndata[1] = n1;
  ndata[2] = n2;
  real *real_data;
  t_complex *complex_data;
  MPI_Comm comm[2];
  comm[0] = MPI_COMM_NULL;
  comm[1] = MPI_COMM_NULL;
  gmx_bool reprod = false;
  int nthreads = 1;
  fprintf(stderr, "test_3dfft_init 1\n");
  gmx_parallel_3dfft_init(&t->pfft_setup, ndata, &real_data, &complex_data,
			  comm, reprod, nthreads);
  fprintf(stderr, "test_3dfft_init 2\n");
  gmx_parallel_3dfft_init_gpu(&t->pfft_setup_gpu, ndata, &real_data, &complex_data,
			      comm, reprod, nthreads);
  fprintf(stderr, "test_3dfft_init 3\n");
  return t;
}
int test_3dfft_size(test_3dfft_t *t, int dim) {
  ivec local_ndata, local_offset, local_size;
  gmx_parallel_3dfft_real_limits(t->pfft_setup,
				 local_ndata, local_offset, local_size);
  return local_size[dim];
}
int test_3dfft_csize(test_3dfft_t *t, int dim) {
  ivec complex_order;
  ivec local_ndata, local_offset, local_size;
  gmx_parallel_3dfft_complex_limits(t->pfft_setup, complex_order,
				    local_ndata, local_offset, local_size);
  return local_size[dim];
}
real *test_3dfft_real(test_3dfft_t *t) {
  return t->pfft_setup_gpu->real_data;
}
t_complex *test_3dfft_complex(test_3dfft_t *t) {
  return t->pfft_setup_gpu->complex_data;
}
void test_3dfft_execute_cpu(test_3dfft_t *t, gmx_bool inverse) {
  gmx_fft_direction dir = inverse ? GMX_FFT_COMPLEX_TO_REAL : GMX_FFT_REAL_TO_COMPLEX;
  int thread = 0;
  gmx_wallcycle_t         wcycle = NULL;
  gmx_parallel_3dfft_execute(t->pfft_setup, dir, thread, wcycle);
}
void test_3dfft_execute_gpu(test_3dfft_t *t, gmx_bool inverse) {
  gmx_fft_direction dir = inverse ? GMX_FFT_COMPLEX_TO_REAL : GMX_FFT_REAL_TO_COMPLEX;
  int thread = 0;
  gmx_wallcycle_t         wcycle = NULL;
  gmx_parallel_3dfft_execute_gpu(t->pfft_setup_gpu, dir, thread, wcycle);
}
void test_3dfft_destroy(test_3dfft_t *t) {
  gmx_parallel_3dfft_destroy(t->pfft_setup);
  gmx_parallel_3dfft_destroy_gpu(t->pfft_setup_gpu);
  delete t;
}


int
    gmx_parallel_3dfft_init_gpu   (gmx_parallel_3dfft_gpu_t *    pfft_setup,
                               ivec                      ndata,
                               real **real_data,
                               t_complex **complex_data,
                               MPI_Comm                  comm[2],
                               gmx_bool                  bReproducible,
                               int                       nthreads)
{
  gmx_parallel_3dfft_gpu_t setup = new gmx_parallel_3dfft_gpu();

  // FIXME: this copies the already setup pointer, to check them after execute
  setup->real_data = *real_data;
  setup->complex_data = *complex_data;

  setup->comm[0] = comm[0];
  setup->comm[1] = comm[1];
  setup->bReproducible = bReproducible;
  setup->nthreads = nthreads;

  /*
    // (local pme and fft differs only by overlap (and pme > fft))
    pmeidx = ix*(local_pme[YY]*local_pme[ZZ])+iy*(local_pme[ZZ])+iz;
    fftidx = ix*(local_fft[YY]*local_fft[ZZ])+iy*(local_fft[ZZ])+iz;
    fftgrid[fftidx] = pmegrid[pmeidx];
    // TODO: align cufft minor dim to 128 bytes
   */
  setup->n[0] = ndata[0];
  setup->n[1] = ndata[1];
  setup->n[2] = ndata[2];
  int x = setup->n[0], y = setup->n[1], z = setup->n[2];

  cudaMalloc((void **) &setup->rdata, x * y * z * sizeof(cufftReal));
  cudaMalloc((void **) &setup->cdata, x * y * (z/2+1) * 2 * sizeof(cufftComplex));

  *pfft_setup = setup;
  //fprintf(stderr, "3dfft_init_gpu\n");

  // FIX: double plans?
  int rank = 3, batch = 1;
  if (cufftPlanMany(&setup->plan, rank, setup->n,
		    NULL, 0, 0,
		    NULL, 0, 0,
		    CUFFT_R2C,
		    batch)
      != CUFFT_SUCCESS) {
    fprintf(stderr, "PLAN_MANY FAIL!!!\n");
    setup = NULL; // FIX
  }


  return 0;
}

int
gmx_parallel_3dfft_real_limits_gpu(gmx_parallel_3dfft_gpu_t      pfft_setup,
                               ivec                      local_ndata,
                               ivec                      local_offset,
                               ivec                      local_size)
{
  //fprintf(stderr, "3dfft_real_limits_gpu\n");
  gmx_parallel_3dfft_gpu_t setup = pfft_setup;
  setup->local_ndata[0] = local_ndata[0];
  setup->local_ndata[1] = local_ndata[1];
  setup->local_ndata[2] = local_ndata[2];
  setup->local_offset[0] = local_offset[0];
  setup->local_offset[1] = local_offset[1];
  setup->local_offset[2] = local_offset[2];
  setup->local_size[0] = local_size[0];
  setup->local_size[1] = local_size[1];
  setup->local_size[2] = local_size[2];
  return 0;
}


int
gmx_parallel_3dfft_complex_limits_gpu(gmx_parallel_3dfft_gpu_t      pfft_setup,
                                  ivec                      complex_order,
                                  ivec                      local_ndata,
                                  ivec                      local_offset,
                                  ivec                      local_size)
{
  //fprintf(stderr, "3dfft_complex_limits_gpu\n");
  gmx_parallel_3dfft_gpu_t setup = pfft_setup;
  setup->complex_order[0] = complex_order[0];
  setup->complex_order[1] = complex_order[1];
  setup->complex_order[2] = complex_order[2];
  setup->local_ndata[0] = local_ndata[0];
  setup->local_ndata[1] = local_ndata[1];
  setup->local_ndata[2] = local_ndata[2];
  setup->local_offset[0] = local_offset[0];
  setup->local_offset[1] = local_offset[1];
  setup->local_offset[2] = local_offset[2];
  setup->local_size[0] = local_size[0];
  setup->local_size[1] = local_size[1];
  setup->local_size[2] = local_size[2];
  return 0;
}

__global__ void transpose_xyz_yzx_kernel(int nx, int ny, int nz,
					 cufftComplex *cdata,
					 bool forward)
{
  // transpose cdata to be contiguous in a y z x loop
  // z-dim has nz/2+1 elems
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z;
  if (x < nx && y < ny && z < nz) {
    int idx1 = (x * ny + y) * (nz/2+1) + z;
    int idx2 = ((ny + y) * (nz/2+1) + z) * nx + x;
    if (forward) {
      cdata[idx2] = cdata[idx1];
    } else {
      cdata[idx1] = cdata[idx2];
    }
  }
}

void transpose_xyz_yzx(int nx, int ny, int nz,
		       cufftComplex *cdata,
		       bool forward)
{
  int block_size = 32;
  dim3 dimGrid((nx + block_size - 1) / block_size, ny, nz/2+1);
  dim3 dimBlock(block_size, 1, 1);
  transpose_xyz_yzx_kernel<<<dimGrid, dimBlock>>>(nx, ny, nz, cdata, forward);
  //cudaMemcpy(cdata, cdata + nx * ny * (nz/2+1), nx * ny * (nz/2+1) * sizeof(cufftComplex), cudaMemcpyDeviceToDevice);
}

int
gmx_parallel_3dfft_execute_gpu(gmx_parallel_3dfft_gpu_t    pfft_setup,
                           enum gmx_fft_direction  dir,
                           int                     thread,
                           gmx_wallcycle_t         wcycle)
{
  //fprintf(stderr, "3dfft_execute_gpu\n");
  gmx_parallel_3dfft_gpu_t setup = pfft_setup;

  int x = setup->n[0], y = setup->n[1], z = setup->n[2];

  //int rank = 3, batch = 1;
  /*fprintf(stderr, "FFT plan %dx%dx%d %s %p(%d)->%p(%d)\n", x, y, z,
	  dir == GMX_FFT_REAL_TO_COMPLEX ? "CUFFT_R2C" : "CUFFT_C2R",
	  setup->real_data, (int) sizeof(real),
	  setup->complex_data, (int) sizeof(t_complex));*/
  /* FIX plan in advance (in init)
  if (cufftPlanMany(&setup->plan, rank, setup->n,
		    NULL, 0, 0,
		    NULL, 0, 0,
		    dir == GMX_FFT_REAL_TO_COMPLEX ? CUFFT_R2C : CUFFT_C2R,
		    batch)
      != CUFFT_SUCCESS) {
    fprintf(stderr, "PLAN_MANY FAIL!!! %d %p\n", thread, &wcycle);
    setup = NULL; // FIX
  }
  */

  if (dir == GMX_FFT_REAL_TO_COMPLEX) {
    cudaMemcpy(setup->rdata, setup->real_data, x * y * z * sizeof(real), cudaMemcpyHostToDevice);
    cufftExecR2C(setup->plan, setup->rdata, setup->cdata);
    // FIXME: -> y major, z middle, x minor or continuous
    transpose_xyz_yzx(x, y, z, setup->cdata, true);
  } else {
    cudaMemcpy(setup->cdata + x * y * (z/2+1), setup->complex_data, x * y * (z/2+1) * sizeof(t_complex), cudaMemcpyHostToDevice);
    // FIXME: y major, z middle, x minor or continuous ->
    transpose_xyz_yzx(x, y, z, setup->cdata, false);
    cufftExecC2R(setup->plan, setup->cdata, setup->rdata);
  }

  if (check_vs_cpu(fft_gpu_flags)) {
    //print_lock();
    for (int ix = 0; ix < x; ++ix) {
      fprintf(stderr, "plane %d\n", ix);
      for (int iy = 0; iy < y; ++iy) {
	for (int iz = 0; iz < z; ++iz) {
	  int i = (ix * y + iy) * z + iz;
	  fprintf(stderr, " %5.2f", setup->real_data[i]);
	}
	fprintf(stderr, "\n");
      }
    }
#if 0
    real sum = 0;
    for (int i = 0; i < x * y * z; ++i) {
      sum += setup->real_data[i];
    }
    check_real("FFT r all", setup->rdata, setup->real_data, x * y * z, true);
    fprintf(stderr, "sum %f\n", (double) sum);
    for (int i = 0; i < x * y * (z/2 + 1); ++i) {
      //check_real("FFT r", &setup->rdata[i], &setup->real_data[i], 1, true);

      cufftComplex c;
      cudaMemcpy(&c, &setup->cdata[i + x * y * (z/2+1)], sizeof(cufftComplex), cudaMemcpyDeviceToHost);
      fprintf(stderr, "FFT %d %f %f\t", i, (double) c.x, (double) c.y);
      fprintf(stderr, "FFT vs %f %f\n", (double) setup->complex_data[i].re, (double) setup->complex_data[i].im);
      /*
      fprintf(stderr, "FFT %d %p %p\n", i, &setup->cdata[i].x, &setup->cdata[i].y);
      check_real("FFT cr", &setup->cdata[i].x, &setup->complex_data[i].re, 1, true);
      check_real("FFT ci", &setup->cdata[i].y, &setup->complex_data[i].im, 1, true);
      */
    }
#endif
    //print_unlock();
  }

  if (dir == GMX_FFT_REAL_TO_COMPLEX) {
    cudaMemcpy(setup->complex_data, setup->cdata + x * y * (z/2+1), x * y * (z/2+1) * sizeof(t_complex), cudaMemcpyDeviceToHost);
  } else {
    cudaMemcpy(setup->real_data, setup->rdata, x * y * z * sizeof(real), cudaMemcpyDeviceToHost);
  }

  if (check_vs_cpu(fft_gpu_flags)) {
    //print_lock();
    for (int ix = 0; ix < x; ++ix) {
      fprintf(stderr, "plane %d\n", ix);
      for (int iy = 0; iy < y; ++iy) {
	for (int iz = 0; iz < z / 2 + 1; ++iz) {
	  int i = (ix * y + iy) * (z/2+1) + iz;
	  fprintf(stderr, " %5.2f+i%5.2f",
		  setup->complex_data[i].re,
		  setup->complex_data[i].im);
	}
	fprintf(stderr, "\n");
      }
    }
    //print_unlock();
  }

  // FIX destroy plans after
  //cufftDestroy(setup->plan);

  return 0;
}

int
gmx_parallel_3dfft_destroy_gpu(gmx_parallel_3dfft_gpu_t    pfft_setup)
{
  //fprintf(stderr, "3dfft_destroy_gpu\n");
  gmx_parallel_3dfft_gpu_t setup = pfft_setup;

  cufftDestroy(setup->plan); // FIX double plan

  cudaFree((void **)setup->rdata);
  cudaFree((void **)setup->cdata);

  delete setup;
  return 0;
}
