#include "gromacs/utility/real.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"

#include "th_v.h"
#include "check.h"

using namespace thrust;

enum TH_V_ID {
  ID_G2T,
  ID_FSH,
  ID_NN,
  ID_XPTR,
  ID_IDXPTR,
  ID_FPTR,
  ID_END
};

static thread_vectors TH_V(32, ID_END);

extern gpu_flags interpol_gpu_flags;

template <typename T>
static T *raw_off(device_vector<T> &v, int off) {
  return thrust::raw_pointer_cast(&v[off]);
}

__global__ void calc_interpolation_idx_gpu_kernel
(int nx, int ny, int nz,
 int start_ix, int start_iy, int start_iz,
 real rxx, real ryx, real ryy, real rzx, real rzy, real rzz,
 int *g2tx, int *g2ty, int *g2tz,
 real *fshx, real *fshy,
 int *nnx, int *nny, int *nnz,
 real *xptr, real *yptr, real *zptr,
 int *idxxptr, int *idxyptr, int *idxzptr,
 real *fxptr, real *fyptr, real *fzptr,
 int n);

void calc_interpolation_idx_gpu_mid
(int nx, int ny, int nz,
 int start_ix, int start_iy, int start_iz,
 real rxx, real ryx, real ryy, real rzx, real rzy, real rzz,
 int *g2tx, int *g2ty, int *g2tz,
 real *fshx, real *fshy,
 int *nnx, int *nny, int *nnz,
 rvec *xptr_v, ivec *idxptr_v, rvec *fptr_v,
 int start, int end, int thread)
{
  int n = end - start;
  int n32 = (n + 31) / 32 * 32;

  local_vectors lv = TH_V.local(thread);

  device_vector<int> &g2t_d = lv.device<int>(ID_G2T, 3 * n32);
  thrust::copy(g2tx, g2tx + n, g2t_d.begin());
  thrust::copy(g2ty, g2ty + n, g2t_d.begin() + n32);
  thrust::copy(g2tz, g2tz + n, g2t_d.begin() + 2 * n32);

  device_vector<real> &fsh_d = lv.device<real>(ID_FSH, 5 * (nx + ny));
  thrust::copy(fshx, fshx + 5 * nx, fsh_d.begin());
  thrust::copy(fshy, fshy + 5 * ny, fsh_d.begin() + 5 * nx);

  device_vector<int> &nn_d = lv.device<int>(ID_NN, 5 * (nx + ny + nz));
  thrust::copy(nnx, nnx + 5 * nx, nn_d.begin());
  thrust::copy(nny, nny + 5 * ny, nn_d.begin() + 5 * nx);
  thrust::copy(nnz, nnz + 5 * nz, nn_d.begin() + 5 * (nx + ny));

  host_vector<real> &xptr_h = lv.host<real>(ID_XPTR, 3 * n32);
  host_vector<int> &idxptr_h = lv.host<int>(ID_IDXPTR, 3 * n32);
  host_vector<real> &fptr_h = lv.host<real>(ID_FPTR, 3 * n32);
  {
    int ix = 0, iy = n32, iz = 2 * n32;
    for (int i = start; i < end; i++)
    {
      real *xptr = xptr_v[i];
      xptr_h[ix++] = xptr[XX];
      xptr_h[iy++] = xptr[YY];
      xptr_h[iz++] = xptr[ZZ];
    }
  }

  device_vector<real> &xptr_d = lv.device<real>(ID_XPTR, 3 * n32);
  device_vector<int> &idxptr_d = lv.device<int>(ID_IDXPTR, 3 * n32);
  device_vector<real> &fptr_d = lv.device<real>(ID_FPTR, 3 * n32);

  xptr_d = xptr_h;

  int block_size = 32;
  int n_blocks = (n + block_size - 1) / block_size;

  calc_interpolation_idx_gpu_kernel<<<n_blocks, block_size>>>
    (nx, ny, nz, start_ix, start_iy, start_iz,
     rxx, ryx, ryy, rzx, rzy, rzz,

     thrust::raw_pointer_cast(&g2t_d[0]),
     thrust::raw_pointer_cast(&g2t_d[n32]),
     thrust::raw_pointer_cast(&g2t_d[2 * n32]),

     thrust::raw_pointer_cast(&fsh_d[0]),
     thrust::raw_pointer_cast(&fsh_d[5 * nx]),

     thrust::raw_pointer_cast(&nn_d[0]),
     thrust::raw_pointer_cast(&nn_d[5 * nx]),
     thrust::raw_pointer_cast(&nn_d[5 * (nx + ny)]),

     thrust::raw_pointer_cast(&xptr_d[0]),
     thrust::raw_pointer_cast(&xptr_d[n32]),
     thrust::raw_pointer_cast(&xptr_d[2 * n32]),

     thrust::raw_pointer_cast(&idxptr_d[0]),
     thrust::raw_pointer_cast(&idxptr_d[n32]),
     thrust::raw_pointer_cast(&idxptr_d[2 * n32]),

     thrust::raw_pointer_cast(&fptr_d[0]),
     thrust::raw_pointer_cast(&fptr_d[n32]),
     thrust::raw_pointer_cast(&fptr_d[2 * n32]),

     n);

  idxptr_h = idxptr_d;
  fptr_h = fptr_d;
  {
    int ix = 0, iy = n32, iz = 2 * n32;
    for (int i = start; i < end; i++)
    {
      int *idxptr = idxptr_v[i];
      real *fptr   = fptr_v[i];
      if (check_vs_cpu(interpol_gpu_flags)) {
	check_int("calc:idxx", &idxptr_h[ix], &idxptr[XX], 1, false);
	check_int("calc:idxy", &idxptr_h[iy], &idxptr[YY], 1, false);
	check_int("calc:idxz", &idxptr_h[iz], &idxptr[ZZ], 1, false);
	check_real("calc:fx", &fptr_h[ix], &fptr[XX], 1, false);
	check_real("calc:fy", &fptr_h[iy], &fptr[YY], 1, false);
	check_real("calc:fz", &fptr_h[iz], &fptr[ZZ], 1, false);
      }
      idxptr[XX] = idxptr_h[ix];
      idxptr[YY] = idxptr_h[iy];
      idxptr[ZZ] = idxptr_h[iz];
      fptr[XX] = fptr_h[ix];
      fptr[YY] = fptr_h[iy];
      fptr[ZZ] = fptr_h[iz];
      ++ix; ++iy; ++iz;
    }
  }
}

__global__ void calc_interpolation_idx_gpu_kernel
(int nx, int ny, int nz,
 int start_ix, int start_iy, int start_iz,
 real rxx, real ryx, real ryy, real rzx, real rzy, real rzz,
 int *g2tx, int *g2ty, int *g2tz,
 real *fshx, real *fshy,
 int *nnx, int *nny, int *nnz,
 real *xptr, real *yptr, real *zptr,
 int *idxxptr, int *idxyptr, int *idxzptr,
 real *fxptr, real *fyptr, real *fzptr,
 int n)
{

  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < n) {

        /* Fractional coordinates along box vectors, add 2.0 to make 100% sure we are positive for triclinic boxes */
	real tx, ty, tz;
        tx = nx * ( xptr[i] * rxx + yptr[i] * ryx + zptr[i] * rzx + 2.0 );
        ty = ny * (                 yptr[i] * ryy + zptr[i] * rzy + 2.0 );
        tz = nz * (                                 zptr[i] * rzz + 2.0 );

	int tix, tiy, tiz;
        tix = (int)(tx);
        tiy = (int)(ty);
        tiz = (int)(tz);

        /* Because decomposition only occurs in x and y,
         * we never have a fraction correction in z.
         */
        fxptr[i] = tx - tix + fshx[tix];
        fyptr[i] = ty - tiy + fshy[tiy];
        fzptr[i] = tz - tiz;

        idxxptr[i] = nnx[tix];
        idxyptr[i] = nny[tiy];
        idxzptr[i] = nnz[tiz];
    }
}
