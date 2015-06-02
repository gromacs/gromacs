#include "pme.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/math/vectypes.h"

#include <cuda.h>

#include "th_v.h"
#include "check.h"

enum TH_V_ID {
  ID_ATC_F,
  ID_GRID,
  ID_COEFFICIENT,
  ID_I,
  ID_I0, ID_J0, ID_K0,
  ID_THX, ID_THY, ID_THZ,
  ID_DTHX, ID_DTHY, ID_DTHZ,
  ID_END
};

static thread_vectors TH_V(32, ID_END);

typedef real *splinevec[DIM];
extern gpu_flags gather_gpu_flags;


#define DO_FSPLINE(order)                      \
    for (int ithx = 0; (ithx < order); ithx++)              \
    {                                              \
        int index_x = (i0[i]+ithx)*pny*pnz;               \
        real tx      = thx[iorder+ithx];                       \
        real dx      = dthx[iorder+ithx];                      \
                                               \
        for (int ithy = 0; (ithy < order); ithy++)          \
        {                                          \
            int index_xy = index_x+(j0[i]+ithy)*pnz;      \
            real ty       = thy[iorder+ithy];                  \
            real dy       = dthy[iorder+ithy];                 \
            real fxy1     = 0, fz1 = 0;		   \
                                               \
            for (int ithz = 0; (ithz < order); ithz++)      \
            {                                      \
                real gval  = grid[index_xy+(k0[i]+ithz)];  \
                fxy1 += thz[iorder+ithz]*gval;            \
                fz1  += dthz[iorder+ithz]*gval;           \
            }                                      \
            fx += dx*ty*fxy1;                      \
            fy += tx*dy*fxy1;                      \
            fz += tx*ty*fz1;                       \
        }                                          \
    }


static __global__ void gather_f_bsplines_kernel
(real *grid, int order, int n,
 int nx, int ny, int nz, int pnx, int pny, int pnz,
 real rxx, real ryx, real ryy, real rzx, real rzy, real rzz,
 real *thx, real *thy, real *thz, real *dthx, real *dthy, real *dthz,
 real *atc_f, real *coefficient_v, int *i0, int *j0, int *k0)
{
  /* sum forces for local particles */
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    real coefficient = coefficient_v[i];
    real fx     = 0;
    real fy     = 0;
    real fz     = 0;
    int iorder = i*order;
    int idim = i * DIM;

    switch (order)
    {
    case 4:
      DO_FSPLINE(4);
      break;
    case 5:
      DO_FSPLINE(5);
      break;
    default:
      DO_FSPLINE(order);
      break;
    }
    /*if (i < 10)
      printf("GPU gather_f_bsplines after DO_FSPLINE %d %f,%f,%f\n", i, (double)fx, (double)fy, (double)fz);*/

    atc_f[idim + XX] += -coefficient*( fx*nx*rxx );
    atc_f[idim + YY] += -coefficient*( fx*nx*ryx + fy*ny*ryy );
    atc_f[idim + ZZ] += -coefficient*( fx*nx*rzx + fy*ny*rzy + fz*nz*rzz );

    /*printf("kernel coeff=%f f=%f,%f,%f\n",
	   (double) coefficient,
	   (double) fx, (double) fy, (double) fz);*/

    /* Since the energy and not forces are interpolated
     * the net force might not be exactly zero.
     * This can be solved by also interpolating F, but
     * that comes at a cost.
     * A better hack is to remove the net force every
     * step, but that must be done at a higher level
     * since this routine doesn't see all atoms if running
     * in parallel. Don't know how important it is?  EL 990726
     */
  }
}


void gather_f_bsplines_gpu_2_pre
(gmx_bool bClearF,
 int *spline_ind, int spline_n,
 real *atc_coefficient, rvec *atc_f,
 real scale, int thread
 )
{
    // copy atc_f before cpu calcucation
    local_vectors lv = TH_V.local(thread);

    thrust::host_vector<real> &atc_f_h = lv.host<real>(ID_ATC_F, DIM * spline_n);
    thrust::host_vector<int> &i_h = lv.host<int>(ID_I, spline_n);
    int oo = 0;
    for (int ii = 0; ii < spline_n; ii++)
    {
        int i           = spline_ind[ii];
        real coefficient_i = scale*atc_coefficient[i];
        if (bClearF)
        {
            atc_f[i][XX] = 0;
            atc_f[i][YY] = 0;
            atc_f[i][ZZ] = 0;
        }

	if (coefficient_i == 0) {
	  continue;
	}

	atc_f_h[oo * DIM + XX] = atc_f[i][XX];
	atc_f_h[oo * DIM + YY] = atc_f[i][YY];
	atc_f_h[oo * DIM + ZZ] = atc_f[i][ZZ];
	i_h[oo] = i;
	oo++;
    }
}

void gather_f_bsplines_gpu_2
(real *grid, gmx_bool bClearF,
 int order,
 int nx, int ny, int nz, int pnx, int pny, int pnz,
 real rxx, real ryx, real ryy, real rzx, real rzy, real rzz,
 int *spline_ind, int spline_n,
 real *atc_coefficient, rvec *atc_f, ivec *atc_idx,
 splinevec *spline_theta, splinevec *spline_dtheta,
 real scale,
 int thread
 )
{
    int ndatatot = pnx*pny*pnz;

    if (check_vs_cpu(gather_gpu_flags)) {
      fprintf(stderr, "gather_f_bsplines_gpu_2 %dx%dx%d=%d %d\n", pnx, pny, pnz, ndatatot, spline_n);
    }

    local_vectors lv = TH_V.local(thread);

    thrust::device_vector<real> &grid_d = lv.device<real>(ID_GRID, ndatatot);
    thrust::copy(grid, grid + ndatatot, grid_d.begin());

    thrust::host_vector<real> &atc_f_h = lv.host<real>(ID_ATC_F, DIM * spline_n);

    thrust::host_vector<real> &coefficient_h = lv.host<real>(ID_COEFFICIENT, spline_n);

    thrust::host_vector<int> &i_h = lv.host<int>(ID_I, spline_n);
    thrust::host_vector<int> &i0_h = lv.host<int>(ID_I0, spline_n);
    thrust::host_vector<int> &j0_h = lv.host<int>(ID_J0, spline_n);
    thrust::host_vector<int> &k0_h = lv.host<int>(ID_K0, spline_n);

    thrust::host_vector<real> &thx_h = lv.host<real>(ID_THX, order * spline_n);
    thrust::host_vector<real> &thy_h = lv.host<real>(ID_THY, order * spline_n);
    thrust::host_vector<real> &thz_h = lv.host<real>(ID_THZ, order * spline_n);
    thrust::host_vector<real> &dthx_h = lv.host<real>(ID_DTHX, order * spline_n);
    thrust::host_vector<real> &dthy_h = lv.host<real>(ID_DTHY, order * spline_n);
    thrust::host_vector<real> &dthz_h = lv.host<real>(ID_DTHZ, order * spline_n);

    int oo = 0;
    for (int ii = 0; ii < spline_n; ii++)
    {
        int i           = spline_ind[ii];
        real coefficient_i = scale*atc_coefficient[i];
        if (bClearF)
        {
            atc_f[i][XX] = 0;
            atc_f[i][YY] = 0;
            atc_f[i][ZZ] = 0;
        }

	if (coefficient_i == 0) {
	  continue;
	}

	coefficient_h[oo] = coefficient_i;
	int *idxptr = atc_idx[i];
	//atc_f_h force-copying is in gather_f_bsplines_gpu_2_pre()
	i_h[oo] = i;
	i0_h[oo] = idxptr[XX];
	j0_h[oo] = idxptr[YY];
	k0_h[oo] = idxptr[ZZ];
	int iiorder = ii*order;
	int ooorder = oo*order;
	for (int o = 0; o < order; ++o) {
	  thx_h[ooorder + o] = (*spline_theta)[XX][iiorder + o];
	  thy_h[ooorder + o] = (*spline_theta)[YY][iiorder + o];
	  thz_h[ooorder + o] = (*spline_theta)[ZZ][iiorder + o];
	  dthx_h[ooorder + o] = (*spline_dtheta)[XX][iiorder + o];
	  dthy_h[ooorder + o] = (*spline_dtheta)[YY][iiorder + o];
	  dthz_h[ooorder + o] = (*spline_dtheta)[ZZ][iiorder + o];
	}
	++oo;
    }

    int n = oo;

    thrust::device_vector<real> &atc_f_d = lv.device<real>(ID_ATC_F, DIM * n);
    thrust::device_vector<real> &coefficient_d = lv.device<real>(ID_COEFFICIENT, n);
    thrust::device_vector<int> &i0_d = lv.device<int>(ID_I0, n);
    thrust::device_vector<int> &j0_d = lv.device<int>(ID_J0, n);
    thrust::device_vector<int> &k0_d = lv.device<int>(ID_K0, n);
    thrust::device_vector<real> &thx_d = lv.device<real>(ID_THX, order * n);
    thrust::device_vector<real> &thy_d = lv.device<real>(ID_THY, order * n);
    thrust::device_vector<real> &thz_d = lv.device<real>(ID_THZ, order * n);
    thrust::device_vector<real> &dthx_d = lv.device<real>(ID_DTHX, order * n);
    thrust::device_vector<real> &dthy_d = lv.device<real>(ID_DTHY, order * n);
    thrust::device_vector<real> &dthz_d = lv.device<real>(ID_DTHZ, order * n);
    atc_f_d = atc_f_h;
    coefficient_d = coefficient_h;
    i0_d = i0_h;
    j0_d = j0_h;
    k0_d = k0_h;
    thx_d = thx_h;
    thy_d = thy_h;
    thz_d = thz_h;
    dthx_d = dthx_h;
    dthy_d = dthy_h;
    dthz_d = dthz_h;


    int block_size = 64;
    int n_blocks = (n + block_size - 1) / block_size;

    gather_f_bsplines_kernel<<<n_blocks, block_size>>>
      (thrust::raw_pointer_cast(&grid_d[0]),
       order, n,
       nx, ny, nz, pnx, pny, pnz,
       rxx, ryx, ryy, rzx, rzy, rzz,
       thrust::raw_pointer_cast(&thx_d[0]),
       thrust::raw_pointer_cast(&thy_d[0]),
       thrust::raw_pointer_cast(&thz_d[0]),
       thrust::raw_pointer_cast(&dthx_d[0]),
       thrust::raw_pointer_cast(&dthy_d[0]),
       thrust::raw_pointer_cast(&dthz_d[0]),
       thrust::raw_pointer_cast(&atc_f_d[0]),
       thrust::raw_pointer_cast(&coefficient_d[0]),
       thrust::raw_pointer_cast(&i0_d[0]),
       thrust::raw_pointer_cast(&j0_d[0]),
       thrust::raw_pointer_cast(&k0_d[0]));

    atc_f_h = atc_f_d;

    for (int ii = 0; ii < n; ii++)
    {
        int i = i_h[ii];
	if (check_vs_cpu(gather_gpu_flags)) {
	  fprintf(stderr, "check %d=%.2f,%.2f,%.2f ", i,
		  atc_f[i][XX],
		  atc_f[i][YY],
		  atc_f[i][ZZ]
		  );
	  check_real("atc_f", &atc_f_h[ii * DIM], atc_f[i], DIM, false);
	}

	atc_f[i][XX] = atc_f_h[ii * DIM + XX];
	atc_f[i][YY] = atc_f_h[ii * DIM + YY];
	atc_f[i][ZZ] = atc_f_h[ii * DIM + ZZ];
    }
}
