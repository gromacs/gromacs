#include "pme.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"

#include <cuda.h>

#include "th_v.h"
#include "check.h"

enum TH_V_ID {
  ID_GRID,
  ID_PME_BSP_MOD_X,
  ID_PME_BSP_MOD_Y,
  ID_PME_BSP_MOD_Z,
  ID_ENERGY,
  ID_VIRIAL,
  ID_END
};

static thread_vectors TH_V(32, ID_END);
extern gpu_flags solve_gpu_flags;


typedef real *splinevec[DIM];

/* Pascal triangle coefficients used in solve_pme_lj_yzx, only need to do 4 calculations due to symmetry */
static const __constant__ real lb_scale_factor_symm[] = { 2.0/64, 12.0/64, 30.0/64, 20.0/64 };


/*__device__ gmx_inline static void calc_exponentials_q_one(const real f, real &d, real &r, real &e)
{
  d = 1.0/d;
  r = expf(r);
  e = f*r*d;
  }*/

static const real sqrt_M_PI = sqrt(M_PI);
static __constant__ real sqrt_M_PI_d;

/*__device__ gmx_inline static void calc_exponentials_lj_one(real &r, real &tmp2, real &d)
{
  d = 1.0/d;
  r = exp(r);
  real mk = tmp2;
  tmp2 = sqrt_M_PI_d*mk*erfcf(mk);
  }*/


__global__ void solve_pme_yzx_iyz_loop_kernel
(int iyz0, int iyz1, int local_ndata_ZZ, int local_ndata_XX,
 int local_offset_XX, int local_offset_YY, int local_offset_ZZ,
 int local_size_XX, int local_size_YY, int local_size_ZZ,
 int nx, int ny, int nz,
 real rxx, real ryx, real ryy, real rzx, real rzy, real rzz,
 real elfac,
 //splinevec pme_bsp_mod,
 real *pme_bsp_mod_XX, real *pme_bsp_mod_YY, real *pme_bsp_mod_ZZ,
 t_complex *grid,
 real ewaldcoeff, real vol,
 gmx_bool bEnerVir,
 real *energy_v, real *virial_v);


int solve_pme_yzx_gpu(real pme_epsilon_r,
		      int nx, int ny, int nz,
		      ivec complex_order, ivec local_ndata, ivec local_offset, ivec local_size,
		      real rxx, real ryx, real ryy, real rzx, real rzy, real rzz,
		      //real *mhx, real *mhy, real *mhz, real *m2, real *denom, real *tmp1, real *eterm, real *m2inv,
		      splinevec pme_bsp_mod,
		      matrix work_vir_q, real *work_energy_q,
		      t_complex *grid,
		      real ewaldcoeff, real vol,
		      gmx_bool bEnerVir,
		      int nthread, int thread)
{
    /* do recip sum over local cells in grid */
    /* y major, z middle, x minor or continuous */
    //t_complex *p0;
    //int     kx, ky, kz, maxkx, maxky, maxkz;
    int     iyz0, iyz1; //, iyz, iy, iz, kxstart, kxend;
    // real    mx, my, mz;
    //real    factor = M_PI*M_PI/(ewaldcoeff*ewaldcoeff);
    //real    ets2, struct2, vfactor, ets2vf;
    //real    d1, d2;
    real energy = 0;
    //real    by, bz;
    real    virxx = 0, virxy = 0, virxz = 0, viryy = 0, viryz = 0, virzz = 0;
    //real    mhxk, mhyk, mhzk, m2k;
    //real    corner_fac;
    real    elfac;

    elfac = ONE_4PI_EPS0/pme_epsilon_r;

    /* Dimensions should be identical for A/B grid, so we just use A here */
    /*
      TODO: Dimensions are passed in for now. call complex limits elsewhere?
    gmx_parallel_3dfft_complex_limits(pme->pfft_setup[PME_GRID_QA],
                                      complex_order,
                                      local_ndata,
                                      local_offset,
                                      local_size);
    gmx_parallel_3dfft_complex_limits_gpu(pme->pfft_setup_gpu[PME_GRID_QA],
                                      complex_order,
                                      local_ndata,
                                      local_offset,
                                      local_size);
    */


    iyz0 = local_ndata[YY]*local_ndata[ZZ]* thread   /nthread;
    iyz1 = local_ndata[YY]*local_ndata[ZZ]*(thread+1)/nthread;

    const int block_size = 32;
    int n = iyz1 - iyz0;
    int n_blocks = (n + block_size - 1) / block_size;

    cudaMemcpyToSymbol( &sqrt_M_PI_d, &sqrt_M_PI, sizeof(real));
    printf("local_size[XX] %d local_ndata[XX] %d\n",
	   local_size[XX], local_ndata[XX]);
    int grid_size = local_size[YY] * local_size[ZZ] * local_ndata[XX]; // * local_size[XX];
    local_vectors lv = TH_V.local(thread);
    thrust::device_vector<t_complex> &grid_d = lv.device<t_complex>(ID_GRID, grid_size);
    thrust::device_vector<real> &pme_bsp_mod_x_d = lv.device<real>(ID_PME_BSP_MOD_X, nx);
    thrust::device_vector<real> &pme_bsp_mod_y_d = lv.device<real>(ID_PME_BSP_MOD_Y, ny);
    thrust::device_vector<real> &pme_bsp_mod_z_d = lv.device<real>(ID_PME_BSP_MOD_Z, nz);
    thrust::device_vector<real> &energy_d = lv.device<real>(ID_ENERGY, n);
    thrust::device_vector<real> &virial_d = lv.device<real>(ID_VIRIAL, 6 * n);
    fprintf(stderr, "grid copy %p to %p, end val %f,%f thread %d/%d\n",
	    grid, thrust::raw_pointer_cast(&grid_d[0]),
	    (double) grid[grid_size - local_size[XX] + local_ndata[XX] - 1].re,
	    (double) grid[grid_size - local_size[XX] + local_ndata[XX] - 1].im,
	    thread, nthread);
    thrust::copy(grid, grid + grid_size, grid_d.begin());
    thrust::copy(pme_bsp_mod[XX], pme_bsp_mod[XX] + nx, pme_bsp_mod_x_d.begin());
    thrust::copy(pme_bsp_mod[YY], pme_bsp_mod[YY] + ny, pme_bsp_mod_y_d.begin());
    thrust::copy(pme_bsp_mod[ZZ], pme_bsp_mod[ZZ] + nz, pme_bsp_mod_z_d.begin());
    fprintf(stderr, "grid copy after\n");

    solve_pme_yzx_iyz_loop_kernel<<<n_blocks, block_size>>>
      (iyz0, iyz1, local_ndata[ZZ], local_ndata[XX],
       local_offset[XX], local_offset[YY], local_offset[ZZ],
       local_size[XX], local_size[YY], local_size[ZZ],
       nx, ny, nz, rxx, ryx, ryy, rzx, rzy, rzz,
       elfac,
       thrust::raw_pointer_cast(&pme_bsp_mod_x_d[0]),
       thrust::raw_pointer_cast(&pme_bsp_mod_y_d[0]),
       thrust::raw_pointer_cast(&pme_bsp_mod_z_d[0]),
       thrust::raw_pointer_cast(&grid_d[0]), ewaldcoeff, vol, bEnerVir,
       thrust::raw_pointer_cast(&energy_d[0]),
       thrust::raw_pointer_cast(&virial_d[0]));

    if (check_vs_cpu(solve_gpu_flags)) {
      thrust::host_vector<t_complex> &grid_h = lv.host<t_complex>(ID_GRID, grid_size);
      //grid_h = grid_d;
      thrust::copy(grid_d.begin(), grid_d.end(), grid_h.begin());
      check_real("grid_q",
		 (real *) thrust::raw_pointer_cast(&grid_h[0]),
		 (real *) grid, 2 * grid_size, false);
    }
    thrust::copy(grid_d.begin(), grid_d.end(), grid);

    if (bEnerVir)
    {
      thrust::host_vector<real> &energy_h = lv.host<real>(ID_ENERGY, n);
      thrust::host_vector<real> &virial_h = lv.host<real>(ID_VIRIAL, 6*n);
      thrust::copy(energy_d.begin(), energy_d.begin() + n, energy_h.begin());
      thrust::copy(virial_d.begin(), virial_d.begin() + 6*n, virial_h.begin());
      for (int i = 0, j = 0; i < n; ++i) {
	energy += energy_h[i];
	virxx += virial_h[j++];
	viryy += virial_h[j++];
	virzz += virial_h[j++];
	virxy += virial_h[j++];
	virxz += virial_h[j++];
	viryz += virial_h[j++];
      }
      if (check_vs_cpu(solve_gpu_flags)) {
	real t;
	t = 0.25*virxx; check_real("viqxx", &t, &work_vir_q[XX][XX], 1, false);
	t = 0.25*viryy; check_real("viqyy", &t, &work_vir_q[YY][YY], 1, false);
	t = 0.25*virzz; check_real("viqzz", &t, &work_vir_q[ZZ][ZZ], 1, false);
	t = 0.25*virxy; check_real("viqxy", &t, &work_vir_q[XX][YY], 1, false);
	t = 0.25*virxz; check_real("viqxz", &t, &work_vir_q[XX][ZZ], 1, false);
	t = 0.25*viryz; check_real("viqyz", &t, &work_vir_q[YY][ZZ], 1, false);
	t = 0.25*virxy; check_real("virxy", &t, &work_vir_q[YY][XX], 1, false);
	t = 0.25*virxz; check_real("virxz", &t, &work_vir_q[ZZ][XX], 1, false);
	t = 0.25*viryz; check_real("viryz", &t, &work_vir_q[ZZ][YY], 1, false);
	t = 0.5*energy; check_real("enerq", &t, &*work_energy_q, 1, false);
      }

        /* Update virial with local values.
         * The virial is symmetric by definition.
         * this virial seems ok for isotropic scaling, but I'm
         * experiencing problems on semiisotropic membranes.
         * IS THAT COMMENT STILL VALID??? (DvdS, 2001/02/07).
         */
        work_vir_q[XX][XX] = 0.25*virxx;
        work_vir_q[YY][YY] = 0.25*viryy;
        work_vir_q[ZZ][ZZ] = 0.25*virzz;
        work_vir_q[XX][YY] = work_vir_q[YY][XX] = 0.25*virxy;
        work_vir_q[XX][ZZ] = work_vir_q[ZZ][XX] = 0.25*virxz;
        work_vir_q[YY][ZZ] = work_vir_q[ZZ][YY] = 0.25*viryz;

        /* This energy should be corrected for a charged system */
        *work_energy_q = 0.5*energy;
    }

    /* Return the loop count */
    return local_ndata[YY]*local_ndata[XX];
}


__global__ void solve_pme_yzx_iyz_loop_kernel
(int iyz0, int iyz1, int local_ndata_ZZ, int local_ndata_XX,
 int local_offset_XX, int local_offset_YY, int local_offset_ZZ,
 int local_size_XX, int local_size_YY, int local_size_ZZ,
 int nx, int ny, int nz,
 real rxx, real ryx, real ryy, real rzx, real rzy, real rzz,
 real elfac,
 //splinevec pme_bsp_mod,
 real *pme_bsp_mod_XX, real *pme_bsp_mod_YY, real *pme_bsp_mod_ZZ,
 t_complex *grid,
 real ewaldcoeff, real vol,
 gmx_bool bEnerVir,
 real *energy_v, real *virial_v)
{
  const real factor = M_PI*M_PI/(ewaldcoeff*ewaldcoeff);

  int maxkx = (nx+1)/2;
  int maxky = (ny+1)/2;
  //int maxkz = nz/2+1;
  //(void) maxkz; // unused


  real energy = 0;
  real    virxx = 0, virxy = 0, virxz = 0, viryy = 0, viryz = 0, virzz = 0;


  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i == 0) { printf("solve_kernel start\n"); }
  int iyz = iyz0 + i;
  if (iyz < iyz1)
  {
        int iy = iyz/local_ndata_ZZ;
        int iz = iyz - iy*local_ndata_ZZ;

        int ky = iy + local_offset_YY;
	real my;

        if (ky < maxky)
        {
            my = ky;
        }
        else
        {
            my = (ky - ny);
        }

        real by = M_PI*vol*pme_bsp_mod_YY[ky];

        int kz = iz + local_offset_ZZ;

        real mz = kz;

        real bz = pme_bsp_mod_ZZ[kz];

        /* 0.5 correction for corner points */
        real corner_fac = 1;
        if (kz == 0 || kz == (nz+1)/2)
        {
            corner_fac = 0.5;
        }

        t_complex *p0 = grid + iy*local_size_ZZ*local_size_XX + iz*local_size_XX;

	int kxstart;
	// FIXME: avoid warp divergence for the (0,0,0) point (how bad is it?)
        /* We should skip the k-space point (0,0,0) */
        /* Note that since here x is the minor index, local_offset[XX]=0 */
        if (local_offset_XX > 0 || ky > 0 || kz > 0)
        {
            kxstart = local_offset_XX;
        }
        else
        {
            kxstart = local_offset_XX + 1;
            p0++;
        }
        int kxend = local_offset_XX + local_ndata_XX;

	real mx, mhxk, mhyk, mhzk, m2k;
	real ets2, struct2, vfactor, ets2vf;

        if (bEnerVir)
        {
            /* More expensive inner loop, especially because of the storage
             * of the mh elements in array's.
             * Because x is the minor grid index, all mh elements
             * depend on kx for triclinic unit cells.
             */

            // /* Two explicit loops to avoid a conditional inside the loop */
	    // NOTE: on gpu, keep the conditional. shouldn't be too bad?
	    for (int kx = kxstart; kx < kxend; kx++, p0++)
            {
	      if (i == 0) { printf("solve_kernel %d\n", kx); }
	        mx = kx < maxkx ? kx : (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                //mhx[kx]   = mhxk;
                //mhy[kx]   = mhyk;
                //mhz[kx]   = mhzk;
                //m2[kx]    = m2k;
                real denom = m2k*bz*by*pme_bsp_mod_XX[kx];
                real tmp1  = -factor*m2k;

		if (iyz == iyz0 && kx == kxstart)
		  printf
		    ("SOLVE_gpu mhxk %f mhyk %f mhzk %f m2k %f denom %f tmp1 %f\n", (double) mhxk, (double) mhyk, (double) mhzk, (double) m2k, (double) denom, (double) tmp1);

                real m2invk = 1.0/m2k;

		//calc_exponentials_q_one(elfac, denom, tmp1, eterm);
		denom = 1.0/denom;
		tmp1 = exp(tmp1);
		real etermk = elfac*tmp1*denom;

                real d1      = p0->re;
                real d2      = p0->im;

                p0->re  = d1*etermk;
                p0->im  = d2*etermk;

                struct2 = 2.0*(d1*d1+d2*d2);

                real tmp1k = etermk*struct2;

                ets2     = corner_fac*tmp1k;
                vfactor  = (factor*m2k + 1.0)*2.0*m2invk;
                energy  += ets2;

                ets2vf   = ets2*vfactor;
                virxx   += ets2vf*mhxk*mhxk - ets2;
                virxy   += ets2vf*mhxk*mhyk;
                virxz   += ets2vf*mhxk*mhzk;
                viryy   += ets2vf*mhyk*mhyk - ets2;
                viryz   += ets2vf*mhyk*mhzk;
                virzz   += ets2vf*mhzk*mhzk - ets2;
            }
        }
        else
        {
            /* We don't need to calculate the energy and the virial.
             * In this case the triclinic overhead is small.
             */

            /* Two explicit loops to avoid a conditional inside the loop */
	    // NOTE: on gpu, keep the conditional. shouldn't be too bad?
	    for (int kx = kxstart; kx < kxend; kx++, p0++)
            {
	        mx = kx < maxkx ? kx : (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                real denom = m2k*bz*by*pme_bsp_mod_XX[kx];
                real tmp1  = -factor*m2k;

		//calc_exponentials_q_one(elfac, denom, tmp1, eterm);
		denom = 1.0/denom;
		tmp1 = exp(tmp1);
		real etermk = elfac*tmp1*denom;

                real d1      = p0->re;
                real d2      = p0->im;

                p0->re  = d1*etermk;
                p0->im  = d2*etermk;
            }
        }
	energy_v[i] = energy;
	virial_v[6*i+0] = virxx;
	virial_v[6*i+1] = viryy;
	virial_v[6*i+2] = virzz;
	virial_v[6*i+3] = virxy;
	virial_v[6*i+4] = virxz;
	virial_v[6*i+5] = viryz;
    }
}


__global__ void solve_pme_lj_yzx_iyz_loop_kernel
(int iyz0, int iyz1, int local_ndata_ZZ, int local_ndata_XX,
 int local_offset_XX, int local_offset_YY, int local_offset_ZZ,
 int local_size_XX, int local_size_YY, int local_size_ZZ,
 int nx, int ny, int nz,
 real rxx, real ryx, real ryy, real rzx, real rzy, real rzz,
 //real elfac,
 //splinevec pme_bsp_mod,
 real *pme_bsp_mod_XX, real *pme_bsp_mod_YY, real *pme_bsp_mod_ZZ,
 t_complex *grid_v, gmx_bool bLB,
 real ewaldcoeff, real vol,
 gmx_bool bEnerVir,
 real *energy_v, real *virial_v);


int solve_pme_lj_yzx_gpu(int nx, int ny, int nz,
			 ivec complex_order, ivec local_ndata, ivec local_offset, ivec local_size,
			 real rxx, real ryx, real ryy, real rzx, real rzy, real rzz,
			 //real *mhx, real *mhy, real *mhz, real *m2, real *denom, real *tmp1, real *tmp2,
			 splinevec pme_bsp_mod,
			 matrix work_vir_lj, real *work_energy_lj,
			 t_complex **grid, gmx_bool bLB,
			 real ewaldcoeff, real vol,
			 gmx_bool bEnerVir, int nthread, int thread)
{
    /* do recip sum over local cells in grid */
    /* y major, z middle, x minor or continuous */
    //int     ig, gcount;
    //int     kx, ky, kz, maxkx, maxky, maxkz;
    int     /*iy,*/ iyz0, iyz1; //, iyz, iz, kxstart, kxend;
    //real    mx, my, mz;
    //real    factor = M_PI*M_PI/(ewaldcoeff*ewaldcoeff);
    //real    ets2, ets2vf;
    //real    eterm, vterm, d1, d2;
    real energy = 0;
    //real    by, bz;
    real    virxx = 0, virxy = 0, virxz = 0, viryy = 0, viryz = 0, virzz = 0;
    //real    mhxk, mhyk, mhzk, m2k;
    //real    mk;
    //real    corner_fac;

    /* Dimensions should be identical for A/B grid, so we just use A here */
    /* Dimensions are passed in. TODO: call elsewhere?
    gmx_parallel_3dfft_complex_limits(pme->pfft_setup[PME_GRID_C6A],
                                      complex_order,
                                      local_ndata,
                                      local_offset,
                                      local_size);
    gmx_parallel_3dfft_complex_limits_gpu(pme->pfft_setup_gpu[PME_GRID_C6A],
                                      complex_order,
                                      local_ndata,
                                      local_offset,
                                      local_size);
    */

    iyz0 = local_ndata[YY]*local_ndata[ZZ]* thread   /nthread;
    iyz1 = local_ndata[YY]*local_ndata[ZZ]*(thread+1)/nthread;

    cudaMemcpyToSymbol( &sqrt_M_PI_d, &sqrt_M_PI, sizeof(real));

    const int block_size = 32;
    int n = iyz1 - iyz0;
    int n_blocks = (n + block_size - 1) / block_size;

    int grid_size = local_size[YY] * local_size[ZZ] * local_size[XX];
    local_vectors lv = TH_V.local(thread);
    thrust::device_vector<t_complex> &grid_d = lv.device<t_complex>(ID_GRID, 6 * grid_size);
    thrust::device_vector<real> &pme_bsp_mod_x_d = lv.device<real>(ID_PME_BSP_MOD_X, nx);
    thrust::device_vector<real> &pme_bsp_mod_y_d = lv.device<real>(ID_PME_BSP_MOD_Y, ny);
    thrust::device_vector<real> &pme_bsp_mod_z_d = lv.device<real>(ID_PME_BSP_MOD_Z, nz);
    thrust::device_vector<real> &energy_d = lv.device<real>(ID_ENERGY, n);
    thrust::device_vector<real> &virial_d = lv.device<real>(ID_VIRIAL, 6 * n);
    for (int ig = 0; ig < 6; ++ig) {
      thrust::copy(grid[ig], grid[ig] + grid_size,
		   grid_d.begin() + ig * grid_size);
    }
    thrust::copy(pme_bsp_mod[XX], pme_bsp_mod[XX] + nx, pme_bsp_mod_x_d.begin());
    thrust::copy(pme_bsp_mod[YY], pme_bsp_mod[YY] + ny, pme_bsp_mod_y_d.begin());
    thrust::copy(pme_bsp_mod[ZZ], pme_bsp_mod[ZZ] + nz, pme_bsp_mod_z_d.begin());

    solve_pme_lj_yzx_iyz_loop_kernel<<<n_blocks, block_size>>>
      (iyz0, iyz1, local_ndata[ZZ], local_ndata[XX],
       local_offset[XX], local_offset[YY], local_offset[ZZ],
       local_size[XX], local_size[YY], local_size[ZZ],
       nx, ny, nz, rxx, ryx, ryy, rzx, rzy, rzz,
       //elfac,
       //pme_bsp_mod,
       thrust::raw_pointer_cast(&pme_bsp_mod_x_d[0]),
       thrust::raw_pointer_cast(&pme_bsp_mod_y_d[0]),
       thrust::raw_pointer_cast(&pme_bsp_mod_z_d[0]),
       thrust::raw_pointer_cast(&grid_d[0]), bLB, ewaldcoeff, vol, bEnerVir,
       thrust::raw_pointer_cast(&energy_d[0]),
       thrust::raw_pointer_cast(&virial_d[0]));

    for (int ig = 0; ig < 6; ++ig) {
      if (check_vs_cpu(solve_gpu_flags)) {
	check_real("grid_lj",
		   (real *) thrust::raw_pointer_cast(&grid_d[ig * grid_size]),
		   (real *) grid[ig], 2 * grid_size, true);
      }
      thrust::copy(grid_d.begin() + ig * grid_size,
		   grid_d.begin() + (ig + 1) * grid_size, grid[ig]);
    }

    if (bEnerVir)
    {
      thrust::host_vector<real> &energy_h = lv.host<real>(ID_ENERGY, n);
      thrust::host_vector<real> &virial_h = lv.host<real>(ID_VIRIAL, 6 * n);
      thrust::copy(energy_d.begin(), energy_d.end(), energy_h.begin());
      thrust::copy(virial_d.begin(), virial_d.end(), virial_h.begin());
      for (int i = 0, j = 0; i < n; ++i) {
	energy += energy_h[i];
	virxx += virial_h[j++];
	viryy += virial_h[j++];
	virzz += virial_h[j++];
	virxy += virial_h[j++];
	virxz += virial_h[j++];
	viryz += virial_h[j++];
      }
      if (check_vs_cpu(solve_gpu_flags)) {
	real t;
	t = 0.25*virxx; check_real("vlxx", &t, &work_vir_lj[XX][XX], 1, false);
	t = 0.25*viryy; check_real("vlyy", &t, &work_vir_lj[YY][YY], 1, false);
	t = 0.25*virzz; check_real("vlzz", &t, &work_vir_lj[ZZ][ZZ], 1, false);
	t = 0.25*virxy; check_real("vlxy", &t, &work_vir_lj[XX][YY], 1, false);
	t = 0.25*virxz; check_real("vlxz", &t, &work_vir_lj[XX][ZZ], 1, false);
	t = 0.25*viryz; check_real("vlyz", &t, &work_vir_lj[YY][ZZ], 1, false);
	t = 0.25*virxy; check_real("vjxy", &t, &work_vir_lj[YY][XX], 1, false);
	t = 0.25*virxz; check_real("vjxz", &t, &work_vir_lj[ZZ][XX], 1, false);
	t = 0.25*viryz; check_real("vjyz", &t, &work_vir_lj[ZZ][YY], 1, false);
	t = 0.5*energy; check_real("enerl", &t, &*work_energy_lj, 1, false);
      }

        work_vir_lj[XX][XX] = 0.25*virxx;
        work_vir_lj[YY][YY] = 0.25*viryy;
        work_vir_lj[ZZ][ZZ] = 0.25*virzz;
        work_vir_lj[XX][YY] = work_vir_lj[YY][XX] = 0.25*virxy;
        work_vir_lj[XX][ZZ] = work_vir_lj[ZZ][XX] = 0.25*virxz;
        work_vir_lj[YY][ZZ] = work_vir_lj[ZZ][YY] = 0.25*viryz;

        /* This energy should be corrected for a charged system */
        *work_energy_lj = 0.5*energy;
    }
    /* Return the loop count */
    return local_ndata[YY]*local_ndata[XX];
}


__global__ void solve_pme_lj_yzx_iyz_loop_kernel
(int iyz0, int iyz1, int local_ndata_ZZ, int local_ndata_XX,
 int local_offset_XX, int local_offset_YY, int local_offset_ZZ,
 int local_size_XX, int local_size_YY, int local_size_ZZ,
 int nx, int ny, int nz,
 real rxx, real ryx, real ryy, real rzx, real rzy, real rzz,
 //real elfac,
 //splinevec pme_bsp_mod,
 real *pme_bsp_mod_XX, real *pme_bsp_mod_YY, real *pme_bsp_mod_ZZ,
 t_complex *grid_v, gmx_bool bLB,
 real ewaldcoeff, real vol,
 gmx_bool bEnerVir,
 real *energy_v, real *virial_v) {

  const int grid_size = local_size_YY * local_size_ZZ * local_size_XX;
  const real factor = M_PI*M_PI/(ewaldcoeff*ewaldcoeff);

  int maxkx = (nx+1)/2;
  int maxky = (ny+1)/2;
  //int maxkz = nz/2+1;
  //(void) maxkz; // unused


  real energy = 0;
  real    virxx = 0, virxy = 0, virxz = 0, viryy = 0, viryz = 0, virzz = 0;


  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int iyz = iyz0 + i;
  if (iyz < iyz1)
  {
        int iy = iyz/local_ndata_ZZ;
        int iz = iyz - iy*local_ndata_ZZ;

        int ky = iy + local_offset_YY;
	real my;

        if (ky < maxky)
        {
            my = ky;
        }
        else
        {
            my = (ky - ny);
        }

        real by = 3.0*vol*pme_bsp_mod_YY[ky]
            / (M_PI*sqrt(M_PI)*ewaldcoeff*ewaldcoeff*ewaldcoeff);

        int kz = iz + local_offset_ZZ;

        real mz = kz;

        real bz = pme_bsp_mod_ZZ[kz];

        /* 0.5 correction for corner points */
        real corner_fac = 1;
        if (kz == 0 || kz == (nz+1)/2)
        {
            corner_fac = 0.5;
        }

        int kxstart = local_offset_XX;
        int kxend   = local_offset_XX + local_ndata_XX;

	real mx, mhxk, mhyk, mhzk, m2k;

        if (bEnerVir)
        {
	  t_complex *p0 = grid_v/*[0]*/ + iy*local_size_ZZ*local_size_XX + iz*local_size_XX;
            /* More expensive inner loop, especially because of the
             * storage of the mh elements in array's.  Because x is the
             * minor grid index, all mh elements depend on kx for
             * triclinic unit cells.
             */

            // /* Two explicit loops to avoid a conditional inside the loop */
	    // NOTE: on gpu, keep the conditional. shouldn't be too bad?
	    for (int kx = kxstart; kx < kxend; kx++, p0++)
            {
	        mx = kx < maxkx ? kx : (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                real denomk = bz*by*pme_bsp_mod_XX[kx];
                real tmp1k  = -factor*m2k;
                real tmp2k  = sqrt(factor*m2k);

		//calc_exponentials_lj_one(tmp1k, tmp2k, denomk); // r tmp2 d
		denomk = 1.0/denomk;
		tmp1k = exp(tmp1k);
		real mk = tmp2k;
		tmp2k = sqrt_M_PI_d*mk*erfcf(mk);

                m2k   = factor*m2k;
                real eterm = -((1.0 - 2.0*m2k)*tmp1k
                          + 2.0*m2k*tmp2k);
                real vterm    = 3.0*(-tmp1k + tmp2k);
                tmp1k = eterm*denomk;
                tmp2k = vterm*denomk;

		if (!bLB)
		{
		  real d1      = p0->re;
		  real d2      = p0->im;

		  eterm   = tmp1k;
		  vterm   = tmp2k;
		  p0->re  = d1*eterm;
		  p0->im  = d2*eterm;

		  real struct2 = 2.0*(d1*d1+d2*d2);

		  tmp1k = eterm*struct2;
		  tmp2k = vterm*struct2;
		}
		else
		{
		  //real *struct2 = denom;
		  real  str2;

		  real struct2k = 0.0;

		  /* Due to symmetry we only need to calculate 4 of the 7 terms */
		  for (int ig = 0; ig <= 3; ++ig)
		    {
		      //t_complex *p0, *p1;
		      real       scale;

		      t_complex *p0k    = grid_v/*[ig]*/ + ig*grid_size + iy*local_size_ZZ*local_size_XX + iz*local_size_XX + (kx - kxstart);
		      t_complex *p1k    = grid_v/*[6-ig]*/ + (6-ig)*grid_size + iy*local_size_ZZ*local_size_XX + iz*local_size_XX + (kx - kxstart);
		      scale = 2.0*lb_scale_factor_symm[ig];
		      struct2k += scale*(p0k->re*p1k->re + p0k->im*p1k->im);
		    }
		  for (int ig = 0; ig <= 6; ++ig)
		  {
		    //t_complex *p0;

                    t_complex *p0k = grid_v/*[ig]*/ + ig*grid_size + iy*local_size_ZZ*local_size_XX + iz*local_size_XX + (kx - kxstart);

		    real d1     = p0k->re;
		    real d2     = p0k->im;

		    eterm  = tmp1k;
		    p0k->re = d1*eterm;
		    p0k->im = d2*eterm;
		  }

		  eterm    = tmp1k;
		  vterm    = tmp2k;
		  str2     = struct2k;
		  tmp1k = eterm*str2;
		  tmp2k = vterm*str2;
		}

                real ets2     = corner_fac*tmp1k;
                vterm    = 2.0*factor*tmp2k;
                energy  += ets2;
                real ets2vf   = corner_fac*vterm;
                virxx   += ets2vf*mhxk*mhxk - ets2;
                virxy   += ets2vf*mhxk*mhyk;
                virxz   += ets2vf*mhxk*mhzk;
                viryy   += ets2vf*mhyk*mhyk - ets2;
                viryz   += ets2vf*mhyk*mhzk;
                virzz   += ets2vf*mhzk*mhzk - ets2;
            }
        }
        else
        {
            /* We don't need to calculate the energy and the virial.
             *  In this case the triclinic overhead is small.
             */

            /* Two explicit loops to avoid a conditional inside the loop */
	    // NOTE: on gpu, keep the conditional. shouldn't be too bad?
            for (int kx = kxstart; kx < kxend; kx++)
            {
	        mx = kx < maxkx ? kx : (kx - nx);

                mhxk      = mx * rxx;
                mhyk      = mx * ryx + my * ryy;
                mhzk      = mx * rzx + my * rzy + mz * rzz;
                real m2k       = mhxk*mhxk + mhyk*mhyk + mhzk*mhzk;
                real denomk = bz*by*pme_bsp_mod_XX[kx];
                real tmp1k  = -factor*m2k;
                real tmp2k  = sqrt(factor*m2k);

		//calc_exponentials_lj_one(tmp1k, tmp2k, denomk); // r tmp2 d
		denomk = 1.0/denomk;
		tmp1k = exp(tmp1k);
		real mk = tmp2k;
		tmp2k = sqrt_M_PI_d*mk*erfcf(mk);

                m2k    = factor*m2k;
                real eterm  = -((1.0 - 2.0*m2k)*tmp1k
                           + 2.0*m2k*tmp2k);
                tmp1k = eterm*denomk;

		int gcount = (bLB ? 7 : 1);
		for (int ig = 0; ig < gcount; ++ig)
		{
		  //t_complex *p0;

		  t_complex *p0k = grid_v/*[ig]*/ + ig*grid_size + iy*local_size_ZZ*local_size_XX + iz*local_size_XX + (kx - kxstart);

		  real d1      = p0k->re;
		  real d2      = p0k->im;

		  eterm   = tmp1k;

		  p0k->re  = d1*eterm;
		  p0k->im  = d2*eterm;
                }
            }
        }
	energy_v[i] = energy;
	virial_v[0] = virxx;
	virial_v[1] = viryy;
	virial_v[2] = virzz;
	virial_v[3] = virxy;
	virial_v[4] = virxz;
	virial_v[5] = viryz;
    }
}
