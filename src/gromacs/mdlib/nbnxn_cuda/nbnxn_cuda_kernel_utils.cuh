/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,  by the GROMACS development team, led by
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

/* Note that floating-point constants in CUDA code should be suffixed
 * with f (e.g. 0.5f), to stop the compiler producing intermediate
 * code that is in double precision.
 */
#include "config.h"

#include "gromacs/gmxlib/cuda_tools/vectype_ops.cuh"

#ifndef NBNXN_CUDA_KERNEL_UTILS_CUH
#define NBNXN_CUDA_KERNEL_UTILS_CUH


#if defined HAVE_CUDA_TEXOBJ_SUPPORT && __CUDA_ARCH__ >= 300
/* Note: convenience macro, needs to be undef-ed at the end of the file. */
#define USE_TEXOBJ
#endif

#define WARP_SIZE_POW2_EXPONENT     (5)
#define CL_SIZE_POW2_EXPONENT       (3)  /* change this together with GPU_NS_CLUSTER_SIZE !*/
#define CL_SIZE_SQ                  (CL_SIZE * CL_SIZE)
#define FBUF_STRIDE                 (CL_SIZE_SQ)

#define ONE_SIXTH_F     0.16666667f
#define ONE_TWELVETH_F  0.08333333f


/*! Apply force switch,  force + energy version. */
static inline __device__
void calculate_force_switch_F(const  cu_nbparam_t nbparam,
                              float               c6,
                              float               c12,
                              float               inv_r,
                              float               r2,
                              float              *F_invr)
{
    float r, r_switch;

    /* force switch constants */
    float disp_shift_V2 = nbparam.dispersion_shift.c2;
    float disp_shift_V3 = nbparam.dispersion_shift.c3;
    float repu_shift_V2 = nbparam.repulsion_shift.c2;
    float repu_shift_V3 = nbparam.repulsion_shift.c3;

    r         = r2 * inv_r;
    r_switch  = r - nbparam.rvdw_switch;
    r_switch  = r_switch >= 0.0f ? r_switch : 0.0f;

    *F_invr  +=
        -c6*(disp_shift_V2 + disp_shift_V3*r_switch)*r_switch*r_switch*inv_r +
        c12*(-repu_shift_V2 + repu_shift_V3*r_switch)*r_switch*r_switch*inv_r;
}

/*! Apply force switch, force-only version. */
static inline __device__
void calculate_force_switch_F_E(const  cu_nbparam_t nbparam,
                                float               c6,
                                float               c12,
                                float               inv_r,
                                float               r2,
                                float              *F_invr,
                                float              *E_lj)
{
    float r, r_switch;

    /* force switch constants */
    float disp_shift_V2 = nbparam.dispersion_shift.c2;
    float disp_shift_V3 = nbparam.dispersion_shift.c3;
    float repu_shift_V2 = nbparam.repulsion_shift.c2;
    float repu_shift_V3 = nbparam.repulsion_shift.c3;

    float disp_shift_F2 = nbparam.dispersion_shift.c2/3;
    float disp_shift_F3 = nbparam.dispersion_shift.c3/4;
    float repu_shift_F2 = nbparam.repulsion_shift.c2/3;
    float repu_shift_F3 = nbparam.repulsion_shift.c3/4;

    r         = r2 * inv_r;
    r_switch  = r - nbparam.rvdw_switch;
    r_switch  = r_switch >= 0.0f ? r_switch : 0.0f;

    *F_invr  +=
        -c6*(disp_shift_V2 + disp_shift_V3*r_switch)*r_switch*r_switch*inv_r +
        c12*(-repu_shift_V2 + repu_shift_V3*r_switch)*r_switch*r_switch*inv_r;
    *E_lj    +=
        c6*(disp_shift_F2 + disp_shift_F3*r_switch)*r_switch*r_switch*r_switch -
        c12*(repu_shift_F2 + repu_shift_F3*r_switch)*r_switch*r_switch*r_switch;
}

/*! Apply potential switch, force-only version. */
static inline __device__
void calculate_potential_switch_F(const  cu_nbparam_t nbparam,
                                  float               c6,
                                  float               c12,
                                  float               inv_r,
                                  float               r2,
                                  float              *F_invr,
                                  float              *E_lj)
{
    float r, r_switch;
    float sw, dsw;

    /* potential switch constants */
    float switch_V3 = nbparam.vdw_switch.c3;
    float switch_V4 = nbparam.vdw_switch.c4;
    float switch_V5 = nbparam.vdw_switch.c5;
    float switch_F2 = 3*nbparam.vdw_switch.c3;
    float switch_F3 = 4*nbparam.vdw_switch.c4;
    float switch_F4 = 5*nbparam.vdw_switch.c5;

    r        = r2 * inv_r;
    r_switch = r - nbparam.rvdw_switch;

    /* Unlike in the F+E kernel, conditional is faster here */
    if (r_switch > 0.0f)
    {
        sw      = 1.0f + (switch_V3 + (switch_V4 + switch_V5*r_switch)*r_switch)*r_switch*r_switch*r_switch;
        dsw     = (switch_F2 + (switch_F3 + switch_F4*r_switch)*r_switch)*r_switch*r_switch;

        *F_invr = (*F_invr)*sw - inv_r*(*E_lj)*dsw;
    }
}

/*! Apply potential switch, force + energy version. */
static inline __device__
void calculate_potential_switch_F_E(const  cu_nbparam_t nbparam,
                                    float               c6,
                                    float               c12,
                                    float               inv_r,
                                    float               r2,
                                    float              *F_invr,
                                    float              *E_lj)
{
    float r, r_switch;
    float sw, dsw;

    /* potential switch constants */
    float switch_V3 = nbparam.vdw_switch.c3;
    float switch_V4 = nbparam.vdw_switch.c4;
    float switch_V5 = nbparam.vdw_switch.c5;
    float switch_F2 = 3*nbparam.vdw_switch.c3;
    float switch_F3 = 4*nbparam.vdw_switch.c4;
    float switch_F4 = 5*nbparam.vdw_switch.c5;

    r        = r2 * inv_r;
    r_switch = r - nbparam.rvdw_switch;
    r_switch = r_switch >= 0.0f ? r_switch : 0.0f;

    /* Unlike in the F-only kernel, masking is faster here */
    sw       = 1.0f + (switch_V3 + (switch_V4 + switch_V5*r_switch)*r_switch)*r_switch*r_switch*r_switch;
    dsw      = (switch_F2 + (switch_F3 + switch_F4*r_switch)*r_switch)*r_switch*r_switch;

    *F_invr  = (*F_invr)*sw - inv_r*(*E_lj)*dsw;
    *E_lj   *= sw;
}

/*! Calculate LJ-PME grid force contribution with
 *  geometric combination rule.
 */
static inline __device__
void calculate_lj_ewald_comb_geom_F(const cu_nbparam_t nbparam,
                                    int                typei,
                                    int                typej,
                                    float              r2,
                                    float              inv_r2,
                                    float              lje_coeff2,
                                    float              lje_coeff6_6,
                                    float             *F_invr)
{
    float c6grid, inv_r6_nm, cr2, expmcr2, poly;

#ifdef USE_TEXOBJ
    c6grid    = tex1Dfetch<float>(nbparam.nbfp_comb_texobj, 2*typei) * tex1Dfetch<float>(nbparam.nbfp_comb_texobj, 2*typej);
#else
    c6grid    = tex1Dfetch(nbfp_comb_texref, 2*typei) * tex1Dfetch(nbfp_comb_texref, 2*typej);
#endif /* USE_TEXOBJ */

    /* Recalculate inv_r6 without exclusion mask */
    inv_r6_nm = inv_r2*inv_r2*inv_r2;
    cr2       = lje_coeff2*r2;
    expmcr2   = expf(-cr2);
    poly      = 1.0f + cr2 + 0.5f*cr2*cr2;

    /* Subtract the grid force from the total LJ force */
    *F_invr  += c6grid*(inv_r6_nm - expmcr2*(inv_r6_nm*poly + lje_coeff6_6))*inv_r2;
}

/*! Calculate LJ-PME grid force + energy contribution with
 *  geometric combination rule.
 */
static inline __device__
void calculate_lj_ewald_comb_geom_F_E(const cu_nbparam_t nbparam,
                                      int                typei,
                                      int                typej,
                                      float              r2,
                                      float              inv_r2,
                                      float              lje_coeff2,
                                      float              lje_coeff6_6,
                                      float              int_bit,
                                      float             *F_invr,
                                      float             *E_lj)
{
    float c6grid, inv_r6_nm, cr2, expmcr2, poly, sh_mask;

#ifdef USE_TEXOBJ
    c6grid    = tex1Dfetch<float>(nbparam.nbfp_comb_texobj, 2*typei) * tex1Dfetch<float>(nbparam.nbfp_comb_texobj, 2*typej);
#else
    c6grid    = tex1Dfetch(nbfp_comb_texref, 2*typei) * tex1Dfetch(nbfp_comb_texref, 2*typej);
#endif /* USE_TEXOBJ */

    /* Recalculate inv_r6 without exclusion mask */
    inv_r6_nm = inv_r2*inv_r2*inv_r2;
    cr2       = lje_coeff2*r2;
    expmcr2   = expf(-cr2);
    poly      = 1.0f + cr2 + 0.5f*cr2*cr2;

    /* Subtract the grid force from the total LJ force */
    *F_invr  += c6grid*(inv_r6_nm - expmcr2*(inv_r6_nm*poly + lje_coeff6_6))*inv_r2;

    /* Shift should be applied only to real LJ pairs */
    sh_mask   = nbparam.sh_lj_ewald*int_bit;
    *E_lj    += ONE_SIXTH_F*c6grid*(inv_r6_nm*(1.0f - expmcr2*poly) + sh_mask);
}

/*! Calculate LJ-PME grid force + energy contribution (if E_lj != NULL) with
 *  Lorentz-Berthelot combination rule.
 *  We use a single F+E kernel with conditional because the performance impact
 *  of this is pretty small and LB on the CPU is anyway very slow.
 */
static inline __device__
void calculate_lj_ewald_comb_LB_F_E(const cu_nbparam_t nbparam,
                                    int                typei,
                                    int                typej,
                                    float              r2,
                                    float              inv_r2,
                                    float              lje_coeff2,
                                    float              lje_coeff6_6,
                                    float              int_bit,
                                    float             *F_invr,
                                    float             *E_lj)
{
    float c6grid, inv_r6_nm, cr2, expmcr2, poly;
    float sigma, sigma2, epsilon;

    /* sigma and epsilon are scaled to give 6*C6 */
#ifdef USE_TEXOBJ
    sigma   = tex1Dfetch<float>(nbparam.nbfp_comb_texobj, 2*typei    ) + tex1Dfetch<float>(nbparam.nbfp_comb_texobj, 2*typej    );
    epsilon = tex1Dfetch<float>(nbparam.nbfp_comb_texobj, 2*typei + 1) * tex1Dfetch<float>(nbparam.nbfp_comb_texobj, 2*typej + 1);
#else
    sigma   = tex1Dfetch(nbfp_comb_texref, 2*typei    ) + tex1Dfetch(nbfp_comb_texref, 2*typej    );
    epsilon = tex1Dfetch(nbfp_comb_texref, 2*typei + 1) * tex1Dfetch(nbfp_comb_texref, 2*typej + 1);
#endif /* USE_TEXOBJ */
    sigma2  = sigma*sigma;
    c6grid  = epsilon*sigma2*sigma2*sigma2;

    /* Recalculate inv_r6 without exclusion mask */
    inv_r6_nm = inv_r2*inv_r2*inv_r2;
    cr2       = lje_coeff2*r2;
    expmcr2   = expf(-cr2);
    poly      = 1.0f + cr2 + 0.5f*cr2*cr2;

    /* Subtract the grid force from the total LJ force */
    *F_invr  += c6grid*(inv_r6_nm - expmcr2*(inv_r6_nm*poly + lje_coeff6_6))*inv_r2;

    if (E_lj != NULL)
    {
        float sh_mask;

        /* Shift should be applied only to real LJ pairs */
        sh_mask   = nbparam.sh_lj_ewald*int_bit;
        *E_lj    += ONE_SIXTH_F*c6grid*(inv_r6_nm*(1.0f - expmcr2*poly) + sh_mask);
    }
}

/*! Interpolate Ewald coulomb force using the table through the tex_nbfp texture.
 *  Original idea: from the OpenMM project
 */
static inline __device__
float interpolate_coulomb_force_r(float r, float scale)
{
    float   normalized = scale * r;
    int     index      = (int) normalized;
    float   fract2     = normalized - index;
    float   fract1     = 1.0f - fract2;

    return fract1 * tex1Dfetch(coulomb_tab_texref, index)
           + fract2 * tex1Dfetch(coulomb_tab_texref, index + 1);
}

#ifdef HAVE_CUDA_TEXOBJ_SUPPORT
static inline __device__
float interpolate_coulomb_force_r(cudaTextureObject_t texobj_coulomb_tab,
                                  float r, float scale)
{
    float   normalized = scale * r;
    int     index      = (int) normalized;
    float   fract2     = normalized - index;
    float   fract1     = 1.0f - fract2;

    return fract1 * tex1Dfetch<float>(texobj_coulomb_tab, index) +
           fract2 * tex1Dfetch<float>(texobj_coulomb_tab, index + 1);
}
#endif


/*! Calculate analytical Ewald correction term. */
static inline __device__
float pmecorrF(float z2)
{
    const float FN6 = -1.7357322914161492954e-8f;
    const float FN5 = 1.4703624142580877519e-6f;
    const float FN4 = -0.000053401640219807709149f;
    const float FN3 = 0.0010054721316683106153f;
    const float FN2 = -0.019278317264888380590f;
    const float FN1 = 0.069670166153766424023f;
    const float FN0 = -0.75225204789749321333f;

    const float FD4 = 0.0011193462567257629232f;
    const float FD3 = 0.014866955030185295499f;
    const float FD2 = 0.11583842382862377919f;
    const float FD1 = 0.50736591960530292870f;
    const float FD0 = 1.0f;

    float       z4;
    float       polyFN0, polyFN1, polyFD0, polyFD1;

    z4          = z2*z2;

    polyFD0     = FD4*z4 + FD2;
    polyFD1     = FD3*z4 + FD1;
    polyFD0     = polyFD0*z4 + FD0;
    polyFD0     = polyFD1*z2 + polyFD0;

    polyFD0     = 1.0f/polyFD0;

    polyFN0     = FN6*z4 + FN4;
    polyFN1     = FN5*z4 + FN3;
    polyFN0     = polyFN0*z4 + FN2;
    polyFN1     = polyFN1*z4 + FN1;
    polyFN0     = polyFN0*z4 + FN0;
    polyFN0     = polyFN1*z2 + polyFN0;

    return polyFN0*polyFD0;
}

/*! Final j-force reduction; this generic implementation works with
 *  arbitrary array sizes.
 */
static inline __device__
void reduce_force_j_generic(float *f_buf, float3 *fout,
                            int tidxi, int tidxj, int aidx)
{
    if (tidxi < 3)
    {
        float f = 0.0f;
        for (int j = tidxj * CL_SIZE; j < (tidxj + 1) * CL_SIZE; j++)
        {
            f += f_buf[FBUF_STRIDE * tidxi + j];
        }

        atomicAdd((&fout[aidx].x)+tidxi, f);
    }
}

/*! Final j-force reduction; this implementation only with power of two
 *  array sizes and with sm >= 3.0
 */
#if __CUDA_ARCH__ >= 300
static inline __device__
void reduce_force_j_warp_shfl(float3 f, float3 *fout,
                              int tidxi, int aidx)
{
    f.x += __shfl_down(f.x, 1);
    f.y += __shfl_up  (f.y, 1);
    f.z += __shfl_down(f.z, 1);

    if (tidxi & 1)
    {
        f.x = f.y;
    }

    f.x += __shfl_down(f.x, 2);
    f.z += __shfl_up  (f.z, 2);

    if (tidxi & 2)
    {
        f.x = f.z;
    }

    f.x += __shfl_down(f.x, 4);

    if (tidxi < 3)
    {
        atomicAdd((&fout[aidx].x) + tidxi, f.x);
    }
}
#endif

/*! Final i-force reduction; this generic implementation works with
 *  arbitrary array sizes.
 * TODO: add the tidxi < 3 trick
 */
static inline __device__
void reduce_force_i_generic(float *f_buf, float3 *fout,
                            float *fshift_buf, bool bCalcFshift,
                            int tidxi, int tidxj, int aidx)
{
    if (tidxj < 3)
    {
        float f = 0.0f;
        for (int j = tidxi; j < CL_SIZE_SQ; j += CL_SIZE)
        {
            f += f_buf[tidxj * FBUF_STRIDE + j];
        }

        atomicAdd(&fout[aidx].x + tidxj, f);

        if (bCalcFshift)
        {
            *fshift_buf += f;
        }
    }
}

/*! Final i-force reduction; this implementation works only with power of two
 *  array sizes.
 */
static inline __device__
void reduce_force_i_pow2(volatile float *f_buf, float3 *fout,
                         float *fshift_buf, bool bCalcFshift,
                         int tidxi, int tidxj, int aidx)
{
    int     i, j;
    float   f;

    /* Reduce the initial CL_SIZE values for each i atom to half
     * every step by using CL_SIZE * i threads.
     * Can't just use i as loop variable because than nvcc refuses to unroll.
     */
    i = CL_SIZE/2;
    # pragma unroll 5
    for (j = CL_SIZE_POW2_EXPONENT - 1; j > 0; j--)
    {
        if (tidxj < i)
        {

            f_buf[                  tidxj * CL_SIZE + tidxi] += f_buf[                  (tidxj + i) * CL_SIZE + tidxi];
            f_buf[    FBUF_STRIDE + tidxj * CL_SIZE + tidxi] += f_buf[    FBUF_STRIDE + (tidxj + i) * CL_SIZE + tidxi];
            f_buf[2 * FBUF_STRIDE + tidxj * CL_SIZE + tidxi] += f_buf[2 * FBUF_STRIDE + (tidxj + i) * CL_SIZE + tidxi];
        }
        i >>= 1;
    }

    /* i == 1, last reduction step, writing to global mem */
    if (tidxj < 3)
    {
        /* tidxj*FBUF_STRIDE selects x, y or z */
        f = f_buf[tidxj * FBUF_STRIDE               + tidxi] +
            f_buf[tidxj * FBUF_STRIDE + i * CL_SIZE + tidxi];

        atomicAdd(&(fout[aidx].x) + tidxj, f);

        if (bCalcFshift)
        {
            *fshift_buf += f;
        }
    }

}

/*! Final i-force reduction wrapper; calls the generic or pow2 reduction depending
 *  on whether the size of the array to be reduced is power of two or not.
 */
static inline __device__
void reduce_force_i(float *f_buf, float3 *f,
                    float *fshift_buf, bool bCalcFshift,
                    int tidxi, int tidxj, int ai)
{
    if ((CL_SIZE & (CL_SIZE - 1)))
    {
        reduce_force_i_generic(f_buf, f, fshift_buf, bCalcFshift, tidxi, tidxj, ai);
    }
    else
    {
        reduce_force_i_pow2(f_buf, f, fshift_buf, bCalcFshift, tidxi, tidxj, ai);
    }
}

/*! Final i-force reduction; this implementation works only with power of two
 *  array sizes and with sm >= 3.0
 */
#if __CUDA_ARCH__ >= 300
static inline __device__
void reduce_force_i_warp_shfl(float3 fin, float3 *fout,
                              float *fshift_buf, bool bCalcFshift,
                              int tidxj, int aidx)
{
    fin.x += __shfl_down(fin.x, CL_SIZE);
    fin.y += __shfl_up  (fin.y, CL_SIZE);
    fin.z += __shfl_down(fin.z, CL_SIZE);

    if (tidxj & 1)
    {
        fin.x = fin.y;
    }

    fin.x += __shfl_down(fin.x, 2*CL_SIZE);
    fin.z += __shfl_up  (fin.z, 2*CL_SIZE);

    if (tidxj & 2)
    {
        fin.x = fin.z;
    }

    /* Threads 0,1,2 and 4,5,6 increment x,y,z for their warp */
    if ((tidxj & 3) < 3)
    {
        atomicAdd(&fout[aidx].x + (tidxj & ~4), fin.x);

        if (bCalcFshift)
        {
            *fshift_buf += fin.x;
        }
    }
}
#endif

/*! Energy reduction; this implementation works only with power of two
 *  array sizes.
 */
static inline __device__
void reduce_energy_pow2(volatile float *buf,
                        float *e_lj, float *e_el,
                        unsigned int tidx)
{
    int     i, j;
    float   e1, e2;

    i = WARP_SIZE/2;

    /* Can't just use i as loop variable because than nvcc refuses to unroll. */
# pragma unroll 10
    for (j = WARP_SIZE_POW2_EXPONENT - 1; j > 0; j--)
    {
        if (tidx < i)
        {
            buf[              tidx] += buf[              tidx + i];
            buf[FBUF_STRIDE + tidx] += buf[FBUF_STRIDE + tidx + i];
        }
        i >>= 1;
    }

    // TODO do this on two threads - requires e_lj and e_el to be stored on adjascent
    // memory locations to make sense
    /* last reduction step, writing to global mem */
    if (tidx == 0)
    {
        e1 = buf[              tidx] + buf[              tidx + i];
        e2 = buf[FBUF_STRIDE + tidx] + buf[FBUF_STRIDE + tidx + i];

        atomicAdd(e_lj, e1);
        atomicAdd(e_el, e2);
    }
}

/*! Energy reduction; this implementation works only with power of two
 *  array sizes and with sm >= 3.0
 */
#if __CUDA_ARCH__ >= 300
static inline __device__
void reduce_energy_warp_shfl(float E_lj, float E_el,
                             float *e_lj, float *e_el,
                             int tidx)
{
    int i, sh;

    sh = 1;
#pragma unroll 5
    for (i = 0; i < 5; i++)
    {
        E_lj += __shfl_down(E_lj, sh);
        E_el += __shfl_down(E_el, sh);
        sh   += sh;
    }

    /* The first thread in the warp writes the reduced energies */
    if (tidx == 0 || tidx == WARP_SIZE)
    {
        atomicAdd(e_lj, E_lj);
        atomicAdd(e_el, E_el);
    }
}
#endif /* __CUDA_ARCH__ */

#undef USE_TEXOBJ

#endif /* NBNXN_CUDA_KERNEL_UTILS_CUH */
