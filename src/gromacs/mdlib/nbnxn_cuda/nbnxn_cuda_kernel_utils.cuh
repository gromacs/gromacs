/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2016,2017, by the GROMACS development team, led by
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

/*! \internal \file
 *  \brief
 *  Utility constant and function declaration for the CUDA non-bonded kernels.
 *  This header should be included once at the top level, just before the
 *  kernels are included (has to be preceded by nbnxn_cuda_types.h).
 *
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \ingroup module_mdlib
 */
#include <assert.h>

/* Note that floating-point constants in CUDA code should be suffixed
 * with f (e.g. 0.5f), to stop the compiler producing intermediate
 * code that is in double precision.
 */

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cuda_kernel_utils.cuh"
#include "gromacs/gpu_utils/vectype_ops.cuh"

#include "nbnxn_cuda_types.h"

#ifndef NBNXN_CUDA_KERNEL_UTILS_CUH
#define NBNXN_CUDA_KERNEL_UTILS_CUH

/* Use texture objects if supported by the target hardware (and in host pass). */
#if GMX_PTX_ARCH >= 300 || GMX_PTX_ARCH == 0
/* Note: convenience macro, needs to be undef-ed at the end of the file. */
#define USE_TEXOBJ
#endif

/*! \brief Log of the i and j cluster size.
 *  change this together with c_clSize !*/
static const int          c_clSizeLog2  = 3;
/*! \brief Square of cluster size. */
static const int          c_clSizeSq    = c_clSize*c_clSize;
/*! \brief j-cluster size after split (4 in the current implementation). */
static const int          c_splitClSize = c_clSize/c_nbnxnGpuClusterpairSplit;
/*! \brief Stride in the force accumualation buffer */
static const int          c_fbufStride  = c_clSizeSq;
/*! \brief i-cluster interaction mask for a super-cluster with all c_numClPerSupercl=8 bits set */
static const unsigned     superClInteractionMask = ((1U << c_numClPerSupercl) - 1U);

static const float        c_oneSixth    = 0.16666667f;
static const float        c_oneTwelveth = 0.08333333f;


/*! Convert LJ sigma,epsilon parameters to C6,C12. */
static __forceinline__ __device__
void convert_sigma_epsilon_to_c6_c12(const float  sigma,
                                     const float  epsilon,
                                     float       *c6,
                                     float       *c12)
{
    float sigma2, sigma6;

    sigma2 = sigma * sigma;
    sigma6 = sigma2 *sigma2 * sigma2;
    *c6    = epsilon * sigma6;
    *c12   = *c6 * sigma6;
}

/*! Apply force switch,  force + energy version. */
static __forceinline__ __device__
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
static __forceinline__ __device__
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
static __forceinline__ __device__
void calculate_potential_switch_F(const  cu_nbparam_t nbparam,
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
static __forceinline__ __device__
void calculate_potential_switch_F_E(const  cu_nbparam_t nbparam,
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


/*! \brief Fetch C6 grid contribution coefficients and return the product of these.
 *
 *  Depending on what is supported, it fetches parameters either
 *  using direct load, texture objects, or texrefs.
 */
static __forceinline__ __device__
float calculate_lj_ewald_c6grid(const cu_nbparam_t nbparam,
                                int                typei,
                                int                typej)
{
#if DISABLE_CUDA_TEXTURES
    return LDG(&nbparam.nbfp_comb[2*typei]) * LDG(&nbparam.nbfp_comb[2*typej]);
#else
#ifdef USE_TEXOBJ
    return tex1Dfetch<float>(nbparam.nbfp_comb_texobj, 2*typei) * tex1Dfetch<float>(nbparam.nbfp_comb_texobj, 2*typej);
#else
    return tex1Dfetch(nbfp_comb_texref, 2*typei) * tex1Dfetch(nbfp_comb_texref, 2*typej);
#endif /* USE_TEXOBJ */
#endif /* DISABLE_CUDA_TEXTURES */
}


/*! Calculate LJ-PME grid force contribution with
 *  geometric combination rule.
 */
static __forceinline__ __device__
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

    c6grid = calculate_lj_ewald_c6grid(nbparam, typei, typej);

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
static __forceinline__ __device__
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

    c6grid = calculate_lj_ewald_c6grid(nbparam, typei, typej);

    /* Recalculate inv_r6 without exclusion mask */
    inv_r6_nm = inv_r2*inv_r2*inv_r2;
    cr2       = lje_coeff2*r2;
    expmcr2   = expf(-cr2);
    poly      = 1.0f + cr2 + 0.5f*cr2*cr2;

    /* Subtract the grid force from the total LJ force */
    *F_invr  += c6grid*(inv_r6_nm - expmcr2*(inv_r6_nm*poly + lje_coeff6_6))*inv_r2;

    /* Shift should be applied only to real LJ pairs */
    sh_mask   = nbparam.sh_lj_ewald*int_bit;
    *E_lj    += c_oneSixth*c6grid*(inv_r6_nm*(1.0f - expmcr2*poly) + sh_mask);
}

/*! Fetch per-type LJ parameters.
 *
 *  Depending on what is supported, it fetches parameters either
 *  using direct load, texture objects, or texrefs.
 */
static __forceinline__ __device__
float2 fetch_nbfp_comb_c6_c12(const cu_nbparam_t nbparam,
                              int                type)
{
    float2 c6c12;
#if DISABLE_CUDA_TEXTURES
    /* Force an 8-byte fetch to save a memory instruction. */
    float2 *nbfp_comb = (float2 *)nbparam.nbfp_comb;
    c6c12 = LDG(&nbfp_comb[type]);
#else
    /* NOTE: as we always do 8-byte aligned loads, we could
       fetch float2 here too just as above. */
#ifdef USE_TEXOBJ
    c6c12.x = tex1Dfetch<float>(nbparam.nbfp_comb_texobj, 2*type);
    c6c12.y = tex1Dfetch<float>(nbparam.nbfp_comb_texobj, 2*type + 1);
#else
    c6c12.x = tex1Dfetch(nbfp_comb_texref, 2*type);
    c6c12.y = tex1Dfetch(nbfp_comb_texref, 2*type + 1);
#endif /* USE_TEXOBJ */
#endif /* DISABLE_CUDA_TEXTURES */

    return c6c12;
}


/*! Calculate LJ-PME grid force + energy contribution (if E_lj != NULL) with
 *  Lorentz-Berthelot combination rule.
 *  We use a single F+E kernel with conditional because the performance impact
 *  of this is pretty small and LB on the CPU is anyway very slow.
 */
static __forceinline__ __device__
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
    float2 c6c12_i = fetch_nbfp_comb_c6_c12(nbparam, typei);
    float2 c6c12_j = fetch_nbfp_comb_c6_c12(nbparam, typej);

    sigma   = c6c12_i.x + c6c12_j.x;
    epsilon = c6c12_i.y * c6c12_j.y;

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
        *E_lj    += c_oneSixth*c6grid*(inv_r6_nm*(1.0f - expmcr2*poly) + sh_mask);
    }
}


/*! Fetch two consecutive values from the Ewald correction F*r table.
 *
 *  Depending on what is supported, it fetches parameters either
 *  using direct load, texture objects, or texrefs.
 */
static __forceinline__ __device__
float2 fetch_coulomb_force_r(const cu_nbparam_t nbparam,
                             int                index)
{
    float2 d;

#if DISABLE_CUDA_TEXTURES
    /* Can't do 8-byte fetch because some of the addresses will be misaligned. */
    d.x = LDG(&nbparam.coulomb_tab[index]);
    d.y = LDG(&nbparam.coulomb_tab[index + 1]);
#else
#ifdef USE_TEXOBJ
    d.x = tex1Dfetch<float>(nbparam.coulomb_tab_texobj, index);
    d.y = tex1Dfetch<float>(nbparam.coulomb_tab_texobj, index + 1);
#else
    d.x =  tex1Dfetch(coulomb_tab_texref, index);
    d.y =  tex1Dfetch(coulomb_tab_texref, index + 1);
#endif // USE_TEXOBJ
#endif // DISABLE_CUDA_TEXTURES

    return d;
}

/*! Linear interpolation using exactly two FMA operations.
 *
 *  Implements numeric equivalent of: (1-t)*d0 + t*d1
 *  Note that CUDA does not have fnms, otherwise we'd use
 *  fma(t, d1, fnms(t, d0, d0)
 *  but input modifiers are designed for this and are fast.
 */
template <typename T>
__forceinline__ __host__ __device__
T lerp(T d0, T d1, T t)
{
    return fma(t, d1, fma(-t, d0, d0));
}

/*! Interpolate Ewald coulomb force correction using the F*r table.
 */
static __forceinline__ __device__
float interpolate_coulomb_force_r(const cu_nbparam_t nbparam,
                                  float              r)
{
    float  normalized = nbparam.coulomb_tab_scale * r;
    int    index      = (int) normalized;
    float  fraction   = normalized - index;

    float2 d01 = fetch_coulomb_force_r(nbparam, index);

    return lerp(d01.x, d01.y, fraction);
}

/*! Fetch C6 and C12 from the parameter table.
 *
 *  Depending on what is supported, it fetches parameters either
 *  using direct load, texture objects, or texrefs.
 */
static __forceinline__ __device__
void fetch_nbfp_c6_c12(float               &c6,
                       float               &c12,
                       const cu_nbparam_t   nbparam,
                       int                  baseIndex)
{
#if DISABLE_CUDA_TEXTURES
    /* Force an 8-byte fetch to save a memory instruction. */
    float2 *nbfp = (float2 *)nbparam.nbfp;
    float2  c6c12;
    c6c12 = LDG(&nbfp[baseIndex]);
    c6    = c6c12.x;
    c12   = c6c12.y;
#else
    /* NOTE: as we always do 8-byte aligned loads, we could
       fetch float2 here too just as above. */
#ifdef USE_TEXOBJ
    c6  = tex1Dfetch<float>(nbparam.nbfp_texobj, 2*baseIndex);
    c12 = tex1Dfetch<float>(nbparam.nbfp_texobj, 2*baseIndex + 1);
#else
    c6  = tex1Dfetch(nbfp_texref, 2*baseIndex);
    c12 = tex1Dfetch(nbfp_texref, 2*baseIndex + 1);
#endif
#endif // DISABLE_CUDA_TEXTURES
}


/*! Calculate analytical Ewald correction term. */
static __forceinline__ __device__
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
static __forceinline__ __device__
void reduce_force_j_generic(float *f_buf, float3 *fout,
                            int tidxi, int tidxj, int aidx)
{
    if (tidxi < 3)
    {
        float f = 0.0f;
        for (int j = tidxj * c_clSize; j < (tidxj + 1) * c_clSize; j++)
        {
            f += f_buf[c_fbufStride * tidxi + j];
        }

        atomicAdd((&fout[aidx].x)+tidxi, f);
    }
}

/*! Final j-force reduction; this implementation only with power of two
 *  array sizes and with sm >= 3.0
 */
#if GMX_PTX_ARCH >= 300 || GMX_PTX_ARCH == 0
static __forceinline__ __device__
void reduce_force_j_warp_shfl(float3 f, float3 *fout,
                              int tidxi, int aidx,
                              const unsigned int activemask)
{
    f.x += gmx_shfl_down_sync(activemask, f.x, 1);
    f.y += gmx_shfl_up_sync  (activemask, f.y, 1);
    f.z += gmx_shfl_down_sync(activemask, f.z, 1);

    if (tidxi & 1)
    {
        f.x = f.y;
    }

    f.x += gmx_shfl_down_sync(activemask, f.x, 2);
    f.z += gmx_shfl_up_sync  (activemask, f.z, 2);

    if (tidxi & 2)
    {
        f.x = f.z;
    }

    f.x += gmx_shfl_down_sync(activemask, f.x, 4);

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
static __forceinline__ __device__
void reduce_force_i_generic(float *f_buf, float3 *fout,
                            float *fshift_buf, bool bCalcFshift,
                            int tidxi, int tidxj, int aidx)
{
    if (tidxj < 3)
    {
        float f = 0.0f;
        for (int j = tidxi; j < c_clSizeSq; j += c_clSize)
        {
            f += f_buf[tidxj * c_fbufStride + j];
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
static __forceinline__ __device__
void reduce_force_i_pow2(volatile float *f_buf, float3 *fout,
                         float *fshift_buf, bool bCalcFshift,
                         int tidxi, int tidxj, int aidx)
{
    int     i, j;
    float   f;

    assert(c_clSize == 1 << c_clSizeLog2);

    /* Reduce the initial c_clSize values for each i atom to half
     * every step by using c_clSize * i threads.
     * Can't just use i as loop variable because than nvcc refuses to unroll.
     */
    i = c_clSize/2;
#pragma unroll 5
    for (j = c_clSizeLog2 - 1; j > 0; j--)
    {
        if (tidxj < i)
        {

            f_buf[                   tidxj * c_clSize + tidxi] += f_buf[                   (tidxj + i) * c_clSize + tidxi];
            f_buf[    c_fbufStride + tidxj * c_clSize + tidxi] += f_buf[    c_fbufStride + (tidxj + i) * c_clSize + tidxi];
            f_buf[2 * c_fbufStride + tidxj * c_clSize + tidxi] += f_buf[2 * c_fbufStride + (tidxj + i) * c_clSize + tidxi];
        }
        i >>= 1;
    }

    /* i == 1, last reduction step, writing to global mem */
    if (tidxj < 3)
    {
        /* tidxj*c_fbufStride selects x, y or z */
        f = f_buf[tidxj * c_fbufStride               + tidxi] +
            f_buf[tidxj * c_fbufStride + i * c_clSize + tidxi];

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
static __forceinline__ __device__
void reduce_force_i(float *f_buf, float3 *f,
                    float *fshift_buf, bool bCalcFshift,
                    int tidxi, int tidxj, int ai)
{
    if ((c_clSize & (c_clSize - 1)))
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
#if GMX_PTX_ARCH >= 300 || GMX_PTX_ARCH == 0
static __forceinline__ __device__
void reduce_force_i_warp_shfl(float3 fin, float3 *fout,
                              float *fshift_buf, bool bCalcFshift,
                              int tidxj, int aidx,
                              const unsigned int activemask)
{
    fin.x += gmx_shfl_down_sync(activemask, fin.x, c_clSize);
    fin.y += gmx_shfl_up_sync  (activemask, fin.y, c_clSize);
    fin.z += gmx_shfl_down_sync(activemask, fin.z, c_clSize);

    if (tidxj & 1)
    {
        fin.x = fin.y;
    }

    fin.x += gmx_shfl_down_sync(activemask, fin.x, 2*c_clSize);
    fin.z += gmx_shfl_up_sync  (activemask, fin.z, 2*c_clSize);

    if (tidxj & 2)
    {
        fin.x = fin.z;
    }

    /* Threads 0,1,2 and 4,5,6 increment x,y,z for their warp */
    if ((tidxj & 3) < 3)
    {
        atomicAdd(&fout[aidx].x + (tidxj & 3), fin.x);

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
static __forceinline__ __device__
void reduce_energy_pow2(volatile float *buf,
                        float *e_lj, float *e_el,
                        unsigned int tidx)
{
    float        e1, e2;

    unsigned int i = warp_size/2;

    /* Can't just use i as loop variable because than nvcc refuses to unroll. */
#pragma unroll 10
    for (int j = warp_size_log2 - 1; j > 0; j--)
    {
        if (tidx < i)
        {
            buf[               tidx] += buf[               tidx + i];
            buf[c_fbufStride + tidx] += buf[c_fbufStride + tidx + i];
        }
        i >>= 1;
    }

    // TODO do this on two threads - requires e_lj and e_el to be stored on adjascent
    // memory locations to make sense
    /* last reduction step, writing to global mem */
    if (tidx == 0)
    {
        e1 = buf[               tidx] + buf[               tidx + i];
        e2 = buf[c_fbufStride + tidx] + buf[c_fbufStride + tidx + i];

        atomicAdd(e_lj, e1);
        atomicAdd(e_el, e2);
    }
}

/*! Energy reduction; this implementation works only with power of two
 *  array sizes and with sm >= 3.0
 */
#if GMX_PTX_ARCH >= 300 || GMX_PTX_ARCH == 0
static __forceinline__ __device__
void reduce_energy_warp_shfl(float E_lj, float E_el,
                             float *e_lj, float *e_el,
                             int tidx,
                             const unsigned int activemask)
{
    int i, sh;

    sh = 1;
#pragma unroll 5
    for (i = 0; i < 5; i++)
    {
        E_lj += gmx_shfl_down_sync(activemask, E_lj, sh);
        E_el += gmx_shfl_down_sync(activemask, E_el, sh);
        sh   += sh;
    }

    /* The first thread in the warp writes the reduced energies */
    if (tidx == 0 || tidx == warp_size)
    {
        atomicAdd(e_lj, E_lj);
        atomicAdd(e_el, E_el);
    }
}
#endif /* GMX_PTX_ARCH */

#undef USE_TEXOBJ

#endif /* NBNXN_CUDA_KERNEL_UTILS_CUH */
