/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \internal \file
 *  \brief
 *  Utility constant and function declaration for the CUDA non-bonded kernels.
 *  This header should be included once at the top level, just before the
 *  kernels are included (has to be preceded by nbnxn_cuda_types.h).
 *
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \ingroup module_nbnxm
 */
#include <cassert>

/* Note that floating-point constants in CUDA code should be suffixed
 * with f (e.g. 0.5f), to stop the compiler producing intermediate
 * code that is in double precision.
 */

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cuda_kernel_utils.cuh"
#include "gromacs/gpu_utils/vectype_ops_cuda.h"
#include "gromacs/nbnxm/nbnxm_enums.h"
#include "gromacs/nbnxm/nbnxm_kernel_utils.h"

#include "nbnxm_cuda_types.h"

#ifndef NBNXM_CUDA_KERNEL_UTILS_CUH
#    define NBNXM_CUDA_KERNEL_UTILS_CUH

namespace gmx
{

/* Convenience defines */
//! CUDA device specific pairlist layout
static constexpr PairlistType sc_warpSize32Layout = PairlistType::Hierarchical8x8x8;

/*! \brief Log of the i and j cluster size.
 *  change this together with c_clSize !*/
static const int __device__ c_clusterSizeLog2 = gmx::StaticLog2<c_clusterSize>::value;
/*! \brief Stride in the force accumulation buffer */
static const int __device__ c_fbufStride = c_clusterSizeSq;

/*! Convert LJ sigma,epsilon parameters to C6,C12. */
static __forceinline__ __device__ void
convert_sigma_epsilon_to_c6_c12(const float sigma, const float epsilon, float* c6, float* c12)
{
    auto c6c12 = convertSigmaEpsilonToC6C12(sigma, epsilon);

    *c6  = c6c12.x;
    *c12 = c6c12.y;
}

/*! Apply force switch,  force + energy version. */
static __forceinline__ __device__ void
calculate_force_switch_F(const NBParamGpu nbparam, float c6, float c12, float inv_r, float r2, float* F_invr)
{
    float dummyValue = 0;
    ljForceSwitch<false>(
            nbparam.dispersion_shift, nbparam.repulsion_shift, nbparam.rvdw_switch, c6, c12, inv_r, r2, F_invr, &dummyValue);
}

/*! Apply force switch, force-only version. */
static __forceinline__ __device__ void calculate_force_switch_F_E(const NBParamGpu nbparam,
                                                                  float            c6,
                                                                  float            c12,
                                                                  float            inv_r,
                                                                  float            r2,
                                                                  float*           F_invr,
                                                                  float*           E_lj)
{
    ljForceSwitch<true>(
            nbparam.dispersion_shift, nbparam.repulsion_shift, nbparam.rvdw_switch, c6, c12, inv_r, r2, F_invr, E_lj);
}

/*! Apply potential switch, force-only version. */
static __forceinline__ __device__ void
calculate_potential_switch_F(const NBParamGpu& nbparam, float inv_r, float r2, float* F_invr, float* E_lj)
{
    ljPotentialSwitch<false>(nbparam.vdw_switch, nbparam.rvdw_switch, inv_r, r2, F_invr, E_lj);
}

/*! Apply potential switch, force + energy version. */
static __forceinline__ __device__ void
calculate_potential_switch_F_E(const NBParamGpu nbparam, float inv_r, float r2, float* F_invr, float* E_lj)
{
    ljPotentialSwitch<true>(nbparam.vdw_switch, nbparam.rvdw_switch, inv_r, r2, F_invr, E_lj);
}


/*! \brief Fetch C6 grid contribution coefficients and return the product of these.
 *
 *  Depending on what is supported, it fetches parameters either
 *  using direct load, texture objects, or texrefs.
 */
static __forceinline__ __device__ float calculate_lj_ewald_c6grid(const NBParamGpu nbparam, int typei, int typej)
{
#    if DISABLE_CUDA_TEXTURES
    float c6_i = LDG(&nbparam.nbfp_comb[typei]).x;
    float c6_j = LDG(&nbparam.nbfp_comb[typej]).x;
#    else
    float c6_i = tex1Dfetch<float2>(nbparam.nbfp_comb_texobj, typei).x;
    float c6_j = tex1Dfetch<float2>(nbparam.nbfp_comb_texobj, typej).x;
#    endif /* DISABLE_CUDA_TEXTURES */
    return c6_i * c6_j;
}


/*! Calculate LJ-PME grid force contribution with
 *  geometric combination rule.
 */
static __forceinline__ __device__ void calculate_lj_ewald_comb_geom_F(const NBParamGpu nbparam,
                                                                      int              typei,
                                                                      int              typej,
                                                                      float            r2,
                                                                      float            inv_r2,
                                                                      float            lje_coeff2,
                                                                      float            lje_coeff6_6,
                                                                      float*           F_invr)
{
    float c6grid, inv_r6_nm, cr2, expmcr2, poly;

    c6grid = calculate_lj_ewald_c6grid(nbparam, typei, typej);

    /* Recalculate inv_r6 without exclusion mask */
    inv_r6_nm = inv_r2 * inv_r2 * inv_r2;
    cr2       = lje_coeff2 * r2;
    expmcr2   = expf(-cr2);
    poly      = 1.0F + cr2 + 0.5F * cr2 * cr2;

    /* Subtract the grid force from the total LJ force */
    *F_invr += c6grid * (inv_r6_nm - expmcr2 * (inv_r6_nm * poly + lje_coeff6_6)) * inv_r2;
}


/*! Calculate LJ-PME grid force + energy contribution with
 *  geometric combination rule.
 */
static __forceinline__ __device__ void calculate_lj_ewald_comb_geom_F_E(const NBParamGpu nbparam,
                                                                        int              typei,
                                                                        int              typej,
                                                                        float            r2,
                                                                        float            inv_r2,
                                                                        float            lje_coeff2,
                                                                        float  lje_coeff6_6,
                                                                        float  int_bit,
                                                                        float* F_invr,
                                                                        float* E_lj)
{
    float c6grid, inv_r6_nm, cr2, expmcr2, poly, sh_mask;

    c6grid = calculate_lj_ewald_c6grid(nbparam, typei, typej);

    /* Recalculate inv_r6 without exclusion mask */
    inv_r6_nm = inv_r2 * inv_r2 * inv_r2;
    cr2       = lje_coeff2 * r2;
    expmcr2   = expf(-cr2);
    poly      = 1.0F + cr2 + 0.5F * cr2 * cr2;

    /* Subtract the grid force from the total LJ force */
    *F_invr += c6grid * (inv_r6_nm - expmcr2 * (inv_r6_nm * poly + lje_coeff6_6)) * inv_r2;

    /* Shift should be applied only to real LJ pairs */
    sh_mask = nbparam.sh_lj_ewald * int_bit;
    *E_lj += c_oneSixth * c6grid * (inv_r6_nm * (1.0F - expmcr2 * poly) + sh_mask);
}

/*! Fetch per-type LJ parameters.
 *
 *  Depending on what is supported, it fetches parameters either
 *  using direct load, texture objects, or texrefs.
 */
static __forceinline__ __device__ float2 fetch_nbfp_comb_c6_c12(const NBParamGpu nbparam, int type)
{
#    if DISABLE_CUDA_TEXTURES
    return LDG(&nbparam.nbfp_comb[type]);
#    else
    return tex1Dfetch<float2>(nbparam.nbfp_comb_texobj, type);
#    endif /* DISABLE_CUDA_TEXTURES */
}


/*! Calculate LJ-PME grid force + energy contribution (if E_lj != nullptr) with
 *  Lorentz-Berthelot combination rule.
 *  We use a single F+E kernel with conditional because the performance impact
 *  of this is pretty small and LB on the CPU is anyway very slow.
 */
static __forceinline__ __device__ void calculate_lj_ewald_comb_LB_F_E(const NBParamGpu nbparam,
                                                                      int              typei,
                                                                      int              typej,
                                                                      float            r2,
                                                                      float            inv_r2,
                                                                      float            lje_coeff2,
                                                                      float            lje_coeff6_6,
                                                                      float            int_bit,
                                                                      float*           F_invr,
                                                                      float*           E_lj)
{
    float c6grid, inv_r6_nm, cr2, expmcr2, poly;
    float sigma, sigma2, epsilon;

    /* sigma and epsilon are scaled to give 6*C6 */
    float2 c6c12_i = fetch_nbfp_comb_c6_c12(nbparam, typei);
    float2 c6c12_j = fetch_nbfp_comb_c6_c12(nbparam, typej);

    sigma   = c6c12_i.x + c6c12_j.x;
    epsilon = c6c12_i.y * c6c12_j.y;

    sigma2 = sigma * sigma;
    c6grid = epsilon * sigma2 * sigma2 * sigma2;

    /* Recalculate inv_r6 without exclusion mask */
    inv_r6_nm = inv_r2 * inv_r2 * inv_r2;
    cr2       = lje_coeff2 * r2;
    expmcr2   = expf(-cr2);
    poly      = 1.0F + cr2 + 0.5F * cr2 * cr2;

    /* Subtract the grid force from the total LJ force */
    *F_invr += c6grid * (inv_r6_nm - expmcr2 * (inv_r6_nm * poly + lje_coeff6_6)) * inv_r2;

    if (E_lj != nullptr)
    {
        float sh_mask;

        /* Shift should be applied only to real LJ pairs */
        sh_mask = nbparam.sh_lj_ewald * int_bit;
        *E_lj += c_oneSixth * c6grid * (inv_r6_nm * (1.0F - expmcr2 * poly) + sh_mask);
    }
}


/*! Fetch two consecutive values from the Ewald correction F*r table.
 *
 *  Depending on what is supported, it fetches parameters either
 *  using direct load, texture objects, or texrefs.
 */
static __forceinline__ __device__ float2 fetch_coulomb_force_r(const NBParamGpu nbparam, int index)
{
    float2 d;

#    if DISABLE_CUDA_TEXTURES
    /* Can't do 8-byte fetch because some of the addresses will be misaligned. */
    d.x = LDG(&nbparam.coulomb_tab[index]);
    d.y = LDG(&nbparam.coulomb_tab[index + 1]);
#    else
    d.x = tex1Dfetch<float>(nbparam.coulomb_tab_texobj, index);
    d.y = tex1Dfetch<float>(nbparam.coulomb_tab_texobj, index + 1);
#    endif // DISABLE_CUDA_TEXTURES

    return d;
}

/*! Interpolate Ewald coulomb force correction using the F*r table.
 */
static __forceinline__ __device__ float interpolate_coulomb_force_r(const NBParamGpu nbparam, float r)
{
    float normalized = nbparam.coulomb_tab_scale * r;
    int   index      = static_cast<int>(normalized);
    float fraction   = normalized - index;

    float2 d01 = fetch_coulomb_force_r(nbparam, index);

    return lerp(d01.x, d01.y, fraction);
}

/*! Fetch C6 and C12 from the parameter table.
 *
 *  Depending on what is supported, it fetches parameters either
 *  using direct load, texture objects, or texrefs.
 */
// NOLINTNEXTLINE(google-runtime-references)
static __forceinline__ __device__ void fetch_nbfp_c6_c12(float& c6, float& c12, const NBParamGpu nbparam, int baseIndex)
{
    float2 c6c12;
#    if DISABLE_CUDA_TEXTURES
    c6c12 = LDG(&nbparam.nbfp[baseIndex]);
#    else
    c6c12 = tex1Dfetch<float2>(nbparam.nbfp_texobj, baseIndex);
#    endif // DISABLE_CUDA_TEXTURES
    c6  = c6c12.x;
    c12 = c6c12.y;
}


/*! Final j-force reduction; this generic implementation works with
 *  arbitrary array sizes.
 */
static __forceinline__ __device__ void
reduce_force_j_generic(const float* f_buf, float3* fout, int tidxi, int tidxj, int aidx)
{
    if (tidxi < 3)
    {
        float f = 0.0F;
        for (int j = tidxj * c_clusterSize; j < (tidxj + 1) * c_clusterSize; j++)
        {
            f += f_buf[c_fbufStride * tidxi + j];
        }

        atomicAdd((&fout[aidx].x) + tidxi, f);
    }
}

/*! Final j-force reduction; this implementation only with power of two
 *  array sizes.
 */
static __forceinline__ __device__ void
reduce_force_j_warp_shfl(float3 f, float3* fout, int tidxi, int aidx, const unsigned int activemask)
{
    f.x += __shfl_down_sync(activemask, f.x, 1);
    f.y += __shfl_up_sync(activemask, f.y, 1);
    f.z += __shfl_down_sync(activemask, f.z, 1);

    if (tidxi & 1)
    {
        f.x = f.y;
    }

    f.x += __shfl_down_sync(activemask, f.x, 2);
    f.z += __shfl_up_sync(activemask, f.z, 2);

    if (tidxi & 2)
    {
        f.x = f.z;
    }

    f.x += __shfl_down_sync(activemask, f.x, 4);

    if (tidxi < 3)
    {
        atomicAdd((&fout[aidx].x) + tidxi, f.x);
    }
}

/*! Final i-force reduction; this generic implementation works with
 *  arbitrary array sizes.
 * TODO: add the tidxi < 3 trick
 */
static __forceinline__ __device__ void reduce_force_i_generic(const float* f_buf,
                                                              float3*      fout,
                                                              float*       fshift_buf,
                                                              bool         bCalcFshift,
                                                              int          tidxi,
                                                              int          tidxj,
                                                              int          aidx)
{
    if (tidxj < 3)
    {
        float f = 0.0F;
        for (int j = tidxi; j < c_clusterSizeSq; j += c_clusterSize)
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
static __forceinline__ __device__ void reduce_force_i_pow2(volatile float* f_buf,
                                                           float3*         fout,
                                                           float*          fshift_buf,
                                                           bool            bCalcFshift,
                                                           int             tidxi,
                                                           int             tidxj,
                                                           int             aidx)
{
    int   i, j;
    float f;

    static_assert(c_clusterSize == 1 << c_clusterSizeLog2);

    /* Reduce the initial c_clusterSize values for each i atom to half
     * every step by using c_clusterSize * i threads.
     * Can't just use i as loop variable because than nvcc refuses to unroll.
     */
    i = c_clusterSize / 2;
#    pragma unroll 5
    for (j = c_clusterSizeLog2 - 1; j > 0; j--)
    {
        if (tidxj < i)
        {

            f_buf[tidxj * c_clusterSize + tidxi] += f_buf[(tidxj + i) * c_clusterSize + tidxi];
            f_buf[c_fbufStride + tidxj * c_clusterSize + tidxi] +=
                    f_buf[c_fbufStride + (tidxj + i) * c_clusterSize + tidxi];
            f_buf[2 * c_fbufStride + tidxj * c_clusterSize + tidxi] +=
                    f_buf[2 * c_fbufStride + (tidxj + i) * c_clusterSize + tidxi];
        }
        i >>= 1;
    }

    /* i == 1, last reduction step, writing to global mem */
    if (tidxj < 3)
    {
        /* tidxj*c_fbufStride selects x, y or z */
        f = f_buf[tidxj * c_fbufStride + tidxi] + f_buf[tidxj * c_fbufStride + i * c_clusterSize + tidxi];

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
static __forceinline__ __device__ void
reduce_force_i(float* f_buf, float3* f, float* fshift_buf, bool bCalcFshift, int tidxi, int tidxj, int ai)
{
    if ((c_clusterSize & (c_clusterSize - 1)))
    {
        reduce_force_i_generic(f_buf, f, fshift_buf, bCalcFshift, tidxi, tidxj, ai);
    }
    else
    {
        reduce_force_i_pow2(f_buf, f, fshift_buf, bCalcFshift, tidxi, tidxj, ai);
    }
}

/*! Final i-force reduction; this implementation works only with power of two
 *  array sizes.
 */
static __forceinline__ __device__ void reduce_force_i_warp_shfl(float3             fin,
                                                                float3*            fout,
                                                                float*             fshift_buf,
                                                                bool               bCalcFshift,
                                                                int                tidxj,
                                                                int                aidx,
                                                                const unsigned int activemask)
{
    fin.x += __shfl_down_sync(activemask, fin.x, c_clusterSize);
    fin.y += __shfl_up_sync(activemask, fin.y, c_clusterSize);
    fin.z += __shfl_down_sync(activemask, fin.z, c_clusterSize);

    if (tidxj & 1)
    {
        fin.x = fin.y;
    }

    fin.x += __shfl_down_sync(activemask, fin.x, 2 * c_clusterSize);
    fin.z += __shfl_up_sync(activemask, fin.z, 2 * c_clusterSize);

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

/*! Energy reduction; this implementation works only with power of two
 *  array sizes.
 */
static __forceinline__ __device__ void
reduce_energy_pow2(volatile float* buf, float* e_lj, float* e_el, unsigned int tidx)
{
    float e1, e2;

    unsigned int i = warp_size / 2;

    /* Can't just use i as loop variable because than nvcc refuses to unroll. */
#    pragma unroll 10
    for (int j = warp_size_log2 - 1; j > 0; j--)
    {
        if (tidx < i)
        {
            buf[tidx] += buf[tidx + i];
            buf[c_fbufStride + tidx] += buf[c_fbufStride + tidx + i];
        }
        i >>= 1;
    }

    // TODO do this on two threads - requires e_lj and e_el to be stored on adjascent
    // memory locations to make sense
    /* last reduction step, writing to global mem */
    if (tidx == 0)
    {
        e1 = buf[tidx] + buf[tidx + i];
        e2 = buf[c_fbufStride + tidx] + buf[c_fbufStride + tidx + i];

        atomicAdd(e_lj, e1);
        atomicAdd(e_el, e2);
    }
}

/*! Energy reduction; this implementation works only with power of two
 *  array sizes.
 */
static __forceinline__ __device__ void
reduce_energy_warp_shfl(float E_lj, float E_el, float* e_lj, float* e_el, int tidx, const unsigned int activemask)
{
    int i, sh;

    sh = 1;
#    pragma unroll 5
    for (i = 0; i < 5; i++)
    {
        E_lj += __shfl_down_sync(activemask, E_lj, sh);
        E_el += __shfl_down_sync(activemask, E_el, sh);
        sh += sh;
    }

    /* The first thread in the warp writes the reduced energies */
    if (tidx == 0 || tidx == warp_size)
    {
        atomicAdd(e_lj, E_lj);
        atomicAdd(e_el, E_el);
    }
}

} // namespace gmx

#endif /* NBNXN_CUDA_KERNEL_UTILS_CUH */
