/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#include "../../gmxlib/cuda_tools/vectype_ops.cuh"

#ifndef NBNXN_CUDA_KERNEL_UTILS_CUH
#define NBNXN_CUDA_KERNEL_UTILS_CUH

#define WARP_SIZE_POW2_EXPONENT     (5)
#define CL_SIZE_POW2_EXPONENT       (3)  /* change this together with GPU_NS_CLUSTER_SIZE !*/
#define CL_SIZE_SQ                  (CL_SIZE * CL_SIZE)
#define FBUF_STRIDE                 (CL_SIZE_SQ)

/*! i-cluster interaction mask for a super-cluster with all NCL_PER_SUPERCL bits set */
const unsigned supercl_interaction_mask = ((1U << NCL_PER_SUPERCL) - 1U);

/*! Interpolate Ewald coulomb force using the table through the tex_nbfp texture.
 *  Original idea: OpenMM
 */
static inline __device__
float interpolate_coulomb_force_r(float r, float scale)
{
    float   normalized = scale * r;
    int     index = (int) normalized;
    float   fract2 = normalized - index;
    float   fract1 = 1.0f - fract2;

    return  fract1 * tex1Dfetch(coulomb_tab_texref, index)
            + fract2 * tex1Dfetch(coulomb_tab_texref, index + 1);
}

#ifdef TEXOBJ_SUPPORTED
static inline __device__
float interpolate_coulomb_force_r(cudaTextureObject_t texobj_coulomb_tab,
                                  float r, float scale)
{
    float   normalized = scale * r;
    int     index = (int) normalized;
    float   fract2 = normalized - index;
    float   fract1 = 1.0f - fract2;

    return  fract1 * tex1Dfetch<float>(texobj_coulomb_tab, index) +
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
    float       polyFN0,polyFN1,polyFD0,polyFD1;

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
    if (tidxi == 0)
    {
        float3 f = make_float3(0.0f);
        for (int j = tidxj * CL_SIZE; j < (tidxj + 1) * CL_SIZE; j++)
        {
            f.x += f_buf[                  j];
            f.y += f_buf[    FBUF_STRIDE + j];
            f.z += f_buf[2 * FBUF_STRIDE + j];
        }

        atomicAdd(&fout[aidx], f);
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
    int i;

#pragma unroll 3
    for (i = 0; i < 3; i++)
    {
        f.x += __shfl_down(f.x, 1<<i);
        f.y += __shfl_down(f.y, 1<<i);
        f.z += __shfl_down(f.z, 1<<i);
    }

    /* Write the reduced j-force on one thread for each j */
    if (tidxi == 0)
    {
        atomicAdd(&fout[aidx], f);
    }
}
#endif

/*! Final i-force reduction; this generic implementation works with
 *  arbitrary array sizes.
 */
static inline __device__
void reduce_force_i_generic(float *f_buf, float3 *fout,
                            float3 *fshift_buf, bool bCalcFshift,
                            int tidxi, int tidxj, int aidx)
{
    if (tidxj == 0)
    {
        float3 f = make_float3(0.0f);
        for (int j = tidxi; j < CL_SIZE_SQ; j += CL_SIZE)
        {
            f.x += f_buf[                  j];
            f.y += f_buf[    FBUF_STRIDE + j];
            f.z += f_buf[2 * FBUF_STRIDE + j];
        }

        atomicAdd(&fout[aidx], f);

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
                         float3 *fshift_buf, bool bCalcFshift,
                         int tidxi, int tidxj, int aidx)
{
    int     i, j;
    float3  f = make_float3(0.0f);

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
    if (tidxj == 0)
    {
        f.x = f_buf[                  tidxj * CL_SIZE + tidxi] + f_buf[                  (tidxj + i) * CL_SIZE + tidxi];
        f.y = f_buf[    FBUF_STRIDE + tidxj * CL_SIZE + tidxi] + f_buf[    FBUF_STRIDE + (tidxj + i) * CL_SIZE + tidxi];
        f.z = f_buf[2 * FBUF_STRIDE + tidxj * CL_SIZE + tidxi] + f_buf[2 * FBUF_STRIDE + (tidxj + i) * CL_SIZE + tidxi];

        atomicAdd(&fout[aidx], f);

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
                    float3 *fshift_buf, bool bCalcFshift,
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
                              float3 *fshift_buf, bool bCalcFshift,
                              int tidxj, int aidx)
{
    int j;

#pragma unroll 2
    for (j = 0; j < 2; j++)
    {
        fin.x += __shfl_down(fin.x,  CL_SIZE<<j);
        fin.y += __shfl_down(fin.y,  CL_SIZE<<j);
        fin.z += __shfl_down(fin.z,  CL_SIZE<<j);
    }

    /* The first thread in the warp writes the reduced force */
    if (tidxj == 0 || tidxj == 4)
    {
        atomicAdd(&fout[aidx], fin);

        if (bCalcFshift)
        {
            fshift_buf->x += fin.x;
            fshift_buf->y += fin.y;
            fshift_buf->z += fin.z;
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
        E_lj += __shfl_down(E_lj,sh);
        E_el += __shfl_down(E_el,sh);
        sh += sh;
    }

    /* The first thread in the warp writes the reduced energies */
    if (tidx == 0 || tidx == WARP_SIZE)
    {
        atomicAdd(e_lj,E_lj);
        atomicAdd(e_el,E_el);
    }
}
#endif /* __CUDA_ARCH__ */

#endif /* NBNXN_CUDA_KERNEL_UTILS_CUH */
