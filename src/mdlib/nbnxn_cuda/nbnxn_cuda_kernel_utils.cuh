/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef NBNXN_CUDA_KERNEL_UTILS_CUH
#define NBNXN_CUDA_KERNEL_UTILS_CUH

#define CLUSTER_SIZE_POW2_EXPONENT (3)  /* change this together with GPU_NS_CLUSTER_SIZE !*/
#define CLUSTER_SIZE_2          (CLUSTER_SIZE * CLUSTER_SIZE)
#define STRIDE_DIM              (CLUSTER_SIZE_2)
#define STRIDE_SI               (3*STRIDE_DIM)

/*! texture ref for nonbonded parameters; bound to cu_nbparam_t.nbfp*/
texture<float, 1, cudaReadModeElementType> tex_nbfp;

/*! texture ref for Ewald coulomb force table; bound to cu_nbparam_t.coulomb_tab */
texture<float, 1, cudaReadModeElementType> tex_coulomb_tab;

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

    return  fract1 * tex1Dfetch(tex_coulomb_tab, index) 
            + fract2 * tex1Dfetch(tex_coulomb_tab, index + 1);
}

/*! Final j-force reduction; this generic implementation works with 
 *  arbitrary array sizes. 
 */
static inline __device__ 
void reduce_force_j_generic(float *f_buf, float4 *fout,
                            int tidxi, int tidxj, int aidx)
{
    if (tidxi == 0)
    {
        float4 f = make_float4(0.0f);
        for (int j = tidxj * CLUSTER_SIZE; j < (tidxj + 1) * CLUSTER_SIZE; j++)
        {
            f.x += f_buf[                 j];
            f.y += f_buf[    STRIDE_DIM + j];
            f.z += f_buf[2 * STRIDE_DIM + j];
        }

        atomicAdd3(&fout[aidx], f);
    }
}

/*! Final i-force reduction; this generic implementation works with 
 *  arbitrary array sizes. 
 */
static inline __device__ 
void reduce_force_i_generic(float *f_buf, float4 *fout,
                            float3 *fshift_buf, gmx_bool calc_fshift,
                            int tidxi, int tidxj, int aidx)
{
    if (tidxj == 0)
    {
        float4 f = make_float4(0.0f);
        for (int j = tidxi; j < CLUSTER_SIZE_2; j += CLUSTER_SIZE)
        {
            f.x += f_buf[                 j];
            f.y += f_buf[    STRIDE_DIM + j];
            f.z += f_buf[2 * STRIDE_DIM + j];
        }

        atomicAdd3(&fout[aidx], f);

        if (calc_fshift)
        {
            *fshift_buf += f;
        }
    }
}

/*! Final i-force reduction; this implementation works only with power of two
 *  array sizes. 
 */
static inline __device__ 
void reduce_force_i_pow2(volatile float *f_buf, float4 *fout,
                         float3 *fshift_buf, gmx_bool calc_fshift,
                         int tidxi, int tidxj, int aidx)
{
    int     i, j; 
    float4  f = make_float4(0.0f);

    /* Reduce the initial CLUSTER_SIZE values for each i atom to half
       every step by using CLUSTER_SIZE * i threads. */
    i = CLUSTER_SIZE/2;
    # pragma unroll 5
    for (j = CLUSTER_SIZE_POW2_EXPONENT - 1; j > 0; j--)
    {
        if (tidxj < i)
        {

            f_buf[                 tidxj * CLUSTER_SIZE + tidxi] += f_buf[                 (tidxj + i) * CLUSTER_SIZE + tidxi];
            f_buf[    STRIDE_DIM + tidxj * CLUSTER_SIZE + tidxi] += f_buf[    STRIDE_DIM + (tidxj + i) * CLUSTER_SIZE + tidxi];
            f_buf[2 * STRIDE_DIM + tidxj * CLUSTER_SIZE + tidxi] += f_buf[2 * STRIDE_DIM + (tidxj + i) * CLUSTER_SIZE + tidxi];
        }
        i >>= 1;
    }

    /* i == 1, last reduction step, writing to global mem */
    if (tidxj == 0)
    {
        f.x = f_buf[                 tidxj * CLUSTER_SIZE + tidxi] + f_buf[                 (tidxj + i) * CLUSTER_SIZE + tidxi];
        f.y = f_buf[    STRIDE_DIM + tidxj * CLUSTER_SIZE + tidxi] + f_buf[    STRIDE_DIM + (tidxj + i) * CLUSTER_SIZE + tidxi];
        f.z = f_buf[2 * STRIDE_DIM + tidxj * CLUSTER_SIZE + tidxi] + f_buf[2 * STRIDE_DIM + (tidxj + i) * CLUSTER_SIZE + tidxi];

        atomicAdd3(&fout[aidx], f);

        if (calc_fshift)
        {
            *fshift_buf += f;
        }
    }
}

/*! Final i-force reduction wrapper; calls the generic or pow2 reduction depending
 *  on whether the size of the array to be reduced is power of two or not.
 */
static inline __device__ 
void reduce_force_i(float *f_buf, float4 *f,
                    float3 *fshift_buf, gmx_bool calc_fshift,
                    int tidxi, int tidxj, int ai)
{
    if ((CLUSTER_SIZE & (CLUSTER_SIZE - 1)))
    {
        reduce_force_i_generic(f_buf, f, fshift_buf, calc_fshift, tidxi, tidxj, ai);
    }
    else
    {
        reduce_force_i_pow2(f_buf, f, fshift_buf, calc_fshift, tidxi, tidxj, ai);
    }
}

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

    i = CLUSTER_SIZE_2/2;

# pragma unroll 10
    for (j = 2 * CLUSTER_SIZE_POW2_EXPONENT - 1; j > 0; j--)
    {
        if (tidx < i)
        {
            buf[             tidx] += buf[             tidx + i];
            buf[STRIDE_DIM + tidx] += buf[STRIDE_DIM + tidx + i];
        }
        i >>= 1;
    }

    /* last reduction step, writing to global mem */
    if (tidx == 0)
    {
        e1 = buf[             tidx] + buf[             tidx + i];
        e2 = buf[STRIDE_DIM + tidx] + buf[STRIDE_DIM + tidx + i];

        atomicAdd(e_lj, e1);
        atomicAdd(e_el, e2); 
    }
}

/*********************************************************************************/
/* Old stuff  */
#if 0
/* This function was used to calculate the Ewald coulomb force, but it's not 
 * used as it's much slower than tabulated f interpolation with texture mem.
 */
inline __device__ float 
coulomb(float q1,
        float q2,
        float r2,
        float inv_r,
        float inv_r2,
        float beta,
        float erfc_tab_scale)
{
    float x      = r2 * inv_r * beta;
    float x2     = x * x; 
    // float inv_x2 = inv_r2 / (beta * beta); 
    float res    =
        q1 * q2 * (erfc(x) * inv_r + beta * exp(-x2)) * inv_r2;
    return res;
}
#endif 

#endif /* NBNXN_CUDA_KERNEL_UTILS_CUH */
