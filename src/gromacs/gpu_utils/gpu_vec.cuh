/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_GPU_UTILS_GPU_VEC_CUH
#define GMX_GPU_UTILS_GPU_VEC_CUH

/* Note that because of the duplicate of ivec, this header (or an
 * OpenCL port of it) cannot be included in a translation unit that
 * also includes the normal vectypes.h */
#define XX 0 /* Defines for indexing in */
#define YY 1 /* vectors                 */
#define ZZ 2
#define DIM 3 /* Dimension of vectors    */
typedef int   ivec[DIM];
typedef float fvec[DIM];

/* maths operations */
/* imported from cpu versions in math/vec.h */
__forceinline__ __device__ void svmul_gpu(float a, const fvec v1, fvec v2)
{
    v2[XX] = a * v1[XX];
    v2[YY] = a * v1[YY];
    v2[ZZ] = a * v1[ZZ];
}


__forceinline__ __device__ void fvec_add_gpu(const fvec a, const fvec b, fvec c)
{
    float x, y, z;

    x = a[XX] + b[XX];
    y = a[YY] + b[YY];
    z = a[ZZ] + b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

__forceinline__ __device__ void ivec_add_gpu(const ivec a, const ivec b, ivec c)
{
    int x, y, z;

    x = a[XX] + b[XX];
    y = a[YY] + b[YY];
    z = a[ZZ] + b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

__forceinline__ __device__ void fvec_inc_atomic(fvec a, const fvec b)
{
    atomicAdd(&a[XX], b[XX]);
    atomicAdd(&a[YY], b[YY]);
    atomicAdd(&a[ZZ], b[ZZ]);
}

__forceinline__ __device__ void fvec_inc_gpu(fvec a, const fvec b)
{
    float x, y, z;

    x = a[XX] + b[XX];
    y = a[YY] + b[YY];
    z = a[ZZ] + b[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}

__forceinline__ __device__ void fvec_dec_atomic(fvec a, const fvec b)
{
    atomicAdd(&a[XX], -1.0f * b[XX]);
    atomicAdd(&a[YY], -1.0f * b[YY]);
    atomicAdd(&a[ZZ], -1.0f * b[ZZ]);
}

__forceinline__ __device__ void fvec_dec_gpu(fvec a, const fvec b)
{
    float x, y, z;

    x = a[XX] - b[XX];
    y = a[YY] - b[YY];
    z = a[ZZ] - b[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}

__forceinline__ __device__ void cprod_gpu(const fvec a, const fvec b, fvec c)
{
    c[XX] = a[YY] * b[ZZ] - a[ZZ] * b[YY];
    c[YY] = a[ZZ] * b[XX] - a[XX] * b[ZZ];
    c[ZZ] = a[XX] * b[YY] - a[YY] * b[XX];
}

__forceinline__ __device__ float iprod_gpu(const fvec a, const fvec b)
{
    return (a[XX] * b[XX] + a[YY] * b[YY] + a[ZZ] * b[ZZ]);
}

__forceinline__ __device__ float norm_gpu(const fvec a)
{
    return sqrt(iprod_gpu(a, a));
}

__forceinline__ __device__ float gmx_angle_gpu(const fvec a, const fvec b)
{
    fvec  w;
    float wlen, s;

    cprod_gpu(a, b, w);

    wlen = norm_gpu(w);
    s    = iprod_gpu(a, b);

    return atan2f(wlen, s); // requires float
}

__forceinline__ __device__ void clear_ivec_gpu(ivec a)
{
    a[XX] = 0;
    a[YY] = 0;
    a[ZZ] = 0;
}
__forceinline__ __device__ void fvec_sub_gpu(const fvec a, const fvec b, fvec c)
{
    float x, y, z;

    x = a[XX] - b[XX];
    y = a[YY] - b[YY];
    z = a[ZZ] - b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

__forceinline__ __device__ float norm2_gpu(const fvec a)
{
    return a[XX] * a[XX] + a[YY] * a[YY] + a[ZZ] * a[ZZ];
}

__forceinline__ __device__ void copy_fvec_gpu(const fvec a, fvec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

__forceinline__ __device__ void copy_ivec_gpu(const ivec a, ivec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

__forceinline__ __device__ float cos_angle_gpu(const fvec a, const fvec b)
{
    /*
     *                  ax*bx + ay*by + az*bz
     * cos-vec (a,b) =  ---------------------
     *                      ||a|| * ||b||
     */
    float cosval;
    int   m;
    float aa, bb, ip, ipa, ipb, ipab;

    ip = ipa = ipb = 0.0f;
    for (m = 0; (m < DIM); m++)
    {
        aa = a[m];
        bb = b[m];
        ip += aa * bb;
        ipa += aa * aa;
        ipb += bb * bb;
    }
    ipab = ipa * ipb;
    if (ipab > 0.0f)
    {
        cosval = ip * rsqrt(ipab);
    }
    else
    {
        cosval = 1.0f;
    }
    if (cosval > 1.0f)
    {
        return 1.0f;
    }
    if (cosval < -1.0f)
    {
        return -1.0f;
    }

    return cosval;
}


__device__ static inline void unitv_gpu(const fvec src, fvec dest)
{
    float linv;

    linv     = rsqrt(norm2_gpu(src));
    dest[XX] = linv * src[XX];
    dest[YY] = linv * src[YY];
    dest[ZZ] = linv * src[ZZ];
}

#endif
