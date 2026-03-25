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

#ifndef GMX_GPU_UTILS_VECTYPE_OPS_SHARED_GPU_TRIG_MATH_H
#define GMX_GPU_UTILS_VECTYPE_OPS_SHARED_GPU_TRIG_MATH_H

#if !defined(__CUDACC__) && !defined(__HIPCC__)
#    error Including header specific for CUDA or HIP device code without compiling the file with the correct compiler
#endif

/**** float3 ****/
static __forceinline__ __host__ __device__ float gmxDeviceNorm(float3 a)
{
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}
static __forceinline__ __host__ __device__ float gmxDeviceNorm2(float3 a)
{
    return (a.x * a.x + a.y * a.y + a.z * a.z);
}

static __forceinline__ __host__ __device__ float gmxDeviceNorm(float4 a)
{
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z + a.w * a.w);
}

/* \brief Compute the scalar product of two vectors.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Scalar product.
 */
static __forceinline__ __device__ float gmxDeviceInternalProd(const float3 a, const float3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/* \brief Compute the vector product of two vectors.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Vector product.
 */
static __forceinline__ __device__ float3 gmxDeviceCrossProd(const float3 a, const float3 b)
{
    float3 c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}


/* \brief Compute the angle between two vectors.
 *
 * Uses atan( |axb| / a.b ) formula.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Angle between vectors in radians.
 */
static __forceinline__ __device__ float gmxDeviceAngle(const float3 a, const float3 b)
{
    float3 w = gmxDeviceCrossProd(a, b);

    float wlen = gmxDeviceNorm(w);
    float s    = gmxDeviceInternalProd(a, b);

    return atan2f(wlen, s); // requires float
}

static __forceinline__ __device__ float gmxDeviceSin(const float a)
{
    return sinf(a);
}

static __forceinline__ __device__ float gmxDeviceCos(const float a)
{
    return cosf(a);
}

/* \brief Cosine of an angle between two vectors.
 *
 * Computes cosine using the following formula:
 *
 *                  ax*bx + ay*by + az*bz
 * cos-vec (a,b) =  ---------------------
 *                      ||a|| * ||b||
 *
 * This function also makes sure that the cosine does not leave the [-1, 1]
 * interval, which can happen due to numerical errors.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Cosine between a and b.
 */
static __forceinline__ __device__ float gmxDeviceCosAngle(const float3 a, const float3 b)
{
    float cosval;

    float ipa  = gmxDeviceNorm2(a);
    float ipb  = gmxDeviceNorm2(b);
    float ip   = gmxDeviceInternalProd(a, b);
    float ipab = ipa * ipb;
    if (ipab > 0.0F)
    {
        cosval = ip * gmxDeviceRSqrt(ipab);
    }
    else
    {
        cosval = 1.0F;
    }
    if (cosval > 1.0F)
    {
        return 1.0F;
    }
    if (cosval < -1.0F)
    {
        return -1.0F;
    }

    return cosval;
}

#endif /* GMX_GPU_UTILS_VECTYPE_OPS_SHARED_CUDA_HIP_TRIG_MATH_H */
