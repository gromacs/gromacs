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

#ifndef GMX_GPU_UTILS_VECTYPE_OPS_CUDA_BASE_MATH_H
#define GMX_GPU_UTILS_VECTYPE_OPS_CUDA_BASE_MATH_H

#if !defined(__CUDACC__)
#    error Including header specific for CUDA device code without compiling the file with the correct compiler
#endif

static __forceinline__ __host__ __device__ float3 operator-(const float3& a)
{
    return make_float3(-a.x, -a.y, -a.z);
}
static __forceinline__ __host__ __device__ float3 operator+(float3 a, float3 b)
{
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
static __forceinline__ __host__ __device__ float3 operator-(float3 a, float3 b)
{
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
static __forceinline__ __host__ __device__ float3 operator*(float3 a, float k)
{
    return make_float3(k * a.x, k * a.y, k * a.z);
}
static __forceinline__ __host__ __device__ float3 operator*(float k, float3 a)
{
    return make_float3(k * a.x, k * a.y, k * a.z);
}
// NOLINTNEXTLINE(google-runtime-references)
static __forceinline__ __host__ __device__ void operator+=(float3& a, float3 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
}
// NOLINTNEXTLINE(google-runtime-references)
static __forceinline__ __host__ __device__ void operator+=(float3& a, float4 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
}
// NOLINTNEXTLINE(google-runtime-references)
static __forceinline__ __host__ __device__ void operator-=(float3& a, float3 b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
}
static __forceinline__ __host__ __device__ float3 operator*(float3 a, float3 b)
{
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}
// NOLINTNEXTLINE(google-runtime-references)
static __forceinline__ __host__ __device__ void operator*=(float3& a, float3 b)
{
    a.x *= b.x;
    a.y *= b.y;
    a.z *= b.z;
}
// NOLINTNEXTLINE(google-runtime-references)
static __forceinline__ __host__ __device__ void operator*=(float3& a, float b)
{
    a.x *= b;
    a.y *= b;
    a.z *= b;
}
static __forceinline__ __host__ __device__ float4 operator+(float4 a, float4 b)
{
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}
static __forceinline__ __host__ __device__ float4 operator+(float4 a, float3 b)
{
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w);
}
static __forceinline__ __host__ __device__ float4 operator-(float4 a, float4 b)
{
    return make_float4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}
static __forceinline__ __host__ __device__ float4 operator*(float4 a, float k)
{
    return make_float4(k * a.x, k * a.y, k * a.z, k * a.w);
}
static __forceinline__ __host__ __device__ void operator+=(float4& a, float4 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
}
// NOLINTNEXTLINE(google-runtime-references)
static __forceinline__ __host__ __device__ void operator+=(float4& a, float3 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
}
// NOLINTNEXTLINE(google-runtime-references)
static __forceinline__ __host__ __device__ void operator-=(float4& a, float3 b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
}
// NOLINTNEXTLINE(google-runtime-references)
static __forceinline__ __device__ float gmxDeviceRSqrt(float& input)
{
    return rsqrt(input);
}

#endif /* GMX_GPU_UTILS_VECTYPE_OPS_CUDA_BASE_MATH_H */
