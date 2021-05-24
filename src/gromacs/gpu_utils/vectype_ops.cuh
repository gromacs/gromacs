/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2015,2016,2019,2020,2021, by the GROMACS development team, led by
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

#ifndef VECTYPE_OPS_CUH
#define VECTYPE_OPS_CUH

/**** float3 ****/
__forceinline__ __host__ __device__ float3 make_float3(float s)
{
    return make_float3(s, s, s);
}
__forceinline__ __host__ __device__ float3 make_float3(float4 a)
{
    return make_float3(a.x, a.y, a.z);
}
__forceinline__ __host__ __device__ float3 operator-(const float3& a)
{
    return make_float3(-a.x, -a.y, -a.z);
}
__forceinline__ __host__ __device__ float3 operator+(float3 a, float3 b)
{
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
__forceinline__ __host__ __device__ float3 operator-(float3 a, float3 b)
{
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
__forceinline__ __host__ __device__ float3 operator*(float3 a, float k)
{
    return make_float3(k * a.x, k * a.y, k * a.z);
}
__forceinline__ __host__ __device__ float3 operator*(float k, float3 a)
{
    return make_float3(k * a.x, k * a.y, k * a.z);
}
// NOLINTNEXTLINE(google-runtime-references)
__forceinline__ __host__ __device__ void operator+=(float3& a, float3 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
}
// NOLINTNEXTLINE(google-runtime-references)
__forceinline__ __host__ __device__ void operator+=(float3& a, float4 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
}
// NOLINTNEXTLINE(google-runtime-references)
__forceinline__ __host__ __device__ void operator-=(float3& a, float3 b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
}
__forceinline__ __host__ __device__ float norm(float3 a)
{
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}
__forceinline__ __host__ __device__ float norm2(float3 a)
{
    return (a.x * a.x + a.y * a.y + a.z * a.z);
}
__forceinline__ __host__ __device__ float dist3(float3 a, float3 b)
{
    return norm(b - a);
}
__forceinline__ __host__ __device__ float3 operator*(float3 a, float3 b)
{
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}
// NOLINTNEXTLINE(google-runtime-references)
__forceinline__ __host__ __device__ void operator*=(float3& a, float3 b)
{
    a.x *= b.x;
    a.y *= b.y;
    a.z *= b.z;
}
// NOLINTNEXTLINE(google-runtime-references)
__forceinline__ __host__ __device__ void operator*=(float3& a, float b)
{
    a.x *= b;
    a.y *= b;
    a.z *= b;
}
__forceinline__ __device__ void atomicAdd(float3* addr, float3 val)
{
    atomicAdd(&addr->x, val.x);
    atomicAdd(&addr->y, val.y);
    atomicAdd(&addr->z, val.z);
}
/****************************************************************/

/**** float4 ****/
__forceinline__ __host__ __device__ float4 make_float4(float s)
{
    return make_float4(s, s, s, s);
}
__forceinline__ __host__ __device__ float4 make_float4(float3 a)
{
    return make_float4(a.x, a.y, a.z, 0.0F);
}
__forceinline__ __host__ __device__ float4 operator+(float4 a, float4 b)
{
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}
__forceinline__ __host__ __device__ float4 operator+(float4 a, float3 b)
{
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w);
}
__forceinline__ __host__ __device__ float4 operator-(float4 a, float4 b)
{
    return make_float4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}
__forceinline__ __host__ __device__ float4 operator*(float4 a, float k)
{
    return make_float4(k * a.x, k * a.y, k * a.z, k * a.w);
}
__forceinline__ __host__ __device__ void operator+=(float4& a, float4 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
}
// NOLINTNEXTLINE(google-runtime-references)
__forceinline__ __host__ __device__ void operator+=(float4& a, float3 b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
}
// NOLINTNEXTLINE(google-runtime-references)
__forceinline__ __host__ __device__ void operator-=(float4& a, float3 b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
}

__forceinline__ __host__ __device__ float norm(float4 a)
{
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z + a.w * a.w);
}

__forceinline__ __host__ __device__ float dist3(float4 a, float4 b)
{
    return norm(b - a);
}

/* \brief Compute the scalar product of two vectors.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Scalar product.
 */
__forceinline__ __device__ float iprod(const float3 a, const float3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/* \brief Compute the vector product of two vectors.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Vector product.
 */
__forceinline__ __device__ float3 cprod(const float3 a, const float3 b)
{
    float3 c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
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
__forceinline__ __device__ float cos_angle(const float3 a, const float3 b)
{
    float cosval;

    float ipa  = norm2(a);
    float ipb  = norm2(b);
    float ip   = iprod(a, b);
    float ipab = ipa * ipb;
    if (ipab > 0.0F)
    {
        cosval = ip * rsqrt(ipab);
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

/* \brief Compute the angle between two vectors.
 *
 * Uses atan( |axb| / a.b ) formula.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Angle between vectors in radians.
 */
__forceinline__ __device__ float gmx_angle(const float3 a, const float3 b)
{
    float3 w = cprod(a, b);

    float wlen = norm(w);
    float s    = iprod(a, b);

    return atan2f(wlen, s); // requires float
}

/* \brief Atomically add components of the vector.
 *
 * Executes atomicAdd one-by-one on all components of the float3 vector.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 */
// NOLINTNEXTLINE(google-runtime-references)
__forceinline__ __device__ void atomicAdd(float3& a, const float3 b)
{
    atomicAdd(&a.x, b.x);
    atomicAdd(&a.y, b.y);
    atomicAdd(&a.z, b.z);
}

#endif /* VECTYPE_OPS_CUH */
