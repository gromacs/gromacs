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

#ifndef GMX_GPU_UTILS_VECTYPE_OPS_CUDA_HIP_SHARED_DECLARATIONS_H
#define GMX_GPU_UTILS_VECTYPE_OPS_CUDA_HIP_SHARED_DECLARATIONS_H

#if !defined(__CUDACC__) && !defined(__HIPCC__)
#    error Including header specific for CUDA or HIP device code without compiling the file with the correct compiler
#endif

/**** float3 ****/
static __forceinline__ __host__ __device__ float3 make_float3(float s)
{
    return make_float3(s, s, s);
}
static __forceinline__ __host__ __device__ float3 make_float3(float4 a)
{
    return make_float3(a.x, a.y, a.z);
}
/**** float4 ****/
static __forceinline__ __host__ __device__ float4 make_float4(float s)
{
    return make_float4(s, s, s, s);
}
static __forceinline__ __host__ __device__ float4 make_float4(float3 a)
{
    return make_float4(a.x, a.y, a.z, 0.0F);
}
/**** atomics ****/
static __forceinline__ __device__ void atomicAdd(float3* addr, float3 val)
{
    atomicAdd(&addr->x, val.x);
    atomicAdd(&addr->y, val.y);
    atomicAdd(&addr->z, val.z);
}
// NOLINTNEXTLINE(google-runtime-references)
static __forceinline__ __device__ void atomicAdd(float3& a, const float3 b)
{
    atomicAdd(&a.x, b.x);
    atomicAdd(&a.y, b.y);
    atomicAdd(&a.z, b.z);
}

#endif /* GMX_GPU_UTILS_VECTYPE_OPS_CUDA_HIP_SHARED_DECLARATIONS_H */
