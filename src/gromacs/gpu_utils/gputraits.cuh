/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
#ifndef GMX_GPU_UTILS_GPUTRAITS_CUH
#define GMX_GPU_UTILS_GPUTRAITS_CUH

/*! \libinternal \file
 *  \brief Declares the CUDA type traits.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_gpu_utils
 */
#include <cuda_runtime.h>

#include "gromacs/utility/vectypes.h"

#define GMX_HOST_ATTRIBUTE __host__
#define GMX_DEVICE_ATTRIBUTE __device__
#define GMX_HOSTDEVICE_ATTRIBUTE GMX_HOST_ATTRIBUTE GMX_DEVICE_ATTRIBUTE
#if !defined(NDEBUG)
#    define GMX_DEVICE_ASSERT(condition) assert(condition)
#else
#    define GMX_DEVICE_ASSERT(condition)
#endif

//! Device texture for fast read-only data fetching
using DeviceTexture = cudaTextureObject_t;

//! \brief Single GPU call timing event - meaningless in CUDA
using CommandEvent = void;

//! Convenience alias for 2-wide float
using Float2 = float2;

//! Convenience alias for 3-wide float
using Float3 = gmx::RVec;

//! Convenience alias for 3-wide float in shared device kernels
using DeviceFloat3 = float3;

//! Convenience alias for 4-wide float
using Float4 = float4;

//! Convenience alias for 4-wide float in shared device kernels.
struct DeviceFloat4
{
    __device__ __forceinline__ DeviceFloat4(float4 in) : storage_(in) {}

    template<typename Index>
    __device__ __forceinline__ float operator[](Index i) const
    {
        switch (i)
        {
            case 0: return storage_.x;
            case 1: return storage_.y;
            case 2: return storage_.z;
            default: GMX_DEVICE_ASSERT(i == 3); return storage_.w;
        }
    }
    __device__ __forceinline__ operator float4() const { return storage_; }

    alignas(16) float4 storage_;
};

//! Convenience alias for int3 in shared device kernels
struct DeviceInt3
{
    __device__ __forceinline__ DeviceInt3(int x, int y, int z) : storage_{ x, y, z } {}
    template<typename Index>
    __device__ __forceinline__ float operator[](Index i) const
    {
        switch (i)
        {
            case 0: return storage_.x;
            case 1: return storage_.y;
            default: GMX_DEVICE_ASSERT(i == 2); return storage_.z;
        }
    }
    int3 storage_;
};

//! Convenience alias for int4 in shared device kernels
struct DeviceInt4
{
    __device__ __forceinline__ DeviceInt4(int4 in) : storage_(in) {}
    template<typename Index>
    __device__ __forceinline__ float operator[](Index i) const
    {
        switch (i)
        {
            case 0: return storage_.x;
            case 1: return storage_.y;
            case 2: return storage_.z;
            default: GMX_DEVICE_ASSERT(i == 3); return storage_.w;
        }
    }
    int4 storage_;
};

static __device__ __forceinline__ DeviceInt4 loadInt4(const int* input, const int index)
{
    return { *(reinterpret_cast<const int4*>(input + 4 * index)) };
}

//! Convenience alias for global device memory
template<typename T>
using DeviceGlobalPtr = T*;
//! Convenience alias for local device memory
template<typename T>
using DeviceLocalPtr = T*;
//! Convenience alias for private device memory
template<typename T>
using DevicePrivatePtr = T*;

/*! \internal \brief
 * GPU kernels scheduling description. This is same in OpenCL/CUDA.
 * Provides reasonable defaults, one typically only needs to set the GPU stream
 * and non-1 work sizes.
 */
struct KernelLaunchConfig
{
    //! Block counts
    size_t gridSize[3] = { 1, 1, 1 };
    //! Per-block thread counts
    size_t blockSize[3] = { 1, 1, 1 };
    //! Shared memory size in bytes
    size_t sharedMemorySize = 0;
};

//! Sets whether device code can use arrays that are embedded in structs.
static constexpr bool c_canEmbedBuffers = true;

#endif
