/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
#ifndef GMX_GPU_UTILS_CUDA_KERNEL_UTILS_CUH
#define GMX_GPU_UTILS_CUDA_KERNEL_UTILS_CUH

/*! \file
 *  \brief CUDA device util functions (callable from GPU kernels).
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 */

#include "config.h"

#include <cassert>

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/utility/basedefinitions.h"

/*! \brief Bitmask corresponding to all threads active in a warp.
 *  NOTE that here too we assume 32-wide warps,
 * same as warp_size = 32 in the included common header.
 */
static const unsigned int c_fullWarpMask = 0xffffffff;

/* Below are backward-compatibility wrappers for CUDA 9 warp-wide intrinsics. */

/*! \brief Compatibility wrapper around the CUDA __syncwarp() instrinsic.  */
static __forceinline__ __device__
void gmx_syncwarp(const unsigned int activeMask = c_fullWarpMask)
{
#if GMX_CUDA_VERSION < 9000
    /* no sync needed on pre-Volta. */
    GMX_UNUSED_VALUE(activeMask);
#else
    __syncwarp(activeMask);
#endif
}

/*! \brief Compatibility wrapper around the CUDA __ballot()/__ballot_sync() instrinsic.  */
static __forceinline__ __device__
unsigned int gmx_ballot_sync(const unsigned int activeMask,
                             const int          pred)
{
#if GMX_CUDA_VERSION < 9000
    GMX_UNUSED_VALUE(activeMask);
    return __ballot(pred);
#else
    return __ballot_sync(activeMask, pred);
#endif
}

/*! \brief Compatibility wrapper around the CUDA __any()/__any_sync() instrinsic.  */
static __forceinline__ __device__
int gmx_any_sync(const unsigned int activeMask,
                 const int          pred)
{
#if GMX_CUDA_VERSION < 9000
    GMX_UNUSED_VALUE(activeMask);
    return __any(pred);
#else
    return __any_sync(activeMask, pred);
#endif
}

/*! \brief Compatibility wrapper around the CUDA __shfl_up()/__shfl_up_sync() instrinsic.  */
template <typename T>
static __forceinline__ __device__
T gmx_shfl_up_sync(const unsigned int activeMask,
                   const T            var,
                   unsigned int       offset,
                   int                width = warp_size)
{
#if GMX_CUDA_VERSION < 9000
    GMX_UNUSED_VALUE(activeMask);
    return __shfl_up(var, offset, width);
#else
    return __shfl_up_sync(activeMask, var, offset, width);
#endif
}

/*! \brief Compatibility wrapper around the CUDA __shfl_down()/__shfl_down_sync() instrinsic.  */
template <typename T>
static __forceinline__ __device__
T gmx_shfl_down_sync(const unsigned int activeMask,
                     const T            var,
                     unsigned int       offset,
                     int                width = warp_size)
{
#if GMX_CUDA_VERSION < 9000
    GMX_UNUSED_VALUE(activeMask);
    return __shfl_down(var, offset, width);
#else
    return __shfl_down_sync(activeMask, var, offset, width);
#endif
}

/* Tabulated/global memory routines below. */

/*! Load directly or using __ldg() when supported. */
template<typename T>
__device__ __forceinline__ T LDG(const T* ptr)
{
#if GMX_PTX_ARCH >= 350
    /* CC >=3.5 supports constant loads through texture or L1 */
    return __ldg(ptr);
#else
    /* Device does not have LDG support, fall back to direct load. */
    return *ptr;
#endif
}

/*! \brief Fetch the value by \p index from the texture object or reference.
 * Fetching from the object is the preferred behaviour on CC >= 3.0.
 *
 * \tparam[in] T        Raw data type
 * \param[in] texObj    Table texture object
 * \param[in] texRef    Table texture reference
 * \param[in] index     Non-negative element index
 * \returns             The value from the table at \p index
 */
template <typename T>
static __forceinline__ __device__
T fetchFromTexture(const cudaTextureObject_t texObj,
                   const struct texture<T, 1, cudaReadModeElementType> texRef,
                   int index)
{
    assert(index >= 0);
    assert(!c_disableCudaTextures);
    T result;
#if GMX_PTX_ARCH >= 300  // Preferring texture objects on any new arch
    GMX_UNUSED_VALUE(texRef);
    result = tex1Dfetch<T>(texObj, index);
#else
    GMX_UNUSED_VALUE(texObj);
    result = tex1Dfetch(texRef, index);
#endif
    return result;
}

/*! \brief Fetch the value by \p index from the parameter lookup table.
 *
 *  Depending on what is supported, it fetches parameters either
 *  using direct load, texture objects, or texture references.
 *
 * \tparam[in] T        Raw data type
 * \param[in] d_ptr     Device pointer to the raw table memory
 * \param[in] texObj    Table texture object
 * \param[in] texRef    Table texture reference
 * \param[in] index     Non-negative element index
 * \returns             The value from the table at \p index
 */
template <typename T>
static __forceinline__ __device__
T fetchFromParamLookupTable(const T                  *d_ptr,
                            const cudaTextureObject_t texObj,
                            const struct texture<T, 1, cudaReadModeElementType> texRef,
                            int index)
{
    assert(index >= 0);
    T result;
#if DISABLE_CUDA_TEXTURES
    GMX_UNUSED_VALUE(texObj);
    GMX_UNUSED_VALUE(texRef);
    result = LDG(d_ptr + index);
#else
    GMX_UNUSED_VALUE(d_ptr);
    result = fetchFromTexture<T>(texObj, texRef, index);
#endif
    return result;
}


#endif /* GMX_GPU_UTILS_CUDA_KERNEL_UTILS_CUH */
