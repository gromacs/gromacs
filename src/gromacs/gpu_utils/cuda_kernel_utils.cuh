/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019, by the GROMACS development team, led by
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

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"

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

/*! \brief Fetch the value by \p index from the texture object.
 *
 * \tparam[in] T        Raw data type
 * \param[in] texObj    Table texture object
 * \param[in] index     Non-negative element index
 * \returns             The value from the table at \p index
 */
template<typename T>
static __forceinline__ __device__ T fetchFromTexture(const cudaTextureObject_t texObj, int index)
{
    assert(index >= 0);
    assert(!c_disableCudaTextures);
    return tex1Dfetch<T>(texObj, index);
}

/*! \brief Fetch the value by \p index from the parameter lookup table.
 *
 *  Depending on what is supported, it fetches parameters either
 *  using direct load or texture objects.
 *
 * \tparam[in] T        Raw data type
 * \param[in] d_ptr     Device pointer to the raw table memory
 * \param[in] texObj    Table texture object
 * \param[in] index     Non-negative element index
 * \returns             The value from the table at \p index
 */
template<typename T>
static __forceinline__ __device__ T fetchFromParamLookupTable(const T*                  d_ptr,
                                                              const cudaTextureObject_t texObj,
                                                              int                       index)
{
    assert(index >= 0);
    T result;
#if DISABLE_CUDA_TEXTURES
    GMX_UNUSED_VALUE(texObj);
    result = LDG(d_ptr + index);
#else
    GMX_UNUSED_VALUE(d_ptr);
    result = fetchFromTexture<T>(texObj, index);
#endif
    return result;
}


#endif /* GMX_GPU_UTILS_CUDA_KERNEL_UTILS_CUH */
