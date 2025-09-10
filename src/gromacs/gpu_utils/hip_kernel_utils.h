/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
#ifndef GMX_GPU_UTILS_HIP_KERNEL_UTILS_H
#define GMX_GPU_UTILS_HIP_KERNEL_UTILS_H

/*! \file
 *  \brief HIP device util functions (callable from GPU kernels).
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 */

#include <hip/hip_runtime.h>

#include "gromacs/utility/basedefinitions.h"

/*! Load directly or using __ldg() when supported. */
template<typename T>
__device__ __forceinline__ T LDG(const T* ptr)
{
    return __ldg(ptr);
}

/*! \brief Fetch the value by \p index from the parameter lookup table.
 *
 *  Depending on what is supported, it fetches parameters either
 *  using direct load or texture objects.
 *
 * \tparam T            Raw data type
 * \param[in] d_ptr     Device pointer to the raw table memory
 * \param[in] texObj    Table texture object
 * \param[in] index     Non-negative element index
 * \returns             The value from the table at \p index
 */
template<typename T>
static __forceinline__ __device__ T fetchFromParamLookupTable(const T*                 d_ptr,
                                                              const hipTextureObject_t texObj,
                                                              int                      index)
{
    assert(index >= 0);
    T result;
    GMX_UNUSED_VALUE(texObj);
    result = LDG(d_ptr + index);
    return result;
}

#define LAUNCH_BOUNDS_EXACT(WORK_GROUP_SIZE, WAVES_PER_EU)                        \
    __attribute__((amdgpu_flat_work_group_size(WORK_GROUP_SIZE, WORK_GROUP_SIZE), \
                   amdgpu_waves_per_eu(WAVES_PER_EU, WAVES_PER_EU)))

#define LAUNCH_BOUNDS_EXACT_SINGLE(WORK_GROUP_SIZE) \
    __attribute__((amdgpu_flat_work_group_size(WORK_GROUP_SIZE, WORK_GROUP_SIZE)))

#define GMX_HIP_MAX_BLOCKS_PER_MP 16
#define GMX_HIP_MAX_THREADS_PER_MP 1024

#endif /* GMX_GPU_UTILS_HIP_KERNEL_UTILS_H */
