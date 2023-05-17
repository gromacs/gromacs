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
#ifndef GMX_GPU_UTILS_GPUTRAITS_SYCL_H
#define GMX_GPU_UTILS_GPUTRAITS_SYCL_H

/*! \libinternal \file
 *  \brief Declares the SYCL type traits.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_gpu_utils
 */

#include <cstddef>

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/math/vectypes.h"

//! Type of device texture object. In SYCL, that would be \c sycl::image, but it's not used.
using DeviceTexture = void*;

//! \brief Single GPU call timing event, not used with SYCL
using CommandEvent = void*;

// TODO: Issue #3312
//! Convenience alias.
using Float4 = sycl::float4;
//! Convenience alias. Not using sycl::float3 due to alignment issues.
using Float3 = gmx::RVec;
//! Convenience alias for sycl::float2
using Float2 = sycl::float2;

/*! \internal \brief
 * GPU kernels scheduling description.
 * One typically only needs to set non-1 work sizes.
 *
 * \note This struct uses CUDA/OpenCL layout, with the first dimension being contiguous.
 *       It is different from the SYCL standard, where the last dimension is contiguous.
 *       The transpose is to be performed internally in ISyclKernelFunctor::launch.
 * \note \c sharedMemorySize is ignored in SYCL.
 */
struct KernelLaunchConfig
{
    //! Work groups (CUDA blocks) counts
    size_t gridSize[3] = { 1, 1, 1 };
    //! Per work group (CUDA block) thread counts
    size_t blockSize[3] = { 1, 1, 1 };
    //! Shared memory size in bytes
    size_t sharedMemorySize = 0;
};

/*! \brief Sets whether device code can use arrays that are embedded in structs.
 *
 * That is not technically true for SYCL, since DeviceBuffer holds not only the memory pointer,
 * but also the context.
 *
 * But our \c prepareGpuKernelArguments and \c launchGpuKernel functions deal
 * with that, so we can pass embedded buffers to them, which is what this
 * constant actually controls.
 */
#define c_canEmbedBuffers true

#endif
