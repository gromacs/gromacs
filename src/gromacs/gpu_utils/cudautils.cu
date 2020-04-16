/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015,2016,2017 The GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "cudautils.cuh"

#include <cassert>
#include <cstdlib>

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/utility/gmxassert.h"

/*** Generic CUDA data operation wrappers ***/

// TODO: template on transferKind to avoid runtime conditionals
int cu_copy_D2H(void* h_dest, void* d_src, size_t bytes, GpuApiCallBehavior transferKind, cudaStream_t s = nullptr)
{
    cudaError_t stat;

    if (h_dest == nullptr || d_src == nullptr || bytes == 0)
    {
        return -1;
    }

    switch (transferKind)
    {
        case GpuApiCallBehavior::Async:
            GMX_ASSERT(isHostMemoryPinned(h_dest), "Destination buffer was not pinned for CUDA");
            stat = cudaMemcpyAsync(h_dest, d_src, bytes, cudaMemcpyDeviceToHost, s);
            CU_RET_ERR(stat, "DtoH cudaMemcpyAsync failed");
            break;

        case GpuApiCallBehavior::Sync:
            stat = cudaMemcpy(h_dest, d_src, bytes, cudaMemcpyDeviceToHost);
            CU_RET_ERR(stat, "DtoH cudaMemcpy failed");
            break;

        default: throw;
    }

    return 0;
}

int cu_copy_D2H_sync(void* h_dest, void* d_src, size_t bytes)
{
    return cu_copy_D2H(h_dest, d_src, bytes, GpuApiCallBehavior::Sync);
}

/*!
 *  The copy is launched in stream s or if not specified, in stream 0.
 */
int cu_copy_D2H_async(void* h_dest, void* d_src, size_t bytes, cudaStream_t s = nullptr)
{
    return cu_copy_D2H(h_dest, d_src, bytes, GpuApiCallBehavior::Async, s);
}

// TODO: template on transferKind to avoid runtime conditionals
int cu_copy_H2D(void* d_dest, const void* h_src, size_t bytes, GpuApiCallBehavior transferKind, cudaStream_t s = nullptr)
{
    cudaError_t stat;

    if (d_dest == nullptr || h_src == nullptr || bytes == 0)
    {
        return -1;
    }

    switch (transferKind)
    {
        case GpuApiCallBehavior::Async:
            GMX_ASSERT(isHostMemoryPinned(h_src), "Source buffer was not pinned for CUDA");
            stat = cudaMemcpyAsync(d_dest, h_src, bytes, cudaMemcpyHostToDevice, s);
            CU_RET_ERR(stat, "HtoD cudaMemcpyAsync failed");
            break;

        case GpuApiCallBehavior::Sync:
            stat = cudaMemcpy(d_dest, h_src, bytes, cudaMemcpyHostToDevice);
            CU_RET_ERR(stat, "HtoD cudaMemcpy failed");
            break;

        default: throw;
    }

    return 0;
}

int cu_copy_H2D_sync(void* d_dest, const void* h_src, size_t bytes)
{
    return cu_copy_H2D(d_dest, h_src, bytes, GpuApiCallBehavior::Sync);
}

/*!
 *  The copy is launched in stream s or if not specified, in stream 0.
 */
int cu_copy_H2D_async(void* d_dest, const void* h_src, size_t bytes, cudaStream_t s = nullptr)
{
    return cu_copy_H2D(d_dest, h_src, bytes, GpuApiCallBehavior::Async, s);
}
