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
/*! \internal \file
 * \brief Implements gmx::HostAllocationPolicy for allocating memory
 * suitable for GPU transfers.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "hostallocator.h"

#include "config.h"

#include <cstdlib>

#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

void *
HostAllocationPolicy::malloc(std::size_t bytes)
{
#if GMX_GPU == GMX_GPU_CUDA
    void *buffer = nullptr;
    if (allocateForGpu_ == Impl::AllocateForGpu)
    {
        if (bytes != 0)
        {
            cudaError_t stat = cudaMallocHost(&buffer, bytes, cudaHostAllocDefault);
            if (stat != cudaSuccess)
            {
                buffer = nullptr;
            }
        }
    }
    else
    {
        buffer = AlignedAllocationPolicy::malloc(bytes);
    }
    return buffer;
#else
    GMX_UNUSED_VALUE(allocateForGpu_);
    // TODO if/when this is properly supported for OpenCL, it will
    // probably need to be aligned to a 4kb page, like CUDA does.
    return AlignedAllocationPolicy::malloc(bytes);
#endif
}

void
HostAllocationPolicy::free(void *buffer)
{
    if (buffer == nullptr)
    {
        return;
    }
#if GMX_GPU == GMX_GPU_CUDA
    if (allocateForGpu_ == Impl::AllocateForGpu)
    {
        cudaFreeHost(buffer);
    }
    else
    {
        AlignedAllocationPolicy::free(buffer);
    }
#else
    GMX_UNUSED_VALUE(allocateForGpu_);
    AlignedAllocationPolicy::free(buffer);
#endif
}

} // namespace gmx
