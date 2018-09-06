/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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
 * suitable for e.g. GPU transfers on CUDA.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "hostallocator.h"

#include "config.h"

#include <cstddef>

#include <memory>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/pinning.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

HostAllocationPolicy::HostAllocationPolicy(PinningPolicy pinningPolicy)
    : pinningPolicy_(pinningPolicy)
{
}

std::size_t HostAllocationPolicy::alignment()
{
    return (pinningPolicy_ == PinningPolicy::PinnedIfSupported ?
            PageAlignedAllocationPolicy::alignment() :
            AlignedAllocationPolicy::alignment());
}

void *HostAllocationPolicy::malloc(std::size_t bytes) const noexcept
{
    if (pinningPolicy_ == PinningPolicy::PinnedIfSupported)
    {
        void *p = PageAlignedAllocationPolicy::malloc(bytes);
        pin(p, bytes);
        return p;
    }
    else
    {
        return AlignedAllocationPolicy::malloc(bytes);
    }
}

void HostAllocationPolicy::free(void *buffer) const noexcept
{
    if (buffer == nullptr)
    {
        // Nothing to do
        return;
    }
    unpin(buffer);
    if (pinningPolicy_ == PinningPolicy::PinnedIfSupported)
    {
        PageAlignedAllocationPolicy::free(buffer);
    }
    else
    {
        AlignedAllocationPolicy::free(buffer);
    }
}

void HostAllocationPolicy::pin(void gmx_unused *p, size_t gmx_unused n) const noexcept
{
#if GMX_GPU == GMX_GPU_CUDA
    // I believe this if statement isn't required for calls from malloc. But it
    // is required for the unit tests. Which might not make much sense to test
    // cases which can't actually happen for an allocator.
    if (p == nullptr || n == 0 || isHostMemoryPinned(p))
    {
        return;
    }
    pinBuffer(p, n);
#endif
}

void HostAllocationPolicy::unpin(void gmx_unused *p) const noexcept
{
#if GMX_GPU == GMX_GPU_CUDA
    if (isHostMemoryPinned(p))
    {
        unpinBuffer(p);
    }
#endif
}

} // namespace gmx
