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
    if (GMX_GPU != GMX_GPU_CUDA)
    {
        GMX_RELEASE_ASSERT(pinningPolicy == PinningPolicy::CannotBePinned,
                           "A suitable build of GROMACS (e.g. with CUDA) is required for a "
                           "HostAllocationPolicy to be set to a mode that produces pinning.");
    }
}

std::size_t HostAllocationPolicy::alignment()
{
    return (pinningPolicy_ == PinningPolicy::CanBePinned ?
            PageAlignedAllocationPolicy::alignment() :
            AlignedAllocationPolicy::alignment());
}

/*! \brief  Pin memory starting at p of size n
 *
 * For every pin unpin has to be called.
 * This is guranteed because for every p!=null&&pinningPolicy_==CanBePinned malloc/free call these.
 * For very std compliant container the allocator object can't be changed independent of the buffer
 * (for move it is propaged) and thus the allocator (and thus pinningPolicy_) can't change between
 * malloc and free.
 */
static void pin(void gmx_unused *p, size_t gmx_unused n) noexcept
{
#if GMX_GPU == GMX_GPU_CUDA
    //always pin (even for size 0) so that we can always unpin without any checks
    n = std::max<size_t>(1, n); //C++11 3.7.4.1 gurantees that every pointer is different thus at least 1 byte
    pinBuffer(p, n);
#endif
}

//! Unpin memory starting at p
static void unpin(void gmx_unused *p) noexcept
{
#if GMX_GPU == GMX_GPU_CUDA
    unpinBuffer(p);
#endif
}

void *HostAllocationPolicy::malloc(std::size_t bytes) const noexcept
{
    if (pinningPolicy_ == PinningPolicy::CanBePinned)
    {
        void *p = PageAlignedAllocationPolicy::malloc(bytes);
        if (p)
        {
            pin(p, bytes);
        }
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
    if (pinningPolicy_ == PinningPolicy::CanBePinned)
    {
        unpin(buffer);
        PageAlignedAllocationPolicy::free(buffer);
    }
    else
    {
        AlignedAllocationPolicy::free(buffer);
    }
}

} // namespace gmx
