/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
/*! \internal \file
 * \brief Implements gmx::HostAllocationPolicy for allocating memory
 * suitable for e.g. GPU transfers on CUDA.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "hostallocator.h"

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/pmalloc.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

HostAllocationPolicy::HostAllocationPolicy() = default;

HostAllocationPolicy::HostAllocationPolicy(const bool propagateDuringContainerCopyConstruction) :
    propagateDuringContainerCopyConstruction_(propagateDuringContainerCopyConstruction)
{
}

HostAllocationPolicy::HostAllocationPolicy(const DeviceContext& deviceContext,
                                           PinningPolicy        pinningPolicy,
                                           const bool propagateDuringContainerCopyConstruction) :
    context_(&deviceContext),
    pinningPolicy_(pinningPolicy),
    propagateDuringContainerCopyConstruction_(propagateDuringContainerCopyConstruction)
{
}

std::size_t HostAllocationPolicy::alignment() const noexcept
{
    return (pinningPolicy_ == PinningPolicy::PinnedIfSupported ? PageAlignedAllocationPolicy::alignment()
                                                               : AlignedAllocationPolicy::alignment());
}

void* HostAllocationPolicy::malloc(std::size_t bytes) const noexcept
{
    if (pinningPolicy_ == PinningPolicy::PinnedIfSupported)
    {
        void* p;
        pmalloc(&p, bytes, &context());
        return p;
    }
    else
    {
        return AlignedAllocationPolicy::malloc(bytes);
    }
}

void HostAllocationPolicy::free(void* buffer) const noexcept
{
    if (buffer == nullptr)
    {
        // Nothing to do
        return;
    }
    if (pinningPolicy_ == PinningPolicy::PinnedIfSupported)
    {
        pfree(buffer, &context());
    }
    else
    {
        AlignedAllocationPolicy::free(buffer);
    }
}

bool HostAllocationPolicy::operator==(const HostAllocationPolicy& b) const
{
    // Currently GROMACS only uses one context per rank, so the
    // context is always equal when valid, but it can be null.
    return (this->context_ == b.context_) && (this->pinningPolicy() == b.pinningPolicy());
}

const DeviceContext& HostAllocationPolicy::context() const
{
    GMX_RELEASE_ASSERT(
            context_ != nullptr,
            "Cannot return a DeviceContext from a HostAllocationPolicy that was not made with one");
    return *context_;
}

HostAllocationPolicy makeHostAllocationPolicy(const bool                 pinBuffers,
                                              const DeviceStreamManager* deviceStreamManager,
                                              const bool propagateDuringContainerCopyConstruction)
{
    GMX_RELEASE_ASSERT(!pinBuffers || deviceStreamManager != nullptr,
                       "Must have a valid deviceStreamManager to allocate to use GPU transfers");
    return pinBuffers && deviceStreamManager
                   ? HostAllocationPolicy{ deviceStreamManager->context(),
                                           PinningPolicy::PinnedIfSupported,
                                           propagateDuringContainerCopyConstruction }
                   : HostAllocationPolicy{};
}

} // namespace gmx
