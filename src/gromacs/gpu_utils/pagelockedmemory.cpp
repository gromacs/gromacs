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
 * \brief Defines PageLockedMemory for handling locking previously
 * allocated page-aligned memory to a physical memory page for use in
 * asynchronous GPU transfers.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gmxpre.h"

#include "pagelockedmemory.h"

#include "config.h"

#include <cinttypes>

#include <utility>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{

// Doxygen is too dumb for extern templates
#if !defined DOXYGEN

template <typename T>
PageLockedMemory::PageLockedMemory(ConstArrayRef<T> arrayRef) : memory_(nullptr)
{
    if (arrayRef.empty())
    {
        return;
    }
    memory_ = arrayRef.data();
#if GMX_GPU == GMX_GPU_CUDA

    GMX_ASSERT((reinterpret_cast<intptr_t>(memory_) % PageAlignedAllocationPolicy::alignment()) == 0,
               "Host memory is not page-aligned for page locking for GPU transfers");
    // TODO In Debug mode, check for pointer having been already registered?
    cudaError_t stat = cudaHostRegister(const_cast<void *>(memory_), arrayRef.size() * sizeof(T), cudaHostRegisterDefault);
    if (stat != cudaSuccess)
    {
        GMX_THROW(InternalError("Could not register the host memory for page locking for GPU transfers"));
    }
#elif GMX_GPU == GMX_GPU_OPENCL
    // There's nothing to do, page locking is not supported.
#else
    GMX_THROW(NotImplementedError("lockHostMemoryToPage is not implemented for this GROMACS build, so cannot lock allocated memory"));
#endif
}

// Explicit instantiations. Add as more are required.
template PageLockedMemory::PageLockedMemory(ConstArrayRef<real>);
template PageLockedMemory::PageLockedMemory(ConstArrayRef<RVec>);

#endif // DOXYGEN

PageLockedMemory::~PageLockedMemory()
{
#if GMX_GPU == GMX_GPU_CUDA
    // Swallow any error - there's nothing we can or should do about
    // being unable to un-register successfully, and any error code
    // might be from a previous asynchronous operation, too.
    cudaHostUnregister(const_cast<void *>(memory_));
#else
    // No need to do anything or give any error.
#endif
}

void PageLockedMemory::swap(PageLockedMemory &other)
{
    std::swap(memory_, other.memory_);
}

const void *PageLockedMemory::memory() const
{
    return memory_;
};

} // namespace gmx
