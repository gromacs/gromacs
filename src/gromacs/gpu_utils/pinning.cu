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
 * \brief Implements functions for pinning memory to be suitable for
 * efficient GPU transfers on CUDA.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "pinning.h"

#include <cstddef>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

//! Is \c ptr aligned on a boundary that is a multiple of \c bytes.
gmx_unused static inline bool isAligned(const void *ptr, size_t bytes)
{
    return (reinterpret_cast<intptr_t>(ptr) % bytes) == 0;
}

void pinBuffer(void *pointer, std::size_t numBytes) noexcept
{
    const char *errorMessage = "Could not register the host memory for page locking for GPU transfers.";

    GMX_ASSERT(isAligned(pointer, PageAlignedAllocationPolicy::alignment()),
               formatString("%s Host memory needs to be page aligned.", errorMessage).c_str());

    ensureNoPendingCudaError(errorMessage);
    cudaError_t stat = cudaHostRegister(pointer, numBytes, cudaHostRegisterDefault);

    // These errors can only arise from a coding error somewhere.
    GMX_RELEASE_ASSERT(stat != cudaErrorInvalidValue &&
                       stat != cudaErrorNotSupported &&
                       stat != cudaErrorHostMemoryAlreadyRegistered,
                       formatString("%s %s: %s", errorMessage,
                                    cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());

    // We always handle the error, but if it's a type we didn't expect
    // (e.g. because CUDA changes the set of errors it returns) then
    // we should get a descriptive assertion in Debug mode so we know
    // to fix our expectations.
    GMX_ASSERT(stat != cudaErrorMemoryAllocation,
               formatString("%s %s: %s which was an unexpected error", errorMessage,
                            cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());

    // It might be preferable to throw InternalError here, because the
    // failing condition can only happen when GROMACS is used with a
    // CUDA API that can return some other error code. But we can't
    // engineer GROMACS to be forward-compatible with future CUDA
    // versions, so if this proves to be a problem in practice, then
    // GROMACS must be patched, or a supported CUDA version used.
    GMX_RELEASE_ASSERT(stat == cudaSuccess,
                       formatString("%s %s: %s", errorMessage,
                                    cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());
}

void unpinBuffer(void *pointer) noexcept
{
    const char *errorMessage = "Could not unregister pinned host memory used for GPU transfers.";

    GMX_ASSERT(pointer != nullptr,
               formatString("%s pointer should not be nullptr when pinned.", errorMessage).c_str());

    ensureNoPendingCudaError(errorMessage);
    cudaError_t stat = cudaHostUnregister(pointer);
    // These errors can only arise from a coding error somewhere.
    GMX_RELEASE_ASSERT(stat != cudaErrorInvalidValue && stat != cudaErrorHostMemoryNotRegistered,
                       formatString("%s %s: %s", errorMessage,
                                    cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());
    // If there's an error whose type we didn't expect (e.g. because a
    // future CUDA changes the set of errors it returns) then we
    // should assert, because our code is wrong.
    //
    // The approach differs from that in pin() because we might
    // unpin() from a destructor, in which case any attempt to throw
    // an uncaught exception would anyway terminate the program. A
    // release assertion is a better behaviour than that.
    GMX_RELEASE_ASSERT(stat == cudaSuccess,
                       formatString("%s %s: %s which was an unexpected error", errorMessage,
                                    cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());
}

} // namespace gmx
