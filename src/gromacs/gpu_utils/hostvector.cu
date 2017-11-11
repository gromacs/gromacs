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
 * \brief Implements gmx::HostVector std-compatible container suitable
 * for GPU transfers on CUDA.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "hostvector.h"

#include <cstdint>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

template <typename T>
void HostVector<T>::usePinningMode(bool newMode)
{
    isPinned_ = newMode;
    pin();
}

//! Is \c ptr aligned on a boundary that is a multiple of \c bytes.
static inline bool isAligned(const void *ptr, size_t bytes)
{
    return (reinterpret_cast<intptr_t>(ptr) % bytes) == 0;
}

namespace
{

//! Ensure our error handling works as well as possible.
static void ensureNoPendingCudaError(const char *descriptiveErrorMessage)
{
    // Ensure there is no pending error that would otherwise affect
    // the behaviour of future error handling.
    cudaError_t stat = cudaGetLastError();
    // If we would find an error in a release build, we do not know
    // what is appropriate to do about it, so assert only for debug
    // builds.
    GMX_ASSERT(stat == cudaSuccess,
               formatString("%s An unhandled error from a previous CUDA operation was detected.", descriptiveErrorMessage).c_str());
}

}

template <typename T>
void HostVector<T>::pin()
{
    if (!usingPinningMode_ || empty() || isPinned_)
    {
        return;
    }
    const char *descriptiveErrorMessage = "Could not register the host memory for page locking for GPU transfers.";

    void       *pointer = v_.data();
    GMX_ASSERT(isAligned(pointer, PageAlignedAllocationPolicy::alignment()),
               formatString("%s Host memory needs to be page aligned.", descriptiveErrorMessage).c_str());

    ensureNoPendingCudaError(descriptiveErrorMessage);
    // Using the capacity of the vector as the size of the memory
    // registered with CUDA potentially uses more pages than needed,
    // but minimizes the number of times we need to un-register and
    // re-register.
    cudaError_t stat = cudaHostRegister(pointer, sizeof(value_type) * capacity(), cudaHostRegisterDefault);

    // These errors can only arise from a coding error somewhere.
    GMX_RELEASE_ASSERT(stat != cudaErrorInvalidValue && stat != cudaErrorNotSupported && stat != cudaErrorHostMemoryAlreadyRegistered,
                       formatString("%s %s: %s", cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());
    // Remaining errors come from correct use of the CUDA API, but it
    // could not do what we asked of it. That should only be
    // cudaErrorMemoryAllocation.
    if (stat != cudaSuccess)
    {
        // We always handle the error, but if it's a type we didn't
        // expect (e.g. because CUDA changes the set of errors it
        // returns) then we should get an assertion in Debug mode so
        // we know to fix our expectations.
        GMX_ASSERT(stat != cudaErrorMemoryAllocation,
                   formatString("%s %s: %s which was an unexpected error", cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());
        GMX_THROW(InternalError(formatString("%s %s: %s", descriptiveErrorMessage, cudaGetErrorName(stat), cudaGetErrorString(stat))));
    }
    isPinned_ = true;
}

template <typename T>
void HostVector<T>::unpin()
{
    if (!isPinned_)
    {
        return;
    }
    const char *descriptiveErrorMessage = "Could not unregister pinned host memory used for GPU transfers.";

    GMX_ASSERT(usingPinningMode_, formatString("%s HostVector should not be pinned if we're not in pinning mode.", descriptiveErrorMessage).c_str());;
    GMX_ASSERT(!empty(), formatString("%s HostVector should not be empty when pinned.", descriptiveErrorMessage).c_str());

    void       *pointer = v_.data();
    ensureNoPendingCudaError(descriptiveErrorMessage);
    cudaError_t stat = cudaHostUnregister(pointer);
    // These errors can only arise from a coding error somewhere.
    GMX_RELEASE_ASSERT(stat != cudaErrorInvalidValue && stat != cudaErrorHostMemoryNotRegistered,
                       formatString("%s %s: %s", cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());
    // No other errors should be issued, but if we get one, someone should know about it.
    if (stat != cudaSuccess)
    {
        // We always handle the error, but if it's a type we didn't
        // expect (e.g. because CUDA changes the set of errors it
        // returns) then we should get an assertion in Debug mode so
        // we know to fix our expectations.
        GMX_ASSERT(stat != cudaSuccess,
                   formatString("%s %s: %s which was an unexpected error", cudaGetErrorName(stat), cudaGetErrorString(stat)).c_str());
        GMX_THROW(InternalError(formatString("%s %s: %s", descriptiveErrorMessage, cudaGetErrorName(stat), cudaGetErrorString(stat))));
    }
    isPinned_ = false;
}

template <typename T>
void HostVector<T>::reserve(size_type newCapacity)
{
    bool shouldRepin = false;
    if (newCapacity > capacity())
    {
        unpin();
    }
    // Note, reserve copies data if it makes a new allocation.
    v_.reserve(newCapacity);
    if (shouldRepin)
    {
        pin();
    }
}

template <typename T>
void HostVector<T>::resize(size_type newSize)
{
    if (newSize > capacity())
    {
        reserve(newSize);
    }
#ifndef NDEBUG
    auto oldData = data();
#endif
    // Will not reallocate, so no need to manage pinning.
    v_.resize(newSize);
#ifndef NDEBUG
    GMX_ASSERT(data() == oldData, "The vector underlying HostVector should not have reallocated");
#endif
}

template <typename T>
void HostVector<T>::free()
{
    unpin();
    // Use the swap trick to give the allocation to a temporary that
    // will deallocate upon destruction.
    VectorType tmp;
    v_.swap(tmp);
    // There's nothing to pin.
}

//! Extern template instantiations. Extend the range as required.
/**@{*/
template class HostVector<int>;
template class HostVector<real>;
template class HostVector<RVec>;
/**@}*/

} // namespace gmx
