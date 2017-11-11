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

#include <cstdlib>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

template <typename T>
void HostVector<T>::useLockingMode(bool newMode)
{
    isLocked_ = newMode;
    lock();
}

//! Is \c ptr aligned on a boundary that is a multiple of \c bytes.
static inline bool isAligned(const void *ptr, size_t bytes)
{
    return (reinterpret_cast<intptr_t>(ptr) % bytes) == 0;
}

template <typename T>
void HostVector<T>::lock()
{
    if (!usingLockingMode_ || empty() || isLocked_)
    {
        return;
    }

    void       *pointer = v_.data();
    GMX_ASSERT(isAligned(pointer, PageAlignedAllocationPolicy::alignment()), "Host memory needs to be page aligned");
    cudaError_t stat;
    // Using the capacity of the vector as the size of the memory
    // registered with CUDA potentially uses more pages than needed,
    // but minimizes the number of times we need to un-register and
    // re-register.
    stat = cudaHostRegister(pointer, sizeof(value_type) * capacity(), cudaHostRegisterDefault);
    if (stat != cudaSuccess)
    {
        GMX_THROW(InternalError(formatString("Could not register the host memory for page locking for GPU transfers : %s", cudaGetErrorString(stat))));
    }
    isLocked_ = true;
}

template <typename T>
void HostVector<T>::unlock()
{
    if (!isLocked_)
    {
        return;
    }

    GMX_ASSERT(usingLockingMode_, "HostVector should not be locked if we're not in locking mode");
    GMX_ASSERT(!empty(), "HostVector should not be empty when locked");

    void       *pointer = v_.data();
    cudaError_t stat;
    stat = cudaHostUnregister(pointer);
    if (stat != cudaSuccess)
    {
        GMX_THROW(InternalError(formatString("Could not unregister the page-locked host memory used for GPU transfers : %s", cudaGetErrorString(stat))));
    }
    isLocked_ = false;
}

template <typename T>
void HostVector<T>::reserve(size_type newCapacity)
{
    bool shouldRelock = false;
    if (newCapacity > capacity())
    {
        unlock();
    }
    // Note, reserve copies data if it makes a new allocation.
    v_.reserve(newCapacity);
    if (shouldRelock)
    {
        lock();
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
    // Will not reallocate, so no need to manage locking.
    v_.resize(newSize);
#ifndef NDEBUG
    GMX_ASSERT(data() == oldData, "The vector underlying HostVector should not have reallocated");
#endif
}

template <typename T>
void HostVector<T>::free()
{
    unlock();
    // Use the swap trick to give the allocation to a temporary that
    // will deallocate upon destruction.
    VectorType tmp;
    v_.swap(tmp);
    // There's nothing to lock.
}

//! Extern template instantiations.
/**@{*/
template class HostVector<int>;
template class HostVector<real>;
template class HostVector<RVec>;
/**@}*/

} // namespace gmx
