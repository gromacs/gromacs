/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
#ifndef GMX_GPU_UTILS_DEVICEBUFFER_DATATYPE_H
#define GMX_GPU_UTILS_DEVICEBUFFER_DATATYPE_H

/*! \libinternal \file
 *  \brief Declares the DeviceBuffer data type.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *
 *  \inlibraryapi
 */

#include "config.h"

#include <memory>

#include "gromacs/math/vectypes.h"

#if GMX_GPU_CUDA || GMX_GPU_HIP

//! \brief A device-side buffer of ValueTypes
template<typename ValueType>
using DeviceBuffer = ValueType*;

#elif GMX_GPU_OPENCL

#    include "gromacs/gpu_utils/gputraits_ocl.h"

/*! \libinternal \brief
 * A minimal cl_mem wrapper that remembers its allocation type.
 * The only point is making template type deduction possible.
 */
template<typename ValueType>
class TypedClMemory
{
private:
    /*! \brief Underlying data
     *
     *  \todo Make it nullptr when there are no snew()'s around
     */
    cl_mem data_;

public:
    //! \brief Default constructor
    TypedClMemory() : data_(nullptr) {}

    //! \brief Needed for cross-platform compilation
    TypedClMemory(std::nullptr_t nullPtr) : data_(nullPtr) {}

    //! \brief An assignment operator - the purpose is to make allocation/zeroing work
    TypedClMemory& operator=(cl_mem data)
    {
        data_ = data;
        return *this;
    }
    //! \brief Returns underlying cl_mem transparently
    operator cl_mem() { return data_; }
};

//! \libinternal \brief A device-side buffer of ValueTypes
template<typename ValueType>
using DeviceBuffer = TypedClMemory<ValueType>;

#elif GMX_GPU_SYCL

/*! \libinternal \brief
 * A minimal wrapper around SYCL USM memory and assistant data structures.
 */
template<typename ValueType>
struct DeviceBuffer
{
    class SyclBufferWrapper;
    std::unique_ptr<SyclBufferWrapper> buffer_;

    DeviceBuffer();
    DeviceBuffer(std::nullptr_t nullPtr);
    ~DeviceBuffer();
    DeviceBuffer(DeviceBuffer<ValueType> const& src);
    DeviceBuffer(DeviceBuffer<ValueType>&& src) noexcept;
    DeviceBuffer& operator=(DeviceBuffer<ValueType> const& src);
    DeviceBuffer& operator=(DeviceBuffer<ValueType>&& src) noexcept;

    //! Helper function to get the size in bytes of a single element
    static constexpr size_t elementSize() { return sizeof(ValueType); }

    //! Both explicit and implicit casts to void* are used in MPI+CUDA code, this stub is necessary for compilation.
    operator void*() const { throw; }

    //! Get underlying device const pointer
    const ValueType* get_pointer() const;

    //! Get underlying device pointer
    ValueType* get_pointer();

    //! Allow implicit conversion to bool to check buffer status for compatibility with other implementations.
    operator bool() const { return buffer_.get() != nullptr; }

    //! An assignment operator to allow resetting buffer by assigning nullptr to it. Necessary for compilation.
    DeviceBuffer& operator=(std::nullptr_t nullPtr);
};

// Must explicitly instantiate for some types.
extern template struct DeviceBuffer<gmx::RVec>;
extern template struct DeviceBuffer<uint64_t>;

#else

//! \brief A device-side buffer of ValueTypes
template<typename ValueType>
using DeviceBuffer = void*;

#endif

#endif // GMX_GPU_UTILS_DEVICEBUFFER_DATATYPE_H
