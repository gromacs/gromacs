/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

#if GMX_GPU == GMX_GPU_CUDA

//! \brief A device-side buffer of ValueTypes
template<typename ValueType>
using DeviceBuffer = ValueType*;

#elif GMX_GPU == GMX_GPU_OPENCL

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

#else

//! \brief A device-side buffer of ValueTypes
template<typename ValueType>
using DeviceBuffer = void*;

#endif


#endif // GMX_GPU_UTILS_DEVICEBUFFER_DATATYPE_H
