/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
#ifndef GMX_GPU_UTILS_DEVICEBUFFER_H
#define GMX_GPU_UTILS_DEVICEBUFFER_H

/*! \libinternal \file
 *  \brief Implements the logic for handling of DeviceBuffer types in OpenCL/CUDA.
 *  Can only be included on GPU build paths.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 *
 *  \inlibraryapi
 */

#include "config.h"

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h" // TODO: this is only for over_alloc_large

#if GMX_GPU == GMX_GPU_CUDA
#    include "gromacs/gpu_utils/devicebuffer.cuh"
#elif GMX_GPU == GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/devicebuffer_ocl.h"
#else
#    error "devicebuffer.h included on non-GPU build!"
#endif

/*! \brief
 *  Reallocates the device-side buffer.
 *
 *  Reallocates the device-side memory pointed by \p buffer.
 *  Allocation is buffered and therefore freeing is only needed
 *  if the previously allocated space is not enough.
 *  \p currentNumValues and \p currentMaxNumValues are updated.
 *  TODO: \p currentNumValues, \p currentMaxNumValues, \p deviceContext
 *  should all be encapsulated in a host-side class together with the buffer.
 *
 *  \tparam        ValueType            Raw value type of the \p buffer.
 *  \param[in,out] buffer               Pointer to the device-side buffer
 *  \param[in]     numValues            Number of values to accommodate.
 *  \param[in,out] currentNumValues     The pointer to the buffer's number of values.
 *  \param[in,out] currentMaxNumValues  The pointer to the buffer's capacity.
 *  \param[in]     deviceContext        The buffer's device context.
 */
template<typename ValueType>
void reallocateDeviceBuffer(DeviceBuffer<ValueType>* buffer,
                            size_t                   numValues,
                            int*                     currentNumValues,
                            int*                     currentMaxNumValues,
                            const DeviceContext&     deviceContext)
{
    GMX_ASSERT(buffer, "needs a buffer pointer");
    GMX_ASSERT(currentNumValues, "needs a size pointer");
    GMX_ASSERT(currentMaxNumValues, "needs a capacity pointer");

    /* reallocate only if the data does not fit */
    if (static_cast<int>(numValues) > *currentMaxNumValues)
    {
        if (*currentMaxNumValues >= 0)
        {
            freeDeviceBuffer(buffer);
        }

        *currentMaxNumValues = over_alloc_large(numValues);
        allocateDeviceBuffer(buffer, *currentMaxNumValues, deviceContext);
    }
    /* size could have changed without actual reallocation */
    *currentNumValues = numValues;
}

#endif
