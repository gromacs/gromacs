/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef DEVICEBUFFER_H
#define DEVICEBUFFER_H

// common base for devicebuffer.cuh/_ocl.h, only to be included by them

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h" // TODO: this is only for over_alloc_large

/*! \brief Reallocation device buffers
 *
 *  Reallocation of the memory pointed by d_ptr and copying of the data from
 *  the location pointed by h_src host-side pointer is done. Allocation is
 *  buffered and therefore freeing is only needed if the previously allocated
 *  space is not enough.
 *  The H2D copy is launched in command queue s and can be done synchronously or
 *  asynchronously (the default is the latter).
 *  If copy_event is not NULL, on return it will contain an event object
 *  identifying the H2D copy. The event can further be used to queue a wait
 *  for this operation or to query profiling information.
 *  OpenCL equivalent of cu_realloc_buffered.
 */
template <typename ValueType, typename Context>
void reallocateDeviceBuffer(DeviceBuffer<ValueType> *buffer,
                            size_t                   numValues,
                            int                     *currentNumValues,
                            int                     *currentMaxNumValues,
                            Context                  context)
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
        allocateDeviceBuffer(buffer, *currentMaxNumValues, context);
    }
    /* size could have changed without actual reallocation */
    *currentNumValues = numValues;
}

#endif // DEVICEBUFFER_H
