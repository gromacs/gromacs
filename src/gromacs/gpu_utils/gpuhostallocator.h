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
/*! \libinternal \file
 * \brief
 * Declares gmx::GpuHostAllocator that is used to make standard library
 * containers that use page-locked memory.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_GPUHOSTALLOCATOR_H
#define GMX_UTILITY_GPUHOSTALLOCATOR_H

#include <cstddef>

#include "gromacs/utility/allocator.h"

namespace gmx
{

/*! \libinternal
 * \brief Policy class for configuring gmx::Allocator, to manage
 * allocations of memory that is locked to pages on the host side and
 * thus able to be transferred automatically by the CUDA runtime.
 */
class GpuHostAllocationPolicy
{
    public:
        /*! \brief Allocate GPU memory
         *
         *  \param bytes Amount of memory (bytes) to allocate. It is valid to ask for
         *               0 bytes, which will return a non-null pointer that is properly
         *               aligned in page-locked memory (but that you should not use). TODO check this.
         *
         * \return Valid pointer if the allocation worked, otherwise nullptr.
         *
         * The memory will always be TODO
         *
         *  \note Memory allocated with this routine must be released with
         *        gmx::GpuHostAllocationPolicy::free(), and absolutely not the system free().
         */
        static void *
        malloc(std::size_t bytes);
        /*! \brief Free GPU memory
         *
         *  \param p  Memory pointer previously returned from gmx::GpuHostAllocationPolicy::malloc()
         *
         *  \note This routine should only be called with pointers obtained from
         *        gmx:GpuHostAllocationPolicy::malloc(), and absolutely not any pointers obtained
         *        the system malloc().
         */
        static void
        free(void *p);
};

/*! \libinternal \brief GPU memory allocator for memory that is page-locked on the host.
 *
 *  \tparam T          Type of objects to allocate
 *
 * This convenience partial specialization can be used for the
 * optional allocator template parameter in standard library
 * containers, which is necessary e.g. to use transfer the data in
 * those containers between host and device. The memory will always be
 * aligned according to the behavior of GpuHostAllocationPolicy.
 */
template <class T>
using GpuHostAllocator = Allocator<T, GpuHostAllocationPolicy>;

}      // namespace gmx

#endif
