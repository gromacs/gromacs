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
/*! \file
 * \brief Declares allocation policy classes and allocators that are
 * used to make library containers compatible with alignment
 * requirements of particular hardware, e.g. memory operations for
 * SIMD or accelerators.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_ALIGNEDALLOCATOR_H
#define GMX_UTILITY_ALIGNEDALLOCATOR_H

#include <cstddef>

#include "gromacs/utility/allocator.h"

namespace gmx
{

/*! \libinternal \brief Policy class for configuring gmx::Allocator, to manage
 * allocations of aligned memory for SIMD code.
 */
class AlignedAllocationPolicy
{
    public:
        /*! \brief Return the alignment size. */
        static std::size_t
        alignment();
        /*! \brief Allocate memory aligned to alignment() bytes.
         *
         *  \param bytes Amount of memory (bytes) to allocate. It is valid to ask for
         *               0 bytes, which will return a non-null pointer that is properly
         *               aligned and padded (but that you should not use).
         *
         * \return Valid pointer if the allocation worked, otherwise nullptr.
         *
         * The memory will always be aligned to 128 bytes, which is our
         * estimate of the longest cache lines on architectures currently in use.
         * It will also be padded by the same amount at the end of the
         * area, to help avoid false cache sharing.
         *
         *  \note Memory allocated with this routine must be released with
         *        gmx::AlignedAllocationPolicy::free(), and absolutely not the system free().
         */
        static void *
        malloc(std::size_t bytes);
        /*! \brief Free aligned memory
         *
         *  \param p  Memory pointer previously returned from malloc()
         *
         *  \note This routine should only be called with pointers obtained from
         *        gmx::AlignedAllocationPolicy::malloc(), and absolutely not any
         *        pointers obtained the system malloc().
         */
        static void
        free(void *p);
};

/*! \brief Aligned memory allocator.
 *
 *  \tparam T          Type of objects to allocate
 *
 * This convenience partial specialization can be used for the
 * optional allocator template parameter in standard library
 * containers, which is necessary e.g. to use SIMD aligned load and
 * store operations on data in those containers. The memory will
 * always be aligned according to the behavior of
 * AlignedAllocationPolicy.
 */
template <class T>
using AlignedAllocator = Allocator<T, AlignedAllocationPolicy>;


/*! \brief Return the memory page size on this system
 *
 * Implements the "construct on first use" idiom to avoid the static
 * initialization order fiasco where a possible static page-aligned
 * container would be initialized before the alignment variable was.
 *
 * Note that thread-safety is guaranteed by the C++11 language
 * standard. */
std::size_t pageSize();

/*! \libinternal \brief Policy class for configuring gmx::Allocator,
 * to manage allocations of page-aligned memory that can be locked for
 * asynchronous transfer to GPU devices.
 */
class PageAlignedAllocationPolicy
{
    public:
        /*! \brief Return the alignment size of memory pages on this system.
         *
         * Queries sysconf/WinAPI, otherwise guesses 4096. */
        static std::size_t
        alignment();
        /*! \brief Allocate memory aligned to alignment() bytes.
         *
         *  \param bytes Amount of memory (bytes) to allocate. It is valid to ask for
         *               0 bytes, which will return a non-null pointer that is properly
         *               aligned and padded (but that you should not use).
         *
         * \return Valid pointer if the allocation worked, otherwise nullptr.
         *
         *  \note Memory allocated with this routine must be released with
         *        gmx::PageAlignedAllocationPolicy::free(), and absolutely not the system free().
         */
        static void *
        malloc(std::size_t bytes);
        /*! \brief Free aligned memory
         *
         *  \param p  Memory pointer previously returned from malloc()
         *
         *  \note This routine should only be called with pointers obtained from
         *        gmx::PageAlignedAllocationPolicy::malloc(), and absolutely not any
         *        pointers obtained the system malloc().
         */
        static void
        free(void *p);
};

/*! \brief PageAligned memory allocator.
 *
 *  \tparam T          Type of objects to allocate
 *
 * This convenience partial specialization can be used for the
 * optional allocator template parameter in standard library
 * containers, which is necessary for locking memory pages for
 * asynchronous transfer between a GPU device and the host.  The
 * memory will always be aligned according to the behavior of
 * PageAlignedAllocationPolicy.
 */
template <class T>
using PageAlignedAllocator = Allocator<T, PageAlignedAllocationPolicy>;

}      // namespace gmx

#endif // GMX_UTILITY_ALIGNEDALLOCATOR_H
