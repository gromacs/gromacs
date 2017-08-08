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
   \brief
   Declares gmx::DetailedAlignedAllocator
   The class is derived from Erik Lindahls gmx::DetailedAlignedAllocator
   with the difference that alignment and padding are not hard-coded
   but adjustable.

   \todo
   This should be merged with gromacs/utility/alignedallocator.{h,cpp}
   if needed long-term.

   DetailedAlignedAllocator is used in containers to allocate memory
   compatible with SIMD aligned load/store requirements.

   \author Erik Lindahl <erik.lindahl@gmail.com>
   \author R. Thomas Ullmann <tullman@gwdg.de>

   \ingroup module_utility
 */
#include "gmxpre.h"

#include "detailed_aligned_allocator.h"

#include "config.h"

#include <cstdlib>

#if HAVE_MM_MALLOC_H
#    include <mm_malloc.h>
#elif HAVE_MALLOC_H
#    include <malloc.h>
#elif HAVE_XMMINTRIN_H
#    include <xmmintrin.h>
#endif

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

namespace internal
{

void *
alignedMalloc(std::size_t bytes, std::size_t alignment, std::size_t padding)
{
    void   *    p;

    // Pad memory at the end with another alignment bytes to avoid false sharing
    std::size_t pbytes = bytes + padding;

#if HAVE__MM_MALLOC
    p = _mm_malloc(pbytes, alignment);
#elif HAVE_POSIX_MEMALIGN
    if (posix_memalign(&p, alignment, pbytes) != 0)
    {
        p = nullptr;
    }
#elif HAVE_MEMALIGN
    p = memalign(alignment, pbytes);
#elif HAVE__ALIGNED_MALLOC
    p = _aligned_malloc(pbytes, alignment);
#else
    p = gmx::internal::alignedMallocGeneric(pbytes, alignment);
#endif

    return p;
}

}   // end namespace internal

}   // end namespace gmx
