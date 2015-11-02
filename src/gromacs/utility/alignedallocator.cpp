/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * \brief
 * Implements AlignedAllocator.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "alignedallocator.h"

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

/*! \brief Allocate aligned memory in a fully portable way
 *
 *  \param bytes  Amount of memory (bytes) to allocate. The routine will return
 *                nullptr if the allocation fails. However, note that asking for
 *                zero bytes will return a pointer that is non-null and properly
 *                aligned (but obviously you cannot use it, since you promised
 *                not to access data beyond the 0 bytes you asked for).
 *
 *  \param alignment  Alignment specification in bytes, must be a power of 2.
 *
 * \return Nonzero pointer if the allocation worked, otherwise nullptr.
 *  This routine should only be called from alignedMalloc(), which also does
 *  the checking for valid values. This particular function is used for platforms
 *  where we have no control of the alignment of memory returned by the system.
 *  Instead, we increase the amount of memory requested internally such that we
 *  both can create a pointer inside this memory that fulfills the memory
 *  alignment requested, and that we have room to store the original pointer
 *  just before this area.
 *
 *  \note This is an internal routine that should only be called from
 *        gmx::alignedMalloc(). Just like system-provided routines, it provides
 *        memory that is aligned - but not padded.
 */
static void *
alignedMallocGeneric(std::size_t bytes, std::size_t alignment)
{
    // The amount of extra memory (beyound what the user asked for) we need is:
    // - sizeof(void *), to store the original pointer
    // - alignment, to make sure we have an aligned pointer in the area
    void * pMalloc = malloc(bytes + sizeof(void *) + alignment);

    if (pMalloc == nullptr)
    {
        return nullptr;
    }

    // Convert pMalloc to size_t (so we work with raw bytes), add the space we
    // need to save the original pointer, and (alignment-1) bytes, and then mask
    // out the lowest bits.
    std::size_t mask     = ~static_cast<std::size_t>(alignment-1);
    void      * pAligned = reinterpret_cast<void *>((reinterpret_cast<std::size_t>(pMalloc) + sizeof(void *) + alignment - 1) & mask);

    // Store original pointer. Since we allocated at least sizeof(void *) extra
    // space this is always a valid memory location.
    reinterpret_cast<void **>(pAligned)[-1] = pMalloc;

    return pAligned;
}


/*! \brief Free aligned memory
 *
 *  \param p  Memory pointer previously returned from
 *            gmx::internal::alignedFreePortable().
 *
 *  Since this routine relies on the original pointer being stored just before
 *  the memory area p points to, bad things will happen if you call this routine
 *  with a pointer obtained any other way, or if you call the system free()
 *  with a pointer obtained from std::alignedMalloc().
 *
 * \note  This is an internal routine that should only be called from
 *        gmx::alignedFree().
 */
static void
alignedFreeGeneric(void *p)
{
    if (p)
    {
        // Pick up the pointer stored just below p, and use that to call free()
        free( reinterpret_cast<void **>(p)[-1] );
    }
}



void *
alignedMalloc(std::size_t bytes)
{
    // For now we always use 128-byte alignment:
    // 1) IBM Power already has cache lines of 128-bytes, and needs it.
    // 2) x86 has 64 byte cache lines, but since a future AVX-1024 (rumored?)
    //    will need 1024/8=128 byte SIMD alignment, it is safer to use that
    //    already now.
    // 3) The old Pentium4 used 256-byte cache prefetching (but 64-byte lines).
    //    However, it's not worth worrying about performance for P4...
    // 4) ARM & Sparc have 64 byte lines, but will be just fine with
    //    128-byte alignment (nobody knows what the future brings)
    //
    // So, for now we're semi-lazy and just align to 128 bytes!
    //
    // TODO LINCS code is copying this assumption independently (for now)
    std::size_t alignment = 128;

    void   *    p;

    // Pad memory at the end with another alignment bytes to avoid false sharing
    bytes += alignment;

#if HAVE__MM_MALLOC
    p = _mm_malloc( bytes, alignment );
#elif HAVE_POSIX_MEMALIGN
    if (posix_memalign(&p, alignment, bytes) != 0)
    {
        p = nullptr;
    }
#elif HAVE_MEMALIGN
    p = memalign(alignment, bytes);
#elif HAVE__ALIGNED_MALLOC
    p = _aligned_malloc(bytes, alignment);
#else
    p = alignedMallocGeneric(bytes, alignment);
#endif

    return p;
}

void
alignedFree(void *p)
{
    if (p)
    {
#if HAVE__MM_MALLOC
        _mm_free(p);
#elif HAVE_POSIX_MEMALIGN || HAVE_MEMALIGN
        free(p);
#elif HAVE__ALIGNED_MALLOC
        _aligned_free(p);
#else
        alignedFreeGeneric(p);
#endif
    }
}

} // namespace internal

} // namespace gmx
