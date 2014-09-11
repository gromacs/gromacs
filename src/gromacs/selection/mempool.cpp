/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2014, by the GROMACS development team, led by
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
 * Implements functions in mempool.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "mempool.h"

#include <stdlib.h>

#include <new>

#include "gromacs/selection/indexutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

//! Alignment in bytes for all returned blocks.
#define ALIGN_STEP 8

/*! \internal \brief
 * Describes a single block allocated from the memory pool.
 */
typedef struct gmx_sel_mempool_block_t
{
    //! Pointer to the start of the block (as returned to the user).
    void                       *ptr;
    //! Size of the block, including padding required to align next block.
    size_t                      size;
} gmx_sel_mempool_block_t;

/*! \internal \brief
 * Describes a memory pool.
 */
struct gmx_sel_mempool_t
{
    //! Number of bytes currently allocated from the pool.
    size_t                      currsize;
    //! Number of bytes free in the pool, or 0 if \a buffer is NULL.
    size_t                      freesize;
    //! Memory area allocated for the pool, or NULL if not yet reserved.
    char                       *buffer;
    //! Pointer to the first free byte (aligned at ::ALIGN_STEP) in \a buffer.
    char                       *freeptr;
    //! Number of blocks allocated from the pool.
    int                         nblocks;
    //! Array describing the allocated blocks.
    gmx_sel_mempool_block_t    *blockstack;
    //! Number of elements allocated for the \a blockstack array.
    int                         blockstack_nalloc;
    /*! \brief
     * Maximum number of bytes that have been reserved from the pool
     * simultaneously.
     */
    size_t                      maxsize;
};

gmx_sel_mempool_t *
_gmx_sel_mempool_create()
{
    gmx_sel_mempool_t *mp;

    snew(mp, 1);
    mp->currsize          = 0;
    mp->freesize          = 0;
    mp->buffer            = NULL;
    mp->freeptr           = NULL;
    mp->nblocks           = 0;
    mp->blockstack        = NULL;
    mp->blockstack_nalloc = 0;
    mp->maxsize           = 0;
    return mp;
}

void
_gmx_sel_mempool_destroy(gmx_sel_mempool_t *mp)
{
    if (!mp->buffer)
    {
        int  i;

        for (i = 0; i < mp->nblocks; ++i)
        {
            sfree(mp->blockstack[i].ptr);
        }
    }
    sfree(mp->buffer);
    sfree(mp->blockstack);
    sfree(mp);
}

void *
_gmx_sel_mempool_alloc(gmx_sel_mempool_t *mp, size_t size)
{
    void   *ptr = NULL;
    size_t  size_walign;

    size_walign = ((size + ALIGN_STEP - 1) / ALIGN_STEP) * ALIGN_STEP;
    if (mp->buffer)
    {
        if (mp->freesize < size)
        {
            GMX_THROW(gmx::InternalError("Out of memory pool memory"));
        }
        ptr           = mp->freeptr;
        mp->freeptr  += size_walign;
        mp->freesize -= size_walign;
        mp->currsize += size_walign;
    }
    else
    {
        ptr = malloc(size);
        if (!ptr)
        {
            throw std::bad_alloc();
        }
        mp->currsize += size_walign;
        if (mp->currsize > mp->maxsize)
        {
            mp->maxsize = mp->currsize;
        }
    }

    if (mp->nblocks >= mp->blockstack_nalloc)
    {
        mp->blockstack_nalloc = mp->nblocks + 10;
        srenew(mp->blockstack, mp->blockstack_nalloc);
    }
    mp->blockstack[mp->nblocks].ptr  = ptr;
    mp->blockstack[mp->nblocks].size = size_walign;
    mp->nblocks++;

    return ptr;
}

void
_gmx_sel_mempool_free(gmx_sel_mempool_t *mp, void *ptr)
{
    int size;

    if (ptr == NULL)
    {
        return;
    }
    GMX_RELEASE_ASSERT(mp->nblocks > 0 && mp->blockstack[mp->nblocks - 1].ptr == ptr,
                       "Invalid order of memory pool free calls");
    mp->nblocks--;
    size          = mp->blockstack[mp->nblocks].size;
    mp->currsize -= size;
    if (mp->buffer)
    {
        mp->freeptr   = (char *)ptr;
        mp->freesize += size;
    }
    else
    {
        sfree(ptr);
    }
}

void
_gmx_sel_mempool_reserve(gmx_sel_mempool_t *mp, size_t size)
{
    GMX_RELEASE_ASSERT(mp->nblocks == 0,
                       "Cannot reserve memory pool when there is something allocated");
    GMX_RELEASE_ASSERT(!mp->buffer, "Cannot reserve memory pool twice");
    if (size == 0)
    {
        size = mp->maxsize;
    }
    mp->buffer = (char *)malloc(size);
    if (!mp->buffer)
    {
        throw std::bad_alloc();
    }
    mp->freesize = size;
    mp->freeptr  = mp->buffer;
}

void
_gmx_sel_mempool_alloc_group(gmx_sel_mempool_t *mp, gmx_ana_index_t *g,
                             int isize)
{
    void *ptr = _gmx_sel_mempool_alloc(mp, sizeof(*g->index)*isize);
    g->index = static_cast<int *>(ptr);
}

void
_gmx_sel_mempool_free_group(gmx_sel_mempool_t *mp, gmx_ana_index_t *g)
{
    _gmx_sel_mempool_free(mp, g->index);
    g->index = NULL;
}
