/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * \brief Memory pooling for selection evaluation.
 *
 * \todo
 * Document these functions.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <stdlib.h>

#include <gmx_fatal.h>
#include <smalloc.h>

#include <indexutil.h>

#include "mempool.h"

#define ALIGN_STEP 8

typedef struct gmx_sel_mempool_block_t
{
    void                       *ptr;
    size_t                      size;
} gmx_sel_mempool_block_t;

struct gmx_sel_mempool_t
{
    size_t                      currsize;
    size_t                      freesize;
    char                       *buffer;
    char                       *freeptr;
    int                         nblocks;
    gmx_sel_mempool_block_t    *blockstack;
    int                         blockstack_nalloc;
    size_t                      maxsize;
};

int
_gmx_sel_mempool_create(gmx_sel_mempool_t **mpp)
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
    *mpp                  = mp;
    return 0;
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

int
_gmx_sel_mempool_alloc(gmx_sel_mempool_t *mp, void **ptrp, size_t size)
{
    void   *ptr = NULL;
    size_t  size_walign;

    *ptrp       = NULL;
    size_walign = ((size + ALIGN_STEP - 1) / ALIGN_STEP) * ALIGN_STEP;
    if (mp->buffer)
    {
        if (mp->freesize < size)
        {
            gmx_bug("out of memory pool memory");
            return ENOMEM;
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
            gmx_mem("out of memory");
            return ENOMEM;
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

    *ptrp = ptr;
    return 0;
}

void
_gmx_sel_mempool_free(gmx_sel_mempool_t *mp, void *ptr)
{
    int size;

    if (ptr == NULL)
    {
        return;
    }
    assert(mp->nblocks > 0 && mp->blockstack[mp->nblocks - 1].ptr == ptr);
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

int
_gmx_sel_mempool_reserve(gmx_sel_mempool_t *mp, size_t size)
{
    assert(mp->nblocks == 0 && !mp->buffer);
    if (size == 0)
    {
        size = mp->maxsize;
    }
    mp->buffer = (char *)malloc(size);
    if (!mp->buffer)
    {
        gmx_mem("out of memory");
        return ENOMEM;
    }
    mp->freesize = size;
    mp->freeptr  = mp->buffer;
    return 0;
}

int
_gmx_sel_mempool_alloc_group(gmx_sel_mempool_t *mp, gmx_ana_index_t *g,
                             int isize)
{
    return _gmx_sel_mempool_alloc(mp, (void **)&g->index,
                                  sizeof(*g->index)*isize);
}

void
_gmx_sel_mempool_free_group(gmx_sel_mempool_t *mp, gmx_ana_index_t *g)
{
    _gmx_sel_mempool_free(mp, g->index);
    g->index = NULL;
}
