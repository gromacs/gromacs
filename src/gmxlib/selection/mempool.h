/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief Declarations for memory pooling functions.
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 */
#ifndef GMX_SELECTION_MEMPOOL_H
#define GMX_SELECTION_MEMPOOL_H

struct gmx_ana_index_t;

typedef struct gmx_sel_mempool_t gmx_sel_mempool_t;

/** Create an empty memory pool. */
int
_gmx_sel_mempool_create(gmx_sel_mempool_t **mpp);
/** Destroy a memory pool. */
void
_gmx_sel_mempool_destroy(gmx_sel_mempool_t *mp);

/** Allocate memory from a memory pool. */
int
_gmx_sel_mempool_alloc(gmx_sel_mempool_t *mp, void **ptrp, size_t size);
/** Release memory allocated from a memory pool. */
void
_gmx_sel_mempool_free(gmx_sel_mempool_t *mp, void *ptr);
/** Set the size of a memory pool. */
int
_gmx_sel_mempool_reserve(gmx_sel_mempool_t *mp, size_t size);

/** Convenience function for allocating an index group from a memory pool. */
int
_gmx_sel_mempool_alloc_group(gmx_sel_mempool_t *mp, struct gmx_ana_index_t *g,
                             int isize);
/** Convenience function for freeing an index group from a memory pool. */
void
_gmx_sel_mempool_free_group(gmx_sel_mempool_t *mp, struct gmx_ana_index_t *g);

#endif
