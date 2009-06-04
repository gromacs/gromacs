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
/*! \file
 * \brief Implementation of functions in position.h.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

#include <smalloc.h>
#include <typedefs.h>
#include <vec.h>

#include <indexutil.h>
#include <position.h>

/*!
 * \param[out] pos      Output structure.
 *
 * Any contents of \p pos are discarded without freeing.
 */
void
gmx_ana_pos_clear(gmx_ana_pos_t *pos)
{
    pos->nr = 0;
    pos->x  = NULL;
    gmx_ana_indexmap_clear(&pos->m);
    pos->g  = NULL;
}

/*!
 * \param[in,out] pos   Position data structure.
 * \param[in]     n     Maximum number of positions.
 * \param[in]     isize Maximum number of atoms.
 *
 * Ensures that enough memory is allocated in \p pos to calculate \p n
 * positions from \p isize atoms.
 */
void
gmx_ana_pos_reserve(gmx_ana_pos_t *pos, int n, int isize)
{
    if (pos->nr < n)
    {
        pos->nr = n;
        srenew(pos->x, n);
    }
    if (isize > 0)
    {
        gmx_ana_indexmap_reserve(&pos->m, n, isize);
    }
}

/*!
 * \param[out]    pos  Position data structure to initialize.
 * \param[in]     x    Position vector to use.
 */
void
gmx_ana_pos_init_const(gmx_ana_pos_t *pos, rvec x)
{
    gmx_ana_pos_clear(pos);
    pos->nr = 1;
    snew(pos->x, 1);
    copy_rvec(x, pos->x[0]);
    gmx_ana_indexmap_init(&pos->m, NULL, NULL, INDEX_UNKNOWN);
    pos->g  = NULL;
}

/*!
 * \param[in,out] pos   Position data structure.
 *
 * Frees any memory allocated within \p pos.
 * The pointer \p pos itself is not freed.
 *
 * \see gmx_ana_pos_free()
 */
void
gmx_ana_pos_deinit(gmx_ana_pos_t *pos)
{
    pos->nr = 0;
    sfree(pos->x); pos->x = NULL;
    gmx_ana_indexmap_deinit(&pos->m);
}

/*!
 * \param[in,out] pos   Position data structure.
 *
 * Frees any memory allocated for \p pos.
 * The pointer \p pos is also freed, and is invalid after the call.
 *
 * \see gmx_ana_pos_deinit()
 */
void
gmx_ana_pos_free(gmx_ana_pos_t *pos)
{
    gmx_ana_pos_deinit(pos);
    sfree(pos);
}

/*!
 * \param[in,out] dest   Destination positions.
 * \param[in]     src    Source positions.
 * \param[in]     bFirst If TRUE, memory is allocated for \p dest and a full
 *   copy is made; otherwise, only variable parts are copied, and no memory
 *   is allocated.
 *
 * \p dest should have been initialized somehow (calloc() is enough).
 */
void
gmx_ana_pos_copy(gmx_ana_pos_t *dest, gmx_ana_pos_t *src, bool bFirst)
{
    if (bFirst)
    {
        gmx_ana_pos_reserve(dest, src->nr, 0);
    }
    dest->nr = src->nr;
    memcpy(dest->x, src->x, dest->nr*sizeof(*dest->x));
    gmx_ana_indexmap_copy(&dest->m, &src->m, bFirst);
    dest->g = src->g;
}
