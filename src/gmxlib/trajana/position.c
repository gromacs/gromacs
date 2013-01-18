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
 * \brief Implementation of functions in position.h.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
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
    pos->v  = NULL;
    pos->f  = NULL;
    gmx_ana_indexmap_clear(&pos->m);
    pos->g        = NULL;
    pos->nalloc_x = 0;
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
    if (pos->nalloc_x < n)
    {
        pos->nalloc_x = n;
        srenew(pos->x, n);
        if (pos->v)
        {
            srenew(pos->v, n);
        }
        if (pos->f)
        {
            srenew(pos->f, n);
        }
    }
    if (isize > 0)
    {
        gmx_ana_indexmap_reserve(&pos->m, n, isize);
    }
}

/*!
 * \param[in,out] pos   Position data structure.
 *
 * Currently, this function can only be called after gmx_ana_pos_reserve()
 * has been called at least once with a \p n > 0.
 */
void
gmx_ana_pos_reserve_velocities(gmx_ana_pos_t *pos)
{
    assert(pos->nalloc_x > 0);
    if (!pos->v)
    {
        snew(pos->v, pos->nalloc_x);
    }
}

/*!
 * \param[in,out] pos   Position data structure.
 *
 * Currently, this function can only be called after gmx_ana_pos_reserve()
 * has been called at least once with a \p n > 0.
 */
void
gmx_ana_pos_reserve_forces(gmx_ana_pos_t *pos)
{
    assert(pos->nalloc_x > 0);
    if (!pos->f)
    {
        snew(pos->f, pos->nalloc_x);
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
    snew(pos->v, 1);
    snew(pos->f, 1);
    pos->nalloc_x = 1;
    copy_rvec(x, pos->x[0]);
    clear_rvec(pos->v[0]);
    clear_rvec(pos->f[0]);
    gmx_ana_indexmap_init(&pos->m, NULL, NULL, INDEX_UNKNOWN);
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
    pos->nr               = 0;
    sfree(pos->x); pos->x = NULL;
    sfree(pos->v); pos->v = NULL;
    sfree(pos->f); pos->f = NULL;
    pos->nalloc_x         = 0;
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
gmx_ana_pos_copy(gmx_ana_pos_t *dest, gmx_ana_pos_t *src, gmx_bool bFirst)
{
    if (bFirst)
    {
        gmx_ana_pos_reserve(dest, src->nr, 0);
        if (src->v)
        {
            gmx_ana_pos_reserve_velocities(dest);
        }
        if (src->f)
        {
            gmx_ana_pos_reserve_forces(dest);
        }
    }
    dest->nr = src->nr;
    memcpy(dest->x, src->x, dest->nr*sizeof(*dest->x));
    if (dest->v)
    {
        assert(src->v);
        memcpy(dest->v, src->v, dest->nr*sizeof(*dest->v));
    }
    if (dest->f)
    {
        assert(src->f);
        memcpy(dest->f, src->f, dest->nr*sizeof(*dest->f));
    }
    gmx_ana_indexmap_copy(&dest->m, &src->m, bFirst);
    dest->g = src->g;
}

/*!
 * \param[in,out] pos  Position data structure.
 * \param[in]     nr   Number of positions.
 */
void
gmx_ana_pos_set_nr(gmx_ana_pos_t *pos, int nr)
{
    pos->nr = nr;
}

/*!
 * \param[in,out] pos  Position data structure.
 * \param         g    Evaluation group.
 *
 * The old group, if any, is discarded.
 * Note that only a pointer to \p g is stored; it is the responsibility of
 * the caller to ensure that \p g is not freed while it can be accessed
 * through \p pos.
 */
void
gmx_ana_pos_set_evalgrp(gmx_ana_pos_t *pos, gmx_ana_index_t *g)
{
    pos->g = g;
}

/*!
 * \param[in,out] pos   Position data structure.
 *
 * Sets the number of positions to 0.
 */
void
gmx_ana_pos_empty_init(gmx_ana_pos_t *pos)
{
    pos->nr        = 0;
    pos->m.nr      = 0;
    pos->m.mapb.nr = 0;
    pos->m.b.nr    = 0;
    pos->m.b.nra   = 0;
    /* This should not really be necessary, but do it for safety... */
    pos->m.mapb.index[0] = 0;
    pos->m.b.index[0]    = 0;
    /* This function should only be used to construct all the possible
     * positions, so the result should always be static. */
    pos->m.bStatic    = TRUE;
    pos->m.bMapStatic = TRUE;
}

/*!
 * \param[in,out] pos   Position data structure.
 *
 * Sets the number of positions to 0.
 */
void
gmx_ana_pos_empty(gmx_ana_pos_t *pos)
{
    pos->nr        = 0;
    pos->m.nr      = 0;
    pos->m.mapb.nr = 0;
    /* This should not really be necessary, but do it for safety... */
    pos->m.mapb.index[0] = 0;
    /* We set the flags to TRUE, although really in the empty state they
     * should be FALSE. This makes it possible to update the flags in
     * gmx_ana_pos_append(), and just make a simple check in
     * gmx_ana_pos_append_finish(). */
    pos->m.bStatic    = TRUE;
    pos->m.bMapStatic = TRUE;
}

/*!
 * \param[in,out] dest  Data structure to which the new position is appended.
 * \param[in,out] g     Data structure to which the new atoms are appended.
 * \param[in]     src   Data structure from which the position is copied.
 * \param[in]     i     Index in \p from to copy.
 */
void
gmx_ana_pos_append_init(gmx_ana_pos_t *dest, gmx_ana_index_t *g,
                        gmx_ana_pos_t *src, int i)
{
    int  j, k;

    j = dest->nr;
    copy_rvec(src->x[i], dest->x[j]);
    if (dest->v)
    {
        if (src->v)
        {
            copy_rvec(src->v[i], dest->v[j]);
        }
        else
        {
            clear_rvec(dest->v[j]);
        }
    }
    if (dest->f)
    {
        if (src->f)
        {
            copy_rvec(src->f[i], dest->f[j]);
        }
        else
        {
            clear_rvec(dest->f[j]);
        }
    }
    dest->m.refid[j] = j;
    dest->m.mapid[j] = src->m.mapid[i];
    dest->m.orgid[j] = src->m.orgid[i];
    for (k = src->m.mapb.index[i]; k < src->m.mapb.index[i+1]; ++k)
    {
        g->index[g->isize++]         = src->g->index[k];
        dest->m.b.a[dest->m.b.nra++] = src->m.b.a[k];
    }
    dest->m.mapb.index[j+1] = g->isize;
    dest->m.b.index[j+1]    = g->isize;
    dest->nr++;
    dest->m.nr      = dest->nr;
    dest->m.mapb.nr = dest->nr;
    dest->m.b.nr    = dest->nr;
}

/*!
 * \param[in,out] dest  Data structure to which the new position is appended
 *      (can be NULL, in which case only \p g is updated).
 * \param[in,out] g     Data structure to which the new atoms are appended.
 * \param[in]     src   Data structure from which the position is copied.
 * \param[in]     i     Index in \p src to copy.
 * \param[in]     refid Reference ID in \p out
 *   (all negative values are treated as -1).
 *
 * If \p dest is NULL, the value of \p refid is not used.
 */
void
gmx_ana_pos_append(gmx_ana_pos_t *dest, gmx_ana_index_t *g,
                   gmx_ana_pos_t *src, int i, int refid)
{
    int  j, k;

    for (k = src->m.mapb.index[i]; k < src->m.mapb.index[i+1]; ++k)
    {
        g->index[g->isize++] = src->g->index[k];
    }
    if (dest)
    {
        j = dest->nr;
        if (dest->v)
        {
            if (src->v)
            {
                copy_rvec(src->v[i], dest->v[j]);
            }
            else
            {
                clear_rvec(dest->v[j]);
            }
        }
        if (dest->f)
        {
            if (src->f)
            {
                copy_rvec(src->f[i], dest->f[j]);
            }
            else
            {
                clear_rvec(dest->f[j]);
            }
        }
        copy_rvec(src->x[i], dest->x[j]);
        if (refid < 0)
        {
            dest->m.refid[j] = -1;
            dest->m.bStatic  = FALSE;
            /* If we are using masks, there is no need to alter the
             * mapid field. */
        }
        else
        {
            if (refid != j)
            {
                dest->m.bStatic    = FALSE;
                dest->m.bMapStatic = FALSE;
            }
            dest->m.refid[j] = refid;
            /* Use the original IDs from the output structure to correctly
             * handle user customization. */
            dest->m.mapid[j] = dest->m.orgid[refid];
        }
        dest->m.mapb.index[j+1] = g->isize;
        dest->nr++;
        dest->m.nr      = dest->nr;
        dest->m.mapb.nr = dest->nr;
    }
}

/*!
 * \param[in,out] pos   Position data structure.
 *
 * After gmx_ana_pos_empty(), internal state of the position data structure
 * is not consistent before this function is called. This function should be
 * called after any gmx_ana_pos_append() calls have been made.
 */
void
gmx_ana_pos_append_finish(gmx_ana_pos_t *pos)
{
    if (pos->m.nr != pos->m.b.nr)
    {
        pos->m.bStatic    = FALSE;
        pos->m.bMapStatic = FALSE;
    }
}
