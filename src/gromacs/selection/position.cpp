/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Implements functions in position.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "position.h"

#include <string.h>

#include "gromacs/math/vec.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

gmx_ana_pos_t::gmx_ana_pos_t()
{
    x = NULL;
    v = NULL;
    f = NULL;
    gmx_ana_indexmap_clear(&m);
    nalloc_x = 0;
}

gmx_ana_pos_t::~gmx_ana_pos_t()
{
    sfree(x);
    sfree(v);
    sfree(f);
    gmx_ana_indexmap_deinit(&m);
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
    GMX_RELEASE_ASSERT(n >= 0, "Invalid position allocation count");
    // Always reserve at least one entry to make NULL checks against pos->x
    // and gmx_ana_pos_reserve_velocities/forces() work as expected in the case
    // that there are actually no positions.
    if (n == 0)
    {
        n = 1;
    }
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
    if (isize >= 0)
    {
        gmx_ana_indexmap_reserve(&pos->m, n, isize);
    }
}

/*!
 * \param[in,out] pos   Position data structure.
 *
 * Currently, this function can only be called after gmx_ana_pos_reserve()
 * has been called at least once with a \p n >= 0.
 */
void
gmx_ana_pos_reserve_velocities(gmx_ana_pos_t *pos)
{
    GMX_RELEASE_ASSERT(pos->nalloc_x > 0,
                       "No memory reserved yet for positions");
    if (!pos->v)
    {
        snew(pos->v, pos->nalloc_x);
    }
}

/*!
 * \param[in,out] pos   Position data structure.
 *
 * Currently, this function can only be called after gmx_ana_pos_reserve()
 * has been called at least once with a \p n >= 0.
 */
void
gmx_ana_pos_reserve_forces(gmx_ana_pos_t *pos)
{
    GMX_RELEASE_ASSERT(pos->nalloc_x > 0,
                       "No memory reserved yet for positions");
    if (!pos->f)
    {
        snew(pos->f, pos->nalloc_x);
    }
}

/*!
 * \param[in,out] pos   Position data structure.
 * \param[in]     n     Maximum number of positions.
 * \param[in]     isize Maximum number of atoms.
 * \param[in]     bVelocities Whether to reserve space for velocities.
 * \param[in]     bForces     Whether to reserve space for forces.
 *
 * Ensures that enough memory is allocated in \p pos to calculate \p n
 * positions from \p isize atoms.
 *
 * This method needs to be called instead of gmx_ana_pos_reserve() if the
 * intent is to use gmx_ana_pos_append_init()/gmx_ana_pos_append().
 */
void
gmx_ana_pos_reserve_for_append(gmx_ana_pos_t *pos, int n, int isize,
                               bool bVelocities, bool bForces)
{
    gmx_ana_pos_reserve(pos, n, isize);
    snew(pos->m.mapb.a, isize);
    pos->m.mapb.nalloc_a = isize;
    if (bVelocities)
    {
        gmx_ana_pos_reserve_velocities(pos);
    }
    if (bForces)
    {
        gmx_ana_pos_reserve_forces(pos);
    }
}

/*!
 * \param[out]    pos  Position data structure to initialize.
 * \param[in]     x    Position vector to use.
 */
void
gmx_ana_pos_init_const(gmx_ana_pos_t *pos, const rvec x)
{
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
 * \param[in,out] dest   Destination positions.
 * \param[in]     src    Source positions.
 * \param[in]     bFirst If true, memory is allocated for \p dest and a full
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
        gmx_ana_pos_reserve(dest, src->count(), -1);
        if (src->v)
        {
            gmx_ana_pos_reserve_velocities(dest);
        }
        if (src->f)
        {
            gmx_ana_pos_reserve_forces(dest);
        }
    }
    memcpy(dest->x, src->x, src->count()*sizeof(*dest->x));
    if (dest->v)
    {
        GMX_ASSERT(src->v, "src velocities should be non-null if dest velocities are allocated");
        memcpy(dest->v, src->v, src->count()*sizeof(*dest->v));
    }
    if (dest->f)
    {
        GMX_ASSERT(src->f, "src forces should be non-null if dest forces are allocated");
        memcpy(dest->f, src->f, src->count()*sizeof(*dest->f));
    }
    gmx_ana_indexmap_copy(&dest->m, &src->m, bFirst);
}

/*!
 * \param[in,out] pos  Position data structure.
 * \param[in]     nr   Number of positions.
 */
void
gmx_ana_pos_set_nr(gmx_ana_pos_t *pos, int nr)
{
    // TODO: This puts the mapping in a somewhat inconsistent state.
    pos->m.mapb.nr = nr;
}

/*!
 * \param[in,out] pos   Position data structure.
 *
 * Sets the number of positions to 0.
 */
void
gmx_ana_pos_empty_init(gmx_ana_pos_t *pos)
{
    pos->m.mapb.nr  = 0;
    pos->m.mapb.nra = 0;
    pos->m.b.nr     = 0;
    pos->m.b.nra    = 0;
    /* Initializing these should not really be necessary, but do it for
     * safety... */
    pos->m.mapb.index[0] = 0;
    pos->m.b.index[0]    = 0;
    /* This function should only be used to construct all the possible
     * positions, so the result should always be static. */
    pos->m.bStatic       = true;
}

/*!
 * \param[in,out] pos   Position data structure.
 *
 * Sets the number of positions to 0.
 */
void
gmx_ana_pos_empty(gmx_ana_pos_t *pos)
{
    pos->m.mapb.nr  = 0;
    pos->m.mapb.nra = 0;
    /* This should not really be necessary, but do it for safety... */
    pos->m.mapb.index[0] = 0;
    /* We set the flag to true, although really in the empty state it
     * should be false. This makes it possible to update the flag in
     * gmx_ana_pos_append(), and just make a simple check in
     * gmx_ana_pos_append_finish(). */
    pos->m.bStatic       = true;
}

/*!
 * \param[in,out] dest  Data structure to which the new position is appended.
 * \param[in]     src   Data structure from which the position is copied.
 * \param[in]     i     Index in \p from to copy.
 */
void
gmx_ana_pos_append_init(gmx_ana_pos_t *dest, gmx_ana_pos_t *src, int i)
{
    int  j, k;

    j = dest->count();
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
        dest->m.mapb.a[dest->m.mapb.nra++] = src->m.mapb.a[k];
        dest->m.b.a[dest->m.b.nra++]       = src->m.b.a[k];
    }
    dest->m.mapb.index[j+1] = dest->m.mapb.nra;
    dest->m.b.index[j+1]    = dest->m.mapb.nra;
    dest->m.mapb.nr++;
    dest->m.b.nr++;
}

/*!
 * \param[in,out] dest  Data structure to which the new position is appended.
 * \param[in]     src   Data structure from which the position is copied.
 * \param[in]     i     Index in \p src to copy.
 * \param[in]     refid Reference ID in \p out
 *   (all negative values are treated as -1).
 */
void
gmx_ana_pos_append(gmx_ana_pos_t *dest, gmx_ana_pos_t *src, int i, int refid)
{
    for (int k = src->m.mapb.index[i]; k < src->m.mapb.index[i+1]; ++k)
    {
        dest->m.mapb.a[dest->m.mapb.nra++] = src->m.mapb.a[k];
    }
    const int j = dest->count();
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
        dest->m.bStatic  = false;
        /* If we are using masks, there is no need to alter the
         * mapid field. */
    }
    else
    {
        if (refid != j)
        {
            dest->m.bStatic = false;
        }
        dest->m.refid[j] = refid;
        /* Use the original IDs from the output structure to correctly
         * handle user customization. */
        dest->m.mapid[j] = dest->m.orgid[refid];
    }
    dest->m.mapb.index[j+1] = dest->m.mapb.nra;
    dest->m.mapb.nr++;
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
    if (pos->m.mapb.nr != pos->m.b.nr)
    {
        pos->m.bStatic = false;
    }
}

/*!
 * \param[in,out] g     Data structure to which the new atoms are appended.
 * \param[in]     src   Data structure from which the position is copied.
 * \param[in]     i     Index in \p src to copy.
 */
void
gmx_ana_pos_add_to_group(gmx_ana_index_t *g, gmx_ana_pos_t *src, int i)
{
    for (int k = src->m.mapb.index[i]; k < src->m.mapb.index[i+1]; ++k)
    {
        g->index[g->isize++] = src->m.mapb.a[k];
    }
}
