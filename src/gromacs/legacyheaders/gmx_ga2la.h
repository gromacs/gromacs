/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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
#ifndef _gmx_ga2la_h
#define _gmx_ga2la_h

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/smalloc.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int  la;
    int  cell;
} gmx_laa_t;

typedef struct {
    int  ga;
    int  la;
    int  cell;
    int  next;
} gmx_lal_t;

typedef struct gmx_ga2la {
    gmx_bool      bAll;
    int           mod;
    int           nalloc;
    gmx_laa_t    *laa;
    gmx_lal_t    *lal;
    int           start_space_search;
} t_gmx_ga2la;

/* Clear all the entries in the ga2la list */
static void ga2la_clear(gmx_ga2la_t ga2la)
{
    int i;

    if (ga2la->bAll)
    {
        for (i = 0; i < ga2la->nalloc; i++)
        {
            ga2la->laa[i].cell = -1;
        }
    }
    else
    {
        for (i = 0; i < ga2la->nalloc; i++)
        {
            ga2la->lal[i].ga   = -1;
            ga2la->lal[i].next = -1;
        }
        ga2la->start_space_search = ga2la->mod;
    }
}

static gmx_ga2la_t ga2la_init(int nat_tot, int nat_loc)
{
    gmx_ga2la_t ga2la;

    snew(ga2la, 1);

    /* There are two methods implemented for finding the local atom number
     * belonging to a global atom number:
     * 1) a simple, direct array
     * 2) a list of linked lists indexed with the global number modulo mod.
     * Memory requirements:
     * 1) nat_tot*2 ints
     * 2) nat_loc*(2+1-2(1-e^-1/2))*4 ints
     * where nat_loc is the number of atoms in the home + communicated zones.
     * Method 1 is faster for low parallelization, 2 for high parallelization.
     * We switch to method 2 when it uses less than half the memory method 1.
     */
    ga2la->bAll = (nat_tot < 1024 || 9*nat_loc >= nat_tot);
    if (ga2la->bAll)
    {
        ga2la->nalloc = nat_tot;
        snew(ga2la->laa, ga2la->nalloc);
    }
    else
    {
        /* Make the direct list twice as long as the number of local atoms.
         * The fraction of entries in the list with:
         * 0   size lists: e^-1/f
         * >=1 size lists: 1 - e^-1/f
         * where f is: the direct list length / #local atoms
         * The fraction of atoms not in the direct list is: 1-f(1-e^-1/f).
         */
        ga2la->mod    = 2*nat_loc;
        ga2la->nalloc = over_alloc_dd(ga2la->mod);
        snew(ga2la->lal, ga2la->nalloc);
    }

    ga2la_clear(ga2la);

    return ga2la;
}

/* Set the ga2la entry for global atom a_gl to local atom a_loc and cell. */
static void ga2la_set(gmx_ga2la_t ga2la, int a_gl, int a_loc, int cell)
{
    int ind, ind_prev, i;

    if (ga2la->bAll)
    {
        ga2la->laa[a_gl].la   = a_loc;
        ga2la->laa[a_gl].cell = cell;

        return;
    }

    ind = a_gl % ga2la->mod;

    if (ga2la->lal[ind].ga >= 0)
    {
        /* Search the last entry in the linked list for this index */
        ind_prev = ind;
        while (ga2la->lal[ind_prev].next >= 0)
        {
            ind_prev = ga2la->lal[ind_prev].next;
        }
        /* Search for space in the array */
        ind = ga2la->start_space_search;
        while (ind < ga2la->nalloc && ga2la->lal[ind].ga >= 0)
        {
            ind++;
        }
        /* If we are at the end of the list we need to increase the size */
        if (ind == ga2la->nalloc)
        {
            ga2la->nalloc = over_alloc_dd(ind+1);
            srenew(ga2la->lal, ga2la->nalloc);
            for (i = ind; i < ga2la->nalloc; i++)
            {
                ga2la->lal[i].ga   = -1;
                ga2la->lal[i].next = -1;
            }
        }
        ga2la->lal[ind_prev].next = ind;

        ga2la->start_space_search = ind + 1;
    }
    ga2la->lal[ind].ga   = a_gl;
    ga2la->lal[ind].la   = a_loc;
    ga2la->lal[ind].cell = cell;
}

/* Delete the ga2la entry for global atom a_gl */
static void ga2la_del(gmx_ga2la_t ga2la, int a_gl)
{
    int ind, ind_prev;

    if (ga2la->bAll)
    {
        ga2la->laa[a_gl].cell = -1;

        return;
    }

    ind_prev = -1;
    ind      = a_gl % ga2la->mod;
    do
    {
        if (ga2la->lal[ind].ga == a_gl)
        {
            if (ind_prev >= 0)
            {
                ga2la->lal[ind_prev].next = ga2la->lal[ind].next;

                /* This index is a linked entry, so we free an entry.
                 * Check if we are creating the first empty space.
                 */
                if (ind < ga2la->start_space_search)
                {
                    ga2la->start_space_search = ind;
                }
            }
            ga2la->lal[ind].ga   = -1;
            ga2la->lal[ind].cell = -1;
            ga2la->lal[ind].next = -1;

            return;
        }
        ind_prev = ind;
        ind      = ga2la->lal[ind].next;
    }
    while (ind >= 0);

    return;
}

/* Change the local atom for present ga2la entry for global atom a_gl */
static void ga2la_change_la(gmx_ga2la_t ga2la, int a_gl, int a_loc)
{
    int ind;

    if (ga2la->bAll)
    {
        ga2la->laa[a_gl].la = a_loc;

        return;
    }

    ind = a_gl % ga2la->mod;
    do
    {
        if (ga2la->lal[ind].ga == a_gl)
        {
            ga2la->lal[ind].la = a_loc;

            return;
        }
        ind = ga2la->lal[ind].next;
    }
    while (ind >= 0);

    return;
}

/* Returns if the global atom a_gl available locally.
 * Sets the local atom and cell,
 * cell can be larger than the number of zones,
 * in which case it indicates that it is more than one cell away
 * in zone cell - #zones.
 */
static gmx_bool ga2la_get(const gmx_ga2la_t ga2la, int a_gl, int *a_loc, int *cell)
{
    int ind;

    if (ga2la->bAll)
    {
        *a_loc = ga2la->laa[a_gl].la;
        *cell  = ga2la->laa[a_gl].cell;

        return (ga2la->laa[a_gl].cell >= 0);
    }

    ind = a_gl % ga2la->mod;
    do
    {
        if (ga2la->lal[ind].ga == a_gl)
        {
            *a_loc = ga2la->lal[ind].la;
            *cell  = ga2la->lal[ind].cell;

            return TRUE;
        }
        ind = ga2la->lal[ind].next;
    }
    while (ind >= 0);

    return FALSE;
}

/* Returns if the global atom a_gl is a home atom.
 * Sets the local atom.
 */
static gmx_bool ga2la_get_home(const gmx_ga2la_t ga2la, int a_gl, int *a_loc)
{
    int ind;

    if (ga2la->bAll)
    {
        *a_loc = ga2la->laa[a_gl].la;

        return (ga2la->laa[a_gl].cell == 0);
    }

    ind = a_gl % ga2la->mod;
    do
    {
        if (ga2la->lal[ind].ga == a_gl)
        {
            if (ga2la->lal[ind].cell == 0)
            {
                *a_loc = ga2la->lal[ind].la;

                return TRUE;
            }
            else
            {
                return FALSE;
            }
        }
        ind = ga2la->lal[ind].next;
    }
    while (ind >= 0);

    return FALSE;
}

/* Returns if the global atom a_gl is a home atom.
 */
static gmx_bool ga2la_is_home(const gmx_ga2la_t ga2la, int a_gl)
{
    int ind;

    if (ga2la->bAll)
    {
        return (ga2la->laa[a_gl].cell == 0);
    }

    ind = a_gl % ga2la->mod;
    do
    {
        if (ga2la->lal[ind].ga == a_gl)
        {
            return (ga2la->lal[ind].cell == 0);
        }
        ind = ga2la->lal[ind].next;
    }
    while (ind >= 0);

    return FALSE;
}

#ifdef __cplusplus
}
#endif

#endif /* _gmx_ga2la_h */
