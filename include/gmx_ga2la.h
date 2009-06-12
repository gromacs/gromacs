/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * 
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */
#ifndef _gmx_ga2la_h
#define _gmx_ga2la_h

#ifdef HAVE_CONFIG_H
#include <config.h>
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
    bool      bAll;
    int       mod;
    int       nalloc;
    gmx_laa_t *laa;
    gmx_lal_t *lal;
    int       start_space_search;
} t_gmx_ga2la;

/* Clear all the entries in the ga2la list */
static void ga2la_clear(gmx_ga2la_t ga2la)
{
    int i;

    if (ga2la->bAll)
    {
        for(i=0; i<ga2la->nalloc; i++)
        {
            ga2la->laa[i].cell = -1;
        }
    }
    else
    {
        for(i=0; i<ga2la->nalloc; i++)
        {
            ga2la->lal[i].ga   = -1;
            ga2la->lal[i].next = -1;
        }
        ga2la->start_space_search = ga2la->mod;
    }
}

static gmx_ga2la_t ga2la_init(int nat_tot,int nat_loc)
{
    gmx_ga2la_t ga2la;

    snew(ga2la,1);

    /* There are two methods implemented for finding the local atom number
     * belonging to a global atom number:
     * 1) a simple, direct arrary
     * 2) a list of linked lists indexed with the global number modulo mod.
     * We use whichever method uses less memory.
     * The direct array indexing method is faster for small systems
     * and small systems should (hopefully) have more than a quarter
     * of the atoms locally (local means home + communicated atoms).
     */
    ga2la->bAll = (4*nat_loc >= nat_tot);
    if (ga2la->bAll)
    {
        ga2la->nalloc = nat_tot;
        snew(ga2la->laa,ga2la->nalloc);
    }
    else
    {
        /* Make the direct list twice as large as the number of local atoms.
         * Then most linked lists should have size 0 or 1.
         */
        ga2la->mod = 2*nat_loc;
        ga2la->nalloc = over_alloc_dd(ga2la->mod);
        snew(ga2la->lal,ga2la->nalloc);
    }

    ga2la_clear(ga2la);

    return ga2la;
}

/* Set the ga2la entry for global atom a_gl to local atom a_loc and cell. */
static void ga2la_set(gmx_ga2la_t ga2la,int a_gl,int a_loc,int cell)
{
    int ind,ind_prev,i;

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
        while(ga2la->lal[ind_prev].next >= 0)
        {
            ind_prev = ga2la->lal[ind_prev].next;
        }
        /* Search for space in the array */
        ind = ga2la->start_space_search;
        while (ind < ga2la->nalloc && ga2la->lal[ind].ga >= 0)
        {
            ind++;
        }
        /* If we are a the end of the list we need to increase the size */
        if (ind == ga2la->nalloc)
        {
            ga2la->nalloc = over_alloc_dd(ind+1);
            srenew(ga2la->lal,ga2la->nalloc);
            for(i=ind; i<ga2la->nalloc; i++)
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
static void ga2la_del(gmx_ga2la_t ga2la,int a_gl)
{
    int ind,ind_prev;

    if (ga2la->bAll)
    {
        ga2la->laa[a_gl].cell = -1;
        
        return;
    }

    ind_prev = -1;
    ind = a_gl % ga2la->mod;
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
        ind = ga2la->lal[ind].next;
    }
    while (ind >= 0);

    return;
}

/* Change the local atom for present ga2la entry for global atom a_gl */
static void ga2la_change_la(gmx_ga2la_t ga2la,int a_gl,int a_loc)
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
static bool ga2la_get(const gmx_ga2la_t ga2la,int a_gl,int *a_loc,int *cell)
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
static bool ga2la_home(const gmx_ga2la_t ga2la,int a_gl,int *a_loc)
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

#endif /* _gmx_ga2la_h */
