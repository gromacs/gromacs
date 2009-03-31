/*
 * $Id$
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
/*! \page nbsearch Neighborhood search routines
 *
 * Functions to find particles within a neighborhood of a set of particles
 * are defined in nbsearch.h.
 * The usage is simple: a data structure is allocated with
 * gmx_ana_nbsearch_create(), and the box shape and reference positions for a
 * frame are set using gmx_ana_nbsearch_init() or gmx_ana_nbsearch_pos_init().
 * Searches can then be performed with gmx_ana_nbsearch_is_within() and
 * gmx_ana_nbsearch_mindist(), or with versions that take the \c gmx_ana_pos_t
 * data structure.
 * When the data structure is no longer required, it can be freed with
 * gmx_ana_nbsearch_free().
 *
 * \todo
 * Implement functions for looping through all pairs within a cutoff.
 *
 * \todo
 * Implement grid-based searching
 * (currently everything is implemented using an expensive O(n^2) loop).
 *
 * \todo
 * Implement exclusions.
 */
/*! \file
 * \brief Implementation of functions in nbsearch.h.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <smalloc.h>
#include <typedefs.h>
#include <pbc.h>
#include <vec.h>

#include <nbsearch.h>
#include <position.h>

//! Data structure for neighborhood searches.
struct gmx_ana_nbsearch_t
{
    //! The cutoff squared.
    real           cutoff2;
    //! Maximum number of reference points.
    int            maxnref;
    //! Number of reference points for the current frame.
    int            nref;
    //! Reference point positions.
    rvec          *xref;
    //! PBC data.
    t_pbc         *pbc;
};

/*!
 * \param[out] data   Neighborhood search data structure pointer to initialize.
 * \param[in]  cutoff Cutoff distance for the search
 *   (<=0 stands for no cutoff).
 * \param[in]  maxn   Maximum number of reference particles.
 * \returns    0 on success.
 */
int
gmx_ana_nbsearch_create(gmx_ana_nbsearch_t **data, real cutoff, int maxn)
{
    gmx_ana_nbsearch_t *d;
    
    snew(d, 1);
    if (cutoff <= 0)
    {
        cutoff = HUGE_VAL;
    }
    d->cutoff2 = sqr(cutoff);
    d->maxnref = maxn;
    *data = d;
    return 0;
}

/*!
 * \param     d Data structure to free.
 *
 * After the call, the pointer \p d is no longer valid.
 */
void
gmx_ana_nbsearch_free(gmx_ana_nbsearch_t *d)
{
    sfree(d);
}

/*!
 * \param[in,out] d   Neighborhood search data structure.
 * \param[in]     pbc PBC information for the frame.
 * \param[in]     n   Number of reference positions for the frame.
 * \param[in]     x   \p n reference positions for the frame.
 * \returns       0 on success.
 *
 * Initializes the data structure \p d such that it can be used to search
 * for the neighbors of \p x.
 */
int
gmx_ana_nbsearch_init(gmx_ana_nbsearch_t *d, t_pbc *pbc, int n, rvec x[])
{
    d->pbc  = pbc;
    d->nref = n;
    d->xref = x;
    return 0;
}

/*!
 * \param[in,out] d   Neighborhood search data structure.
 * \param[in]     pbc PBC information for the frame.
 * \param[in]     p   Reference positions for the frame.
 * \returns       0 on success.
 *
 * A convenience wrapper for gmx_ana_nbsearch_init().
 */
int
gmx_ana_nbsearch_pos_init(gmx_ana_nbsearch_t *d, t_pbc *pbc, gmx_ana_pos_t *p)
{
    return gmx_ana_nbsearch_init(d, pbc, p->nr, p->x);
}

/*!
 * \param[in] d   Neighborhood search data structure.
 * \param[in] x   Test position.
 * \returns   TRUE if \p x is within the cutoff of any reference position,
 *   FALSE otherwise.
 */
bool
gmx_ana_nbsearch_is_within(gmx_ana_nbsearch_t *d, rvec x)
{
    int  i;
    rvec dx;

    for (i = 0; i < d->nref; ++i)
    {
        if (d->pbc)
        {
            pbc_dx(d->pbc, x, d->xref[i], dx);
        }
        else
        {
            rvec_sub(x, d->xref[i], dx);
        }
        if (norm2(dx) <= d->cutoff2)
        {
            return TRUE;
        }
    }
    return FALSE;
}

/*!
 * \param[in] d   Neighborhood search data structure.
 * \param[in] p   Test positions.
 * \param[in] i   Use the i'th position in \p p for testing.
 * \returns   TRUE if the test position is within the cutoff of any reference
 *   position, FALSE otherwise.
 */
bool
gmx_ana_nbsearch_pos_is_within(gmx_ana_nbsearch_t *d, gmx_ana_pos_t *p, int i)
{
    return gmx_ana_nbsearch_is_within(d, p->x[i]);
}

/*!
 * \param[in] d   Neighborhood search data structure.
 * \param[in] x   Test position.
 * \returns   The distance to the nearest reference position, or the cutoff
 *   value if there are no reference positions within the cutoff.
 */
real
gmx_ana_nbsearch_mindist(gmx_ana_nbsearch_t *d, rvec x)
{
    int  i;
    rvec dx;
    real d2, mind2;

    mind2 = d->cutoff2;
    for (i = 0; i < d->nref; ++i)
    {
        if (d->pbc)
        {
            pbc_dx(d->pbc, x, d->xref[i], dx);
        }
        else
        {
            rvec_sub(x, d->xref[i], dx);
        }
        d2 = norm2(dx);
        if (d2 < mind2)
        {
            mind2 = d2;
        }
    }
    return sqrt(mind2);
}

/*!
 * \param[in] d   Neighborhood search data structure.
 * \param[in] p   Test positions.
 * \param[in] i   Use the i'th position in \p p for testing.
 * \returns   The distance to the nearest reference position, or the cutoff
 *   value if there are no reference positions within the cutoff.
 */
real
gmx_ana_nbsearch_pos_mindist(gmx_ana_nbsearch_t *d, gmx_ana_pos_t *p, int i)
{
    return gmx_ana_nbsearch_mindist(d, p->x[i]);
}
