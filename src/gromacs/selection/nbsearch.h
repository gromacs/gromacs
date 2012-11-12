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
 * \brief API for neighborhood searching.
 *
 * The API is documented in more detail on a separate page:
 * \ref nbsearch
 *
 * The functions within this file can be used independently of the other parts
 * of the library.
 * The library also uses the functions internally.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_NBSEARCH_H
#define GMX_SELECTION_NBSEARCH_H

#include "../legacyheaders/typedefs.h"

#include "indexutil.h"

struct gmx_ana_pos_t;

/** Data structure for neighborhood searches. */
typedef struct gmx_ana_nbsearch_t gmx_ana_nbsearch_t;

/** Create a new neighborhood search data structure. */
gmx_ana_nbsearch_t *
gmx_ana_nbsearch_create(real cutoff, int maxn);
/** Free memory allocated for neighborhood search. */
void
gmx_ana_nbsearch_free(gmx_ana_nbsearch_t *d);

/** Initializes neighborhood search for a new frame. */
void
gmx_ana_nbsearch_init(gmx_ana_nbsearch_t *d, t_pbc *pbc, int n, const rvec x[]);
/** Initializes neighborhood search for a frame using \c gmx_ana_pos_t.  */
void
gmx_ana_nbsearch_pos_init(gmx_ana_nbsearch_t *d, t_pbc *pbc,
                          const struct gmx_ana_pos_t *p);
/** Sets the exclusions for the next neighborhood search. */
void
gmx_ana_nbsearch_set_excl(gmx_ana_nbsearch_t *d, int nexcl, int excl[]);
/** Check whether a point is within a neighborhood. */
bool
gmx_ana_nbsearch_is_within(gmx_ana_nbsearch_t *d, const rvec x);
/** Check whether a position is within a neighborhood. */
bool
gmx_ana_nbsearch_pos_is_within(gmx_ana_nbsearch_t *d,
                               const struct gmx_ana_pos_t *p, int i);
/** Calculates the minimun distance from the reference points. */
real
gmx_ana_nbsearch_mindist(gmx_ana_nbsearch_t *d, const rvec x);
/** Calculates the minimun distance from the reference points. */
real
gmx_ana_nbsearch_pos_mindist(gmx_ana_nbsearch_t *d,
                             const struct gmx_ana_pos_t *p, int i);
/** Finds the first reference position within the cutoff. */
bool
gmx_ana_nbsearch_first_within(gmx_ana_nbsearch_t *d, const rvec x, int *jp);
/** Finds the first reference position within the cutoff. */
bool
gmx_ana_nbsearch_pos_first_within(gmx_ana_nbsearch_t *d,
                                  const struct gmx_ana_pos_t *p, int i, int *jp);
/** Finds the next reference position within the cutoff. */
bool
gmx_ana_nbsearch_next_within(gmx_ana_nbsearch_t *d, int *jp);

namespace gmx
{

/*
 * C++ wrapper for neighborhood searching.
 *
 */
class NeighborhoodSearch
{
    public:
        NeighborhoodSearch(real cutoff, int maxn)
            : d_(gmx_ana_nbsearch_create(cutoff, maxn))
        {
        }
        ~NeighborhoodSearch() { gmx_ana_nbsearch_free(d_); }

        void init(t_pbc *pbc, int n, const rvec x[])
        { gmx_ana_nbsearch_init(d_, pbc, n, x); }

        void init(t_pbc *pbc, const gmx_ana_pos_t *p)
        { gmx_ana_nbsearch_pos_init(d_, pbc, p); }

        void setExclusions(int nexcl, atom_id *excl)
        { gmx_ana_nbsearch_set_excl(d_, nexcl, excl); }


        bool isWithin(const rvec x)
        { return gmx_ana_nbsearch_is_within(d_, x); }

        bool isWithin(const gmx_ana_pos_t *p, int i)
        { return gmx_ana_nbsearch_pos_is_within(d_, p, i); }

        real minimumDistance(const rvec x)
        { return gmx_ana_nbsearch_mindist(d_, x); }

        real minimumDistance(const gmx_ana_pos_t *p, int i)
        { return gmx_ana_nbsearch_pos_mindist(d_, p, i); }

        bool firstWithin(const rvec x, int *jp)
        { return gmx_ana_nbsearch_first_within(d_, x, jp); }

        bool firstWithin(const gmx_ana_pos_t *p, int i, int *jp)
        { return gmx_ana_nbsearch_pos_first_within(d_, p, i, jp); }

        bool nextWithin(int *jp)
        { return gmx_ana_nbsearch_next_within(d_, jp); }

    private:
        gmx_ana_nbsearch_t *d_;
};

} // namespace gmx

#endif
