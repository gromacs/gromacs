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
 * \brief
 * C++ wrapper for analysis tool neighborhood searching.
 */
#ifndef GMX_TRAJECTORYANALYSIS_NBSEARCH_H
#define GMX_TRAJECTORYANALYSIS_NBSEARCH_H

#include "../selection/nbsearch.h"

namespace gmx
{

/*! \brief
 * C++ wrapper for neighborhood searching in nbsearch.h.
 */
class NeighborhoodSearch
{
    public:
        ~NeighborhoodSearch() { gmx_ana_nbsearch_free(_d); }

        static int create(NeighborhoodSearch **nbp, real cutoff, int maxn)
        {
            gmx_ana_nbsearch_t *d;
            int rc = gmx_ana_nbsearch_create(&d, cutoff, maxn);
            if (rc != 0)
            {
                *nbp = NULL;
                return rc;
            }
            *nbp = new NeighborhoodSearch(d);
            return 0;
        }

        int init(t_pbc *pbc, int n, const rvec x[])
        { return gmx_ana_nbsearch_init(_d, pbc, n, x); }

        int init(t_pbc *pbc, const gmx_ana_pos_t *p)
        { return gmx_ana_nbsearch_pos_init(_d, pbc, p); }

        int setExclusions(int nexcl, atom_id *excl)
        { return gmx_ana_nbsearch_set_excl(_d, nexcl, excl); }


        bool isWithin(const rvec x)
        { return gmx_ana_nbsearch_is_within(_d, x); }

        bool isWithin(const gmx_ana_pos_t *p, int i)
        { return gmx_ana_nbsearch_pos_is_within(_d, p, i); }

        real minimumDistance(const rvec x)
        { return gmx_ana_nbsearch_mindist(_d, x); }

        real minimumDistance(const gmx_ana_pos_t *p, int i)
        { return gmx_ana_nbsearch_pos_mindist(_d, p, i); }

        bool firstWithin(const rvec x, int *jp)
        { return gmx_ana_nbsearch_first_within(_d, x, jp); }

        bool firstWithin(const gmx_ana_pos_t *p, int i, int *jp)
        { return gmx_ana_nbsearch_pos_first_within(_d, p, i, jp); }

        bool nextWithin(int *jp)
        { return gmx_ana_nbsearch_next_within(_d, jp); }

    private:
        NeighborhoodSearch(gmx_ana_nbsearch_t *d) : _d(d) {}

        gmx_ana_nbsearch_t  *_d;
};

} // namespace gmx

#endif
