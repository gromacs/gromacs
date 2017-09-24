/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2006,2007,2008,2009,2010,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 *
 * \brief This file implements functions for domdec to use
 * while managing inter-atomic constraints.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "domdec_vsite.h"

#include <assert.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "domdec_specatomcomm.h"
#include "hash.h"

void dd_move_f_vsites(gmx_domdec_t *dd, rvec *f, rvec *fshift)
{
    if (dd->vsite_comm)
    {
        dd_move_f_specat(dd, dd->vsite_comm, f, fshift);
    }
}

void dd_clear_f_vsites(gmx_domdec_t *dd, rvec *f)
{
    int i;

    if (dd->vsite_comm)
    {
        for (i = dd->vsite_comm->at_start; i < dd->vsite_comm->at_end; i++)
        {
            clear_rvec(f[i]);
        }
    }
}

void dd_move_x_vsites(gmx_domdec_t *dd, const matrix box, rvec *x)
{
    if (dd->vsite_comm)
    {
        dd_move_x_specat(dd, dd->vsite_comm, box, x, nullptr, FALSE);
    }
}

int dd_make_local_vsites(gmx_domdec_t *dd, int at_start, t_ilist *lil)
{
    gmx_domdec_specat_comm_t   *spac;
    ind_req_t                  *ireq;
    gmx_hash_t                 *ga2la_specat;
    int  ftype, nral, i, j, a;
    t_ilist                    *lilf;
    t_iatom                    *iatoms;
    int  at_end;

    spac         = dd->vsite_comm;
    ireq         = &spac->ireq[0];
    ga2la_specat = dd->ga2la_vsite;

    ireq->n = 0;
    /* Loop over all the home vsites */
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nral = NRAL(ftype);
            lilf = &lil[ftype];
            for (i = 0; i < lilf->nr; i += 1+nral)
            {
                iatoms = lilf->iatoms + i;
                /* Check if we have the other atoms */
                for (j = 1; j < 1+nral; j++)
                {
                    if (iatoms[j] < 0)
                    {
                        /* This is not a home atom,
                         * we need to ask our neighbors.
                         */
                        a = -iatoms[j] - 1;
                        /* Check to not ask for the same atom more than once */
                        if (gmx_hash_get_minone(dd->ga2la_vsite, a) == -1)
                        {
                            /* Add this non-home atom to the list */
                            if (ireq->n+1 > ireq->nalloc)
                            {
                                ireq->nalloc = over_alloc_large(ireq->n+1);
                                srenew(ireq->ind, ireq->nalloc);
                            }
                            ireq->ind[ireq->n++] = a;
                            /* Temporarily mark with -2,
                             * we get the index later.
                             */
                            gmx_hash_set(ga2la_specat, a, -2);
                        }
                    }
                }
            }
        }
    }

    at_end = setup_specat_communication(dd, ireq, dd->vsite_comm, ga2la_specat,
                                        at_start, 1, "vsite", "");

    /* Fill in the missing indices */
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nral = NRAL(ftype);
            lilf = &lil[ftype];
            for (i = 0; i < lilf->nr; i += 1+nral)
            {
                iatoms = lilf->iatoms + i;
                for (j = 1; j < 1+nral; j++)
                {
                    if (iatoms[j] < 0)
                    {
                        iatoms[j] = gmx_hash_get_minone(ga2la_specat, -iatoms[j]-1);
                    }
                }
            }
        }
    }

    return at_end;
}

void init_domdec_vsites(gmx_domdec_t *dd, int n_intercg_vsite)
{
    if (debug)
    {
        fprintf(debug, "Begin init_domdec_vsites\n");
    }

    /* Use a hash table for the global to local index.
     * The number of keys is a rough estimate, it will be optimized later.
     */
    dd->ga2la_vsite = gmx_hash_init(std::min(n_intercg_vsite/20,
                                             n_intercg_vsite/(2*dd->nnodes)));

    dd->vsite_comm = specat_comm_init(1);
}
