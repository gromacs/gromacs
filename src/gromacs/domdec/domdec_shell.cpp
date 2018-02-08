/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2006,2007,2008,2009,2010,2012,2013,2014, by the GROMACS development team, led by
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
 * while managing interactions involving shells/Drudes.
 *
 * \author Justin Lemkul <jalemkul@vt.edu> 
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "domdec_shell.h"

#include <assert.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "domdec_specatomcomm.h"
#include "hash.h"

void dd_move_f_shells(gmx_domdec_t *dd, rvec *f, rvec *fshift)
{
    if (dd->shell_comm)
    {
        dd_move_f_specat(dd, dd->shell_comm, f, fshift);
    }
}

void dd_clear_f_shells(gmx_domdec_t *dd, rvec *f)
{
    int i;

    if (dd->shell_comm)
    {
        for (i = dd->shell_comm->at_start; i < dd->shell_comm->at_end; i++)
        {
            clear_rvec(f[i]);
        }
    }
}

void dd_move_x_shells(gmx_domdec_t *dd, matrix box, rvec *x0)
{
    if (dd->shell_comm)
    {
        dd_move_x_specat(dd, dd->shell_comm, box, x0, NULL, FALSE);
    }
}

void dd_move_v_shells(gmx_domdec_t *dd, rvec *v)
{
    if (dd->shell_comm)
    {
        dd_move_v_specat(dd, dd->shell_comm, v);
    }
}

/* This routine sets up proper communication for function types related to
 * polarizable interactions. */
int dd_make_local_shells(gmx_domdec_t *dd, int at_start, t_ilist *lil)
{
    gmx_domdec_specat_comm_t   *spac;
    ind_req_t                  *ireq;
    gmx_hash_t                 *ga2la_specat;
    int                         ftype, nral, i, j, a;
    int                         at_end;
    t_ilist                    *lilf;
    t_iatom                    *iatoms;

    spac = dd->shell_comm;
    ireq = &spac->ireq[0];
    ga2la_specat = dd->ga2la_shell;

    ireq->n = 0;

    /* All Drude bonded functions (anisotropy, Thole, etc) are
     * considered here. Double-counting is avoided using checks
     * within the loop below. */
    for (ftype = F_DRUDEBONDS; ftype <= F_HYPER_POL; ftype++)
    {
        /* we only care about Drudes/shells within this loop, and
         * by choosing these functions we can easily assume the atom
         * indices without any complex decisions within this block of code */
        nral = NRAL(ftype);
        lilf = &lil[ftype];

        for (i = 0; i < lilf->nr; i += 1+nral)
        {
            iatoms = lilf->iatoms + i;

            /* check for other atoms */
            for (j = 1; j < 1 + nral; j++)
            {
                if (iatoms[j] < 0)
                {
                    if (debug)
                    {
                        fprintf(debug, "DD LOCAL SHELLS: non-local iatoms[%d] = %d in %s (i = %d)\n",
                                j, iatoms[j], interaction_function[ftype].longname, i);
                    }
                    /* not a home atom */
                    a = -iatoms[j] - 1;
                    /* returns -1 if key not present in hash */
                    if (gmx_hash_get_minone(dd->ga2la_shell, a) == -1)
                    {
                        /* add non-home atom to list */
                        if (ireq->n+1 > ireq->nalloc)
                        {
                            ireq->nalloc = over_alloc_large(ireq->n+1);
                            srenew(ireq->ind, ireq->nalloc);
                        }
                        ireq->ind[ireq->n++] = a;
                        /* temporarily mark -2, get index later */
                        gmx_hash_set(ga2la_specat, a, -2);
                    }    
                }
            }
        }
    }

    at_end = setup_specat_communication(dd, ireq, dd->shell_comm, ga2la_specat, at_start, 1, "shell", "");

    /* fill in missing indices */
    for (ftype = F_DRUDEBONDS; ftype <= F_HYPER_POL; ftype++)
    {
        /* same as above */
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
                    if (debug)
                    {
                        fprintf(debug, "DD LOCAL SHELLS: after gmx_hash_get_minone iatoms[%d] = %d\n",
                                j, iatoms[j]);
                    }
                }
                /* TODO: fill in shellfc data here. May need to assume order of atoms,
                 * and at this point we can probably only support shell.nnucl = 1 */
            }
        }
    }

    return at_end;
}

void init_domdec_shells(gmx_domdec_t *dd, int n_intercg_shells)
{
    if (debug)
    {
        fprintf(debug, "Begin init_domdec_shells for %d shells\n", n_intercg_shells);
    }

    /* Use a hash table for the global to local index.
     * The number of keys is a rough estimate, it will be optimized later.
     */
    dd->ga2la_shell = gmx_hash_init(std::min(n_intercg_shells/20,
                                             n_intercg_shells/(2*dd->nnodes)));

    dd->shell_comm = specat_comm_init(1);
}

