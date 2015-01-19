/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 * \brief This file defines functions for managing threading of listed
 * interactions.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_listed-forces
 */
#include "gmxpre.h"

#include "manage-threading.h"

#include <assert.h>

#include <algorithm>

#include "gromacs/listed-forces/listed-forces.h"
#include "gromacs/mdlib/forcerec-threading.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

//! Divides listed interactions over threads
static void divide_bondeds_over_threads(t_idef *idef, int nthreads)
{
    int ftype;
    int nat1;
    int t;
    int il_nr_thread;

    idef->nthreads = nthreads;

    if (F_NRE*(nthreads+1) > idef->il_thread_division_nalloc)
    {
        idef->il_thread_division_nalloc = F_NRE*(nthreads+1);
        snew(idef->il_thread_division, idef->il_thread_division_nalloc);
    }

    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (ftype_is_bonded_potential(ftype))
        {
            nat1 = interaction_function[ftype].nratoms + 1;

            for (t = 0; t <= nthreads; t++)
            {
                /* Divide the interactions equally over the threads.
                 * When the different types of bonded interactions
                 * are distributed roughly equally over the threads,
                 * this should lead to well localized output into
                 * the force buffer on each thread.
                 * If this is not the case, a more advanced scheme
                 * (not implemented yet) will do better.
                 */
                il_nr_thread = (((idef->il[ftype].nr/nat1)*t)/nthreads)*nat1;

                /* Ensure that distance restraint pairs with the same label
                 * end up on the same thread.
                 * This is slighlty tricky code, since the next for iteration
                 * may have an initial il_nr_thread lower than the final value
                 * in the previous iteration, but this will anyhow be increased
                 * to the approriate value again by this while loop.
                 */
                while (ftype == F_DISRES &&
                       il_nr_thread > 0 &&
                       il_nr_thread < idef->il[ftype].nr &&
                       idef->iparams[idef->il[ftype].iatoms[il_nr_thread]].disres.label ==
                       idef->iparams[idef->il[ftype].iatoms[il_nr_thread-nat1]].disres.label)
                {
                    il_nr_thread += nat1;
                }

                idef->il_thread_division[ftype*(nthreads+1)+t] = il_nr_thread;
            }
        }
    }
}

//! Construct a reduction mask for which interaction was computed on which thread
static void
calc_bonded_reduction_mask(gmx_bitmask_t *mask,
                           const t_idef *idef,
                           int shift,
                           int t, int nt)
{
    int           ftype, nb, nat1, nb0, nb1, i, a;
    bitmask_clear(mask);

    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (ftype_is_bonded_potential(ftype))
        {
            nb = idef->il[ftype].nr;
            if (nb > 0)
            {
                nat1 = interaction_function[ftype].nratoms + 1;

                /* Divide this interaction equally over the threads.
                 * This is not stored: should match division in calc_bonds.
                 */
                nb0 = idef->il_thread_division[ftype*(nt+1)+t];
                nb1 = idef->il_thread_division[ftype*(nt+1)+t+1];

                for (i = nb0; i < nb1; i += nat1)
                {
                    for (a = 1; a < nat1; a++)
                    {
                        bitmask_set_bit(mask, idef->il[ftype].iatoms[i+a]>>shift);
                    }
                }
            }
        }
    }
}


/*! \brief We divide the force array in a maximum of BITMASK_SIZE (default 32) blocks.
 * Minimum force block reduction size is thus 2^6=64.
 */
const int maxBlockBits = BITMASK_SIZE;

void setup_bonded_threading(t_forcerec   *fr, t_idef *idef)
{
    int t;
    int ctot, c, b;

    assert(fr->nthreads >= 1);

    /* Divide the bonded interaction over the threads */
    divide_bondeds_over_threads(idef, fr->nthreads);

    if (fr->nthreads == 1)
    {
        fr->red_nblock = 0;

        return;
    }

    fr->red_ashift = 6;
    while (fr->natoms_force > (int)(maxBlockBits*(1U<<fr->red_ashift)))
    {
        fr->red_ashift++;
    }
    if (debug)
    {
        fprintf(debug, "bonded force buffer block atom shift %d bits\n",
                fr->red_ashift);
    }

    /* Determine to which blocks each thread's bonded force calculation
     * contributes. Store this is a mask for each thread.
     */
#pragma omp parallel for num_threads(fr->nthreads) schedule(static)
    for (t = 1; t < fr->nthreads; t++)
    {
        calc_bonded_reduction_mask(&fr->f_t[t].red_mask,
                                   idef, fr->red_ashift, t, fr->nthreads);
    }

    /* Determine the maximum number of blocks we need to reduce over */
    fr->red_nblock = 0;
    ctot           = 0;
    for (t = 0; t < fr->nthreads; t++)
    {
        c = 0;
        for (b = 0; b < maxBlockBits; b++)
        {
            if (bitmask_is_set(fr->f_t[t].red_mask, b))
            {
                fr->red_nblock = std::max(fr->red_nblock, b+1);
                c++;
            }
        }
        if (debug)
        {
#if BITMASK_SIZE <= 64 //move into bitmask when it is C++
            std::string flags = gmx::formatString("%x", fr->f_t[t].red_mask);
#else
            std::string flags = gmx::formatAndJoin(fr->f_t[t].red_mask,
                                                   fr->f_t[t].red_mask+BITMASK_ALEN,
                                                   "", gmx::StringFormatter("%x"));
#endif
            fprintf(debug, "thread %d flags %s count %d\n",
                    t, flags.c_str(), c);
        }
        ctot += c;
    }
    if (debug)
    {
        fprintf(debug, "Number of blocks to reduce: %d of size %d\n",
                fr->red_nblock, 1<<fr->red_ashift);
        fprintf(debug, "Reduction density %.2f density/#thread %.2f\n",
                ctot*(1<<fr->red_ashift)/(double)fr->natoms_force,
                ctot*(1<<fr->red_ashift)/(double)(fr->natoms_force*fr->nthreads));
    }
}
