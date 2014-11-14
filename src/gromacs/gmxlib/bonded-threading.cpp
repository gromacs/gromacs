/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/legacyheaders/bonded-threading.h"

#include <assert.h>

#include <algorithm>

#include "gromacs/listed-forces/bonded.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static void divide_bondeds_over_threads(t_idef *idef, int nthreads)
{
    int      ftype[F_NRE];
    t_ilist *il[F_NRE];
    int      nat[F_NRE], ind[F_NRE];
    int      nat_tot, nat_sum;
    int      f, ntype;
    int      t;

    assert(nthreads > 0);

    idef->nthreads = nthreads;

    if (F_NRE*(nthreads+1) > idef->il_thread_division_nalloc)
    {
        idef->il_thread_division_nalloc = F_NRE*(nthreads+1);
        snew(idef->il_thread_division, idef->il_thread_division_nalloc);
    }

    ntype   = 0;
    nat_tot = 0;
    for (f = 0; f < F_NRE; f++)
    {
        if (!ftype_is_bonded_potential(f))
        {
            continue;
        }

        if (idef->il[f].nr == 0)
        {
            /* No interactions, avoid all the integer math below */
            for (t = 0; t <= nthreads; t++)
            {
                idef->il_thread_division[f*(nthreads+1)+t] = 0;
            }
        }
        else if (nthreads <= 4 || f == F_DISRES)
        {
            /* On up to 4 threads, load balancing the bonded work
             * is more important than minimizing the reduction cost.
             */
            int nat1, nr_t;

            nat1 = 1 + NRAL(f);

            for (t = 0; t <= nthreads; t++)
            {
                /* Divide equally over the threads */
                nr_t = (((idef->il[f].nr/nat1)*t)/nthreads)*nat1;

                if (f == F_DISRES)
                {
                    /* Ensure that distance restraint pairs with the same label
                     * end up on the same thread.
                     */
                    while (nr_t > 0 && nr_t < idef->il[f].nr &&
                           idef->iparams[idef->il[f].iatoms[nr_t]].disres.label ==
                           idef->iparams[idef->il[f].iatoms[nr_t-nat1]].disres.label)
                    {
                        nr_t += nat1;
                    }
                }

                idef->il_thread_division[f*(nthreads+1)+t] = nr_t;
            }
        }
        else
        {
            /* Add this ftype to the list to be distributed */
            ftype[ntype] = f;
            il[ntype]    = &idef->il[f];
            nat[ntype]   = NRAL(f);
            ind[ntype]   = 0;
            nat_tot     += il[ntype]->nr/(nat[ntype] + 1)*nat[ntype];

            /* The first index for the thread division is always 0 */
            idef->il_thread_division[ftype[ntype]*(nthreads+1)] = 0;

            ntype++;
        }
    }

    nat_sum = 0;
    for (t = 1; t <= nthreads; t++)
    {
        int nat_thread;

        /* Here we assume that the computational cost is proportional
         * to the number of atoms in the interaction. This is a rough
         * measure, but roughly correct. Usually there are very few
         * interactions anyhow and there are distributed relatively
         * uniformly. Proper and RB dihedrals are often distributed
         * non-uniformly, but their cost is roughly equal.
         */
        nat_thread = (nat_tot*t)/nthreads;

        while (nat_sum < nat_thread)
        {
            /* To divide bonds based on atom order, we compare
             * the index of the first atom in the bonded interaction.
             * This works well, since the domain decomposition generates
             * bondeds in order of the atoms by looking up interactions
             * which are linked to the first atom in each interaction.
             * It usually also works well without DD, since then the atoms
             * in bonded interactions are usually in increasing order.
             */
            int f_min, a_min = 0;

            f_min = -1;
            for (f = 0; f < ntype; f++)
            {
                if (ind[f] < il[f]->nr && (f_min == -1 ||
                                           il[f]->iatoms[ind[f]+1] < a_min))
                {
                    f_min = f;
                    a_min = il[f]->iatoms[ind[f]+1];
                }
            }
            assert(f_min >= 0);

            ind[f_min] += nat[f_min] + 1;
            nat_sum    += nat[f_min];
        }

        for (f = 0; f < ntype; f++)
        {
            idef->il_thread_division[ftype[f]*(nthreads+1) + t] = ind[f];
        }
    }

    for (f = 0; f < ntype; f++)
    {
        assert(ind[f] == il[f]->nr);
    }

}

static unsigned
calc_bonded_reduction_mask(const t_idef *idef,
                           int shift,
                           int t, int nt)
{
    unsigned mask;
    int      ftype, nb, nat1, nb0, nb1, i, a;

    mask = 0;

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
                        mask |= (1U << (idef->il[ftype].iatoms[i+a]>>shift));
                    }
                }
            }
        }
    }

    return mask;
}

void setup_bonded_threading(t_forcerec   *fr, t_idef *idef)
{
#define MAX_BLOCK_BITS 32
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

    /* We divide the force array in a maximum of 32 blocks.
     * Minimum force block reduction size is 2^6=64.
     */
    fr->red_ashift = 6;
    while (fr->natoms_force > (int)(MAX_BLOCK_BITS*(1U<<fr->red_ashift)))
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
        fr->f_t[t].red_mask =
            calc_bonded_reduction_mask(idef, fr->red_ashift, t, fr->nthreads);
    }

    /* Determine the maximum number of blocks we need to reduce over */
    fr->red_nblock = 0;
    ctot           = 0;
    for (t = 0; t < fr->nthreads; t++)
    {
        c = 0;
        for (b = 0; b < MAX_BLOCK_BITS; b++)
        {
            if (fr->f_t[t].red_mask & (1U<<b))
            {
                fr->red_nblock = std::max(fr->red_nblock, b+1);
                c++;
            }
        }
        if (debug)
        {
            fprintf(debug, "thread %d flags %x count %d\n",
                    t, fr->f_t[t].red_mask, c);
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
