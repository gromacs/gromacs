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

/* struct for passing all data required for a function type */
typedef struct {
    int      ftype;
    t_ilist *il;
    int      nat;
    int      ind;
} ilist_data_t;

/* This routine attempts to divide all interactions of the ntype bondeds
 * types stored in ild over the threads such that each thread has roughly
 * equal load and different threads avoid touching the same atoms as much
 * as possible.
 */
static void divide_bondeds_by_locality(int           ntype,
                                       ilist_data_t *ild,
                                       int           nat_tot,
                                       int           nthread,
                                       t_idef       *idef)
{
    int nat_sum;
    int t;

    nat_sum = 0;
    /* Loop over the end bounds of the nthread threads to determine
     * which interactions threads 0 to nthread shall calculate.
     */
    for (t = 1; t <= nthread; t++)
    {
        int nat_thread;

        /* Here we assume that the computational cost is proportional
         * to the number of atoms in the interaction. This is a rough
         * measure, but roughly correct. Usually there are very few
         * interactions anyhow and there are distributed relatively
         * uniformly. Proper and RB dihedrals are often distributed
         * non-uniformly, but their cost is roughly equal.
         */
        nat_thread = (nat_tot*t)/nthread;

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
            int f_min;     /* the ild index with the smallest unassigned bonded atom index */
            int a_min = 0; /* the smallest unassigned bonded atom index */
            int f;

            /* Find out which of the ntype has the smallest atom index */
            f_min = -1;
            for (f = 0; f < ntype; f++)
            {
                if (ild[f].ind < ild[f].il->nr && (f_min == -1 ||
                                                   ild[f].il->iatoms[ild[f].ind + 1] < a_min))
                {
                    f_min = f;
                    a_min = ild[f].il->iatoms[ild[f].ind + 1];
                }
            }
            assert(f_min >= 0);

            /* Assign the interaction with the smallest atom index (of type
             * index f_min) to thread t-1 by increasing ind.
             */
            ild[f_min].ind += ild[f_min].nat + 1;
            nat_sum        += ild[f_min].nat;
        }

        /* Store the bonded end boundaries (at index t) for thread t-1 */
        int f;
        for (f = 0; f < ntype; f++)
        {
            idef->il_thread_division[ild[f].ftype*(nthread + 1) + t] = ild[f].ind;
        }
    }

    int f;
    for (f = 0; f < ntype; f++)
    {
        assert(ild[f].ind == ild[f].il->nr);
    }
}

static void divide_bondeds_over_threads(t_idef *idef, int nthread)
{
    /* The optimal value after which to switch from uniform to localized
     * bonded interaction distribution is 3, 4 or 5 depending on the system
     * and hardware.
     */
    const int    max_nthread_uniform = 4;
    ilist_data_t ild[F_NRE];
    int          nat_tot;
    int          ntype;
    int          f;

    assert(nthread > 0);

    idef->nthreads = nthread;

    if (F_NRE*(nthread + 1) > idef->il_thread_division_nalloc)
    {
        idef->il_thread_division_nalloc = F_NRE*(nthread + 1);
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
            int t;
            for (t = 0; t <= nthread; t++)
            {
                idef->il_thread_division[f*(nthread + 1) + t] = 0;
            }
        }
        else if (nthread <= max_nthread_uniform || f == F_DISRES)
        {
            /* On up to 4 threads, load balancing the bonded work
             * is more important than minimizing the reduction cost.
             */
            int nat1, t;

            /* nat1 = 1 + #atoms(ftype) which is the stride use for iatoms */
            nat1 = 1 + NRAL(f);

            for (t = 0; t <= nthread; t++)
            {
                int nr_t;

                /* Divide equally over the threads */
                nr_t = (((idef->il[f].nr/nat1)*t)/nthread)*nat1;

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

                idef->il_thread_division[f*(nthread + 1) + t] = nr_t;
            }
        }
        else
        {
            /* Add this ftype to the list to be distributed */
            int nat;

            nat              = NRAL(f);
            ild[ntype].ftype = f;
            ild[ntype].il    = &idef->il[f];
            ild[ntype].nat   = nat;
            ild[ntype].ind   = 0;
            nat_tot         += ild[ntype].il->nr/(nat + 1)*nat;

            /* The first index for the thread division is always 0 */
            idef->il_thread_division[f*(nthread + 1)] = 0;

            ntype++;
        }
    }

    if (ntype > 0)
    {
        divide_bondeds_by_locality(ntype, ild, nat_tot, nthread, idef);
    }

    if (debug)
    {
        int f;

        fprintf(debug, "Division of bondeds over threads:\n");
        for (f = 0; f < F_NRE; f++)
        {
            if (ftype_is_bonded_potential(f) && idef->il[f].nr > 0)
            {
                int t;

                fprintf(debug, "%16s", interaction_function[f].name);
                for (t = 0; t < nthread; t++)
                {
                    fprintf(debug, " %4d",
                            (idef->il_thread_division[f*(nthread + 1) + t + 1] -
                             idef->il_thread_division[f*(nthread + 1) + t])/
                            (1 + NRAL(f)));
                }
                fprintf(debug, "\n");
            }
        }
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
