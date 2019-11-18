/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include "manage_threading.h"

#include "config.h"

#include <cassert>
#include <cinttypes>
#include <climits>
#include <cstdlib>

#include <algorithm>
#include <string>

#include "gromacs/listed_forces/gpubonded.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "listed_internal.h"
#include "utilities.h"

/*! \brief struct for passing all data required for a function type */
typedef struct
{
    const t_ilist* il;    /**< pointer to t_ilist entry corresponding to ftype */
    int            ftype; /**< the function type index */
    int            nat;   /**< nr of atoms involved in a single ftype interaction */
} ilist_data_t;

/*! \brief Divides listed interactions over threads
 *
 * This routine attempts to divide all interactions of the numType bondeds
 * types stored in ild over the threads such that each thread has roughly
 * equal load and different threads avoid touching the same atoms as much
 * as possible.
 */
static void divide_bondeds_by_locality(bonded_threading_t* bt, int numType, const ilist_data_t* ild)
{
    int nat_tot, nat_sum;
    int ind[F_NRE];    /* index into the ild[].il->iatoms */
    int at_ind[F_NRE]; /* index of the first atom of the interaction at ind */
    int f, t;

    assert(numType <= F_NRE);

    nat_tot = 0;
    for (f = 0; f < numType; f++)
    {
        /* Sum #bondeds*#atoms_per_bond over all bonded types */
        nat_tot += ild[f].il->nr / (ild[f].nat + 1) * ild[f].nat;
        /* The start bound for thread 0 is 0 for all interactions */
        ind[f] = 0;
        /* Initialize the next atom index array */
        assert(ild[f].il->nr > 0);
        at_ind[f] = ild[f].il->iatoms[1];
    }

    nat_sum = 0;
    /* Loop over the end bounds of the nthreads threads to determine
     * which interactions threads 0 to nthreads shall calculate.
     *
     * NOTE: The cost of these combined loops is #interactions*numType.
     * This code is running single threaded (difficult to parallelize
     * over threads). So the relative cost of this function increases
     * linearly with the number of threads. Since the inner-most loop
     * is cheap and this is done only at DD repartitioning, the cost should
     * be negligble. At high thread count many other parts of the code
     * scale the same way, so it's (currently) not worth improving this.
     */
    for (t = 1; t <= bt->nthreads; t++)
    {
        int nat_thread;

        /* Here we assume that the computational cost is proportional
         * to the number of atoms in the interaction. This is a rough
         * measure, but roughly correct. Usually there are very few
         * interactions anyhow and there are distributed relatively
         * uniformly. Proper and RB dihedrals are often distributed
         * non-uniformly, but their cost is roughly equal.
         */
        nat_thread = (nat_tot * t) / bt->nthreads;

        while (nat_sum < nat_thread)
        {
            /* To divide bonds based on atom order, we compare
             * the index of the first atom in the bonded interaction.
             * This works well, since the domain decomposition generates
             * bondeds in order of the atoms by looking up interactions
             * which are linked to the first atom in each interaction.
             * It usually also works well without DD, since than the atoms
             * in bonded interactions are usually in increasing order.
             * If they are not assigned in increasing order, the balancing
             * is still good, but the memory access and reduction cost will
             * be higher.
             */
            int f_min;

            /* Find out which of the types has the lowest atom index */
            f_min = 0;
            for (f = 1; f < numType; f++)
            {
                if (at_ind[f] < at_ind[f_min])
                {
                    f_min = f;
                }
            }
            assert(f_min >= 0 && f_min < numType);

            /* Assign the interaction with the lowest atom index (of type
             * index f_min) to thread t-1 by increasing ind.
             */
            ind[f_min] += ild[f_min].nat + 1;
            nat_sum += ild[f_min].nat;

            /* Update the first unassigned atom index for this type */
            if (ind[f_min] < ild[f_min].il->nr)
            {
                at_ind[f_min] = ild[f_min].il->iatoms[ind[f_min] + 1];
            }
            else
            {
                /* We have assigned all interactions of this type.
                 * Setting at_ind to INT_MAX ensures this type will not be
                 * chosen in the for loop above during next iterations.
                 */
                at_ind[f_min] = INT_MAX;
            }
        }

        /* Store the bonded end boundaries (at index t) for thread t-1 */
        for (f = 0; f < numType; f++)
        {
            bt->workDivision.setBound(ild[f].ftype, t, ind[f]);
        }
    }

    for (f = 0; f < numType; f++)
    {
        assert(ind[f] == ild[f].il->nr);
    }
}

//! Return whether function type \p ftype in \p idef has perturbed interactions
static bool ftypeHasPerturbedEntries(const t_idef& idef, int ftype)
{
    GMX_ASSERT(idef.ilsort == ilsortNO_FE || idef.ilsort == ilsortFE_SORTED,
               "Perturbed interations should be sorted here");

    const t_ilist& ilist = idef.il[ftype];

    return (idef.ilsort != ilsortNO_FE && ilist.nr_nonperturbed != ilist.nr);
}

//! Divides bonded interactions over threads and GPU
static void divide_bondeds_over_threads(bonded_threading_t* bt, bool useGpuForBondeds, const t_idef& idef)
{
    ilist_data_t ild[F_NRE];

    GMX_ASSERT(bt->nthreads > 0, "Must have positive number of threads");
    const int numThreads = bt->nthreads;

    bt->haveBondeds      = false;
    int    numType       = 0;
    size_t fTypeGpuIndex = 0;
    for (int fType = 0; fType < F_NRE; fType++)
    {
        if (!ftype_is_bonded_potential(fType))
        {
            continue;
        }

        const t_ilist& il                     = idef.il[fType];
        int            nrToAssignToCpuThreads = il.nr;

        if (useGpuForBondeds && fTypeGpuIndex < gmx::fTypesOnGpu.size()
            && gmx::fTypesOnGpu[fTypeGpuIndex] == fType)
        {
            fTypeGpuIndex++;

            /* Perturbation is not implemented in the GPU bonded kernels.
             * But instead of doing all on the CPU, we could do only
             * the actually perturbed interactions on the CPU.
             */
            if (!ftypeHasPerturbedEntries(idef, fType))
            {
                /* We will assign this interaction type to the GPU */
                nrToAssignToCpuThreads = 0;
            }
        }

        if (nrToAssignToCpuThreads > 0)
        {
            bt->haveBondeds = true;
        }

        if (nrToAssignToCpuThreads == 0)
        {
            /* No interactions, avoid all the integer math below */
            for (int t = 0; t <= numThreads; t++)
            {
                bt->workDivision.setBound(fType, t, 0);
            }
        }
        else if (numThreads <= bt->max_nthread_uniform || fType == F_DISRES)
        {
            /* On up to 4 threads, load balancing the bonded work
             * is more important than minimizing the reduction cost.
             */

            const int stride = 1 + NRAL(fType);

            for (int t = 0; t <= numThreads; t++)
            {
                /* Divide equally over the threads */
                int nr_t = (((nrToAssignToCpuThreads / stride) * t) / numThreads) * stride;

                if (fType == F_DISRES)
                {
                    /* Ensure that distance restraint pairs with the same label
                     * end up on the same thread.
                     */
                    while (nr_t > 0 && nr_t < nrToAssignToCpuThreads
                           && idef.iparams[il.iatoms[nr_t]].disres.label
                                      == idef.iparams[il.iatoms[nr_t - stride]].disres.label)
                    {
                        nr_t += stride;
                    }
                }

                bt->workDivision.setBound(fType, t, nr_t);
            }
        }
        else
        {
            /* Add this fType to the list to be distributed */
            int nat            = NRAL(fType);
            ild[numType].ftype = fType;
            ild[numType].il    = &il;
            ild[numType].nat   = nat;

            /* The first index for the thread division is always 0 */
            bt->workDivision.setBound(fType, 0, 0);

            numType++;
        }
    }

    if (numType > 0)
    {
        divide_bondeds_by_locality(bt, numType, ild);
    }

    if (debug)
    {
        int f;

        fprintf(debug, "Division of bondeds over threads:\n");
        for (f = 0; f < F_NRE; f++)
        {
            if (ftype_is_bonded_potential(f) && idef.il[f].nr > 0)
            {
                int t;

                fprintf(debug, "%16s", interaction_function[f].name);
                for (t = 0; t < numThreads; t++)
                {
                    fprintf(debug, " %4d",
                            (bt->workDivision.bound(f, t + 1) - bt->workDivision.bound(f, t))
                                    / (1 + NRAL(f)));
                }
                fprintf(debug, "\n");
            }
        }
    }
}

//! Construct a reduction mask for which parts (blocks) of the force array are touched on which thread task
static void calc_bonded_reduction_mask(int                       natoms,
                                       f_thread_t*               f_thread,
                                       const t_idef&             idef,
                                       int                       thread,
                                       const bonded_threading_t& bondedThreading)
{
    static_assert(BITMASK_SIZE == GMX_OPENMP_MAX_THREADS,
                  "For the error message below we assume these two are equal.");

    if (bondedThreading.nthreads > BITMASK_SIZE)
    {
#pragma omp master
        gmx_fatal(FARGS,
                  "You are using %d OpenMP threads, which is larger than GMX_OPENMP_MAX_THREADS "
                  "(%d). Decrease the number of OpenMP threads or rebuild GROMACS with a larger "
                  "value for GMX_OPENMP_MAX_THREADS passed to CMake.",
                  bondedThreading.nthreads, GMX_OPENMP_MAX_THREADS);
#pragma omp barrier
    }
    GMX_ASSERT(bondedThreading.nthreads <= BITMASK_SIZE,
               "We need at least nthreads bits in the mask");

    int nblock = (natoms + reduction_block_size - 1) >> reduction_block_bits;

    if (nblock > f_thread->block_nalloc)
    {
        f_thread->block_nalloc = over_alloc_large(nblock);
        srenew(f_thread->mask, f_thread->block_nalloc);
        srenew(f_thread->block_index, f_thread->block_nalloc);
        // NOTE: It seems f_thread->f does not need to be aligned
        sfree_aligned(f_thread->f);
        snew_aligned(f_thread->f, f_thread->block_nalloc * reduction_block_size, 128);
    }

    gmx_bitmask_t* mask = f_thread->mask;

    for (int b = 0; b < nblock; b++)
    {
        bitmask_clear(&mask[b]);
    }

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (ftype_is_bonded_potential(ftype))
        {
            int nb = idef.il[ftype].nr;
            if (nb > 0)
            {
                int nat1 = interaction_function[ftype].nratoms + 1;

                int nb0 = bondedThreading.workDivision.bound(ftype, thread);
                int nb1 = bondedThreading.workDivision.bound(ftype, thread + 1);

                for (int i = nb0; i < nb1; i += nat1)
                {
                    for (int a = 1; a < nat1; a++)
                    {
                        bitmask_set_bit(&mask[idef.il[ftype].iatoms[i + a] >> reduction_block_bits], thread);
                    }
                }
            }
        }
    }

    /* Make an index of the blocks our thread touches, so we can do fast
     * force buffer clearing.
     */
    f_thread->nblock_used = 0;
    for (int b = 0; b < nblock; b++)
    {
        if (bitmask_is_set(mask[b], thread))
        {
            f_thread->block_index[f_thread->nblock_used++] = b;
        }
    }
}

void setup_bonded_threading(bonded_threading_t* bt, int numAtoms, bool useGpuForBondeds, const t_idef& idef)
{
    int ctot = 0;

    assert(bt->nthreads >= 1);

    /* Divide the bonded interaction over the threads */
    divide_bondeds_over_threads(bt, useGpuForBondeds, idef);

    if (!bt->haveBondeds)
    {
        /* We don't have bondeds, so there is nothing to reduce */
        return;
    }

    /* Determine to which blocks each thread's bonded force calculation
     * contributes. Store this as a mask for each thread.
     */
#pragma omp parallel for num_threads(bt->nthreads) schedule(static)
    for (int t = 0; t < bt->nthreads; t++)
    {
        try
        {
            calc_bonded_reduction_mask(numAtoms, bt->f_t[t].get(), idef, t, *bt);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }

    /* Reduce the masks over the threads and determine which blocks
     * we need to reduce over.
     */
    int nblock_tot = (numAtoms + reduction_block_size - 1) >> reduction_block_bits;
    /* Ensure we have sufficient space for all blocks */
    if (static_cast<size_t>(nblock_tot) > bt->block_index.size())
    {
        bt->block_index.resize(nblock_tot);
    }
    if (static_cast<size_t>(nblock_tot) > bt->mask.size())
    {
        bt->mask.resize(nblock_tot);
    }
    bt->nblock_used = 0;
    for (int b = 0; b < nblock_tot; b++)
    {
        gmx_bitmask_t* mask = &bt->mask[b];

        /* Generate the union over the threads of the bitmask */
        bitmask_clear(mask);
        for (int t = 0; t < bt->nthreads; t++)
        {
            bitmask_union(mask, bt->f_t[t]->mask[b]);
        }
        if (!bitmask_is_zero(*mask))
        {
            bt->block_index[bt->nblock_used++] = b;
        }

        if (debug)
        {
            int c = 0;
            for (int t = 0; t < bt->nthreads; t++)
            {
                if (bitmask_is_set(*mask, t))
                {
                    c++;
                }
            }
            ctot += c;

            if (gmx_debug_at)
            {
                fprintf(debug, "block %d flags %s count %d\n", b, to_hex_string(*mask).c_str(), c);
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "Number of %d atom blocks to reduce: %d\n", reduction_block_size, bt->nblock_used);
        fprintf(debug, "Reduction density %.2f for touched blocks only %.2f\n",
                ctot * reduction_block_size / static_cast<double>(numAtoms),
                ctot / static_cast<double>(bt->nblock_used));
    }
}

void tear_down_bonded_threading(bonded_threading_t* bt)
{
    delete bt;
}

f_thread_t::f_thread_t(int numEnergyGroups) : grpp(numEnergyGroups)
{
    snew(fshift, SHIFTS);
}

f_thread_t::~f_thread_t()
{
    sfree(mask);
    sfree(fshift);
    sfree(block_index);
    sfree_aligned(f);
}

bonded_threading_t::bonded_threading_t(const int numThreads, const int numEnergyGroups) :
    nthreads(numThreads),
    nblock_used(0),
    haveBondeds(false),
    workDivision(nthreads),
    foreignLambdaWorkDivision(1)
{
    f_t.resize(numThreads);
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int t = 0; t < nthreads; t++)
    {
        try
        {
            /* Note that thread 0 uses the global fshift and energy arrays,
             * but to keep the code simple, we initialize all data here.
             */
            f_t[t] = std::make_unique<f_thread_t>(numEnergyGroups);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
}

bonded_threading_t* init_bonded_threading(FILE* fplog, const int nenergrp)
{
    /* These thread local data structures are used for bondeds only.
     *
     * Note that we also use there structures when running single-threaded.
     * This is because the bonded force buffer uses type rvec4, whereas
     * the normal force buffer is uses type rvec. This leads to a little
     * reduction overhead, but the speed gain in the bonded calculations
     * of doing transposeScatterIncr/DecrU with aligment 4 instead of 3
     * is much larger than the reduction overhead.
     */
    bonded_threading_t* bt = new bonded_threading_t(gmx_omp_nthreads_get(emntBonded), nenergrp);

    /* The optimal value after which to switch from uniform to localized
     * bonded interaction distribution is 3, 4 or 5 depending on the system
     * and hardware.
     */
    const int max_nthread_uniform = 4;
    char*     ptr;

    if ((ptr = getenv("GMX_BONDED_NTHREAD_UNIFORM")) != nullptr)
    {
        sscanf(ptr, "%d", &bt->max_nthread_uniform);
        if (fplog != nullptr)
        {
            fprintf(fplog, "\nMax threads for uniform bonded distribution set to %d by env.var.\n",
                    bt->max_nthread_uniform);
        }
    }
    else
    {
        bt->max_nthread_uniform = max_nthread_uniform;
    }

    return bt;
}
