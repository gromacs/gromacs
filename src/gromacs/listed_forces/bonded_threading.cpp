/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief This file defines functions for managing threading of listed
 * interactions.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include "bonded_threading.h"

#include "config.h"

#include <climits>

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/bitmask.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"

#include "listed_internal.h"
#include "utilities.h"

namespace gmx
{

void BondedThreading::divideByLocality(const int numType,
                                       const EnumerationArray<InteractionFunction, IListData>& ild)
{
    int                                        nat_tot, nat_sum;
    EnumerationArray<InteractionFunction, int> ind; /* index into the ild[].il->iatoms */
    EnumerationArray<InteractionFunction, int> at_ind; /* index of the first atom of the interaction at ind */
    int f, t;

    GMX_ASSERT(numType <= static_cast<int>(InteractionFunction::Count), "Too many function types");

    nat_tot = 0;
    for (f = 0; f < numType; f++)
    {
        /* Sum #bondeds*#atoms_per_bond over all bonded types */
        nat_tot += ild[f].il->size() / (ild[f].numAtoms + 1) * ild[f].numAtoms;
        /* The start bound for thread 0 is 0 for all interactions */
        ind[f] = 0;
        /* Initialize the next atom index array */
        assert(!ild[f].il->empty());
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
    for (t = 1; t <= numThreads_; t++)
    {
        int nat_thread;

        /* Here we assume that the computational cost is proportional
         * to the number of atoms in the interaction. This is a rough
         * measure, but roughly correct. Usually there are very few
         * interactions anyhow and there are distributed relatively
         * uniformly. Proper and RB dihedrals are often distributed
         * non-uniformly, but their cost is roughly equal.
         */
        nat_thread = (nat_tot * t) / numThreads_;

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
            GMX_ASSERT(f_min >= 0 && f_min < numType, "Function type out of range");

            /* Assign the interaction with the lowest atom index (of type
             * index f_min) to thread t-1 by increasing ind.
             */
            ind[f_min] += ild[f_min].numAtoms + 1;
            nat_sum += ild[f_min].numAtoms;

            /* Update the first unassigned atom index for this type */
            if (ind[f_min] < ild[f_min].il->size())
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
            workDivision_.setBound(ild[f].ftype, t, ind[f]);
        }
    }

    for (f = 0; f < numType; f++)
    {
        GMX_ASSERT(ind[f] == ild[f].il->size(), "Index is consistent");
    }
}

//! Return whether function type \p ftype in \p idef has perturbed interactions
static bool ftypeHasPerturbedEntries(const InteractionDefinitions& idef, InteractionFunction ftype)
{
    GMX_ASSERT(idef.ilsort == ilsortNO_FE || idef.ilsort == ilsortFE_SORTED,
               "Perturbed interactions should be sorted here");

    const InteractionList& ilist = idef.il[ftype];

    return (idef.ilsort != ilsortNO_FE && idef.numNonperturbedInteractions[ftype] != ilist.size());
}

void BondedThreading::divideWork(const bool useGpuForBondeds, const InteractionDefinitions& idef)
{
    EnumerationArray<InteractionFunction, IListData> ild;

    GMX_ASSERT(numThreads_ > 0, "Must have positive number of threads");
    const int numThreads = numThreads_;

    ArrayRef<const t_iparams> iparams = idef.iparams;

    haveBondeds_         = false;
    int    numType       = 0;
    size_t fTypeGpuIndex = 0;
    for (const auto fType : EnumerationWrapper<InteractionFunction>{})
    {
        if (!ftypeIsListedPotential(fType))
        {
            continue;
        }

        const InteractionList& il                     = idef.il[fType];
        int                    nrToAssignToCpuThreads = il.size();

        if (useGpuForBondeds && fTypeGpuIndex < fTypesOnGpu.size() && fTypesOnGpu[fTypeGpuIndex] == fType)
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
            haveBondeds_ = true;
        }

        if (nrToAssignToCpuThreads == 0)
        {
            /* No interactions, avoid all the integer math below */
            for (int t = 0; t <= numThreads; t++)
            {
                workDivision_.setBound(fType, t, 0);
            }
        }
        else if (numThreads <= maxNthreadUniform_ || fType == InteractionFunction::DistanceRestraints)
        {
            /* On up to 4 threads, load balancing the bonded work
             * is more important than minimizing the reduction cost.
             */

            const int stride = 1 + NRAL(fType);

            for (int t = 0; t <= numThreads; t++)
            {
                /* Divide equally over the threads */
                int nr_t = (((nrToAssignToCpuThreads / stride) * t) / numThreads) * stride;

                if (fType == InteractionFunction::DistanceRestraints)
                {
                    /* Ensure that distance restraint pairs with the same label
                     * end up on the same thread.
                     */
                    while (nr_t > 0 && nr_t < nrToAssignToCpuThreads
                           && iparams[il.iatoms[nr_t]].disres.label
                                      == iparams[il.iatoms[nr_t - stride]].disres.label)
                    {
                        nr_t += stride;
                    }
                }

                workDivision_.setBound(fType, t, nr_t);
            }
        }
        else
        {
            /* Add this fType to the list to be distributed */
            int numAtoms          = NRAL(fType);
            ild[numType].ftype    = fType;
            ild[numType].il       = &il;
            ild[numType].numAtoms = numAtoms;

            /* The first index for the thread division is always 0 */
            workDivision_.setBound(fType, 0, 0);

            numType++;
        }
    }

    if (numType > 0)
    {
        divideByLocality(numType, ild);
    }

    if (debug)
    {
        fprintf(debug, "Division of bondeds over threads:\n");
        for (const auto f : EnumerationWrapper<InteractionFunction>{})
        {
            if (ftypeIsListedPotential(f) && !idef.il[f].empty())
            {
                int t;

                fprintf(debug, "%16s", interaction_function[f].name);
                for (t = 0; t < numThreads; t++)
                {
                    fprintf(debug,
                            " %4d",
                            (workDivision_.bound(f, t + 1) - workDivision_.bound(f, t)) / (1 + NRAL(f)));
                }
                fprintf(debug, "\n");
            }
        }
    }
}

void BondedThreading::calculateReductionMask(const int natoms, const InteractionDefinitions& idef, const int thread)
{
    static_assert(BITMASK_SIZE == GMX_OPENMP_MAX_THREADS,
                  "For the error message below we assume these two are equal.");

    if (numThreads_ > BITMASK_SIZE)
    {
        gmx_fatal(FARGS,
                  "You are using %d OpenMP threads, which is larger than GMX_OPENMP_MAX_THREADS "
                  "(%d). Decrease the number of OpenMP threads or rebuild GROMACS with a larger "
                  "value for GMX_OPENMP_MAX_THREADS passed to CMake.",
                  numThreads_,
                  GMX_OPENMP_MAX_THREADS);
    }
    GMX_ASSERT(numThreads_ <= BITMASK_SIZE, "We need at least nthreads bits in the mask");

    ThreadForceBuffer<rvec4>& forceBufferForThisThread = threadedForceBuffers_.threadForceBuffer(thread);
    forceBufferForThisThread.resizeBufferAndClearMask(natoms);

    for (const auto ftype : EnumerationWrapper<InteractionFunction>{})
    {
        if (ftypeIsListedPotential(ftype))
        {
            int nb = idef.il[ftype].size();
            if (nb > 0)
            {
                int nat1 = interaction_function[ftype].nratoms + 1;

                int nb0 = workDivision_.bound(ftype, thread);
                int nb1 = workDivision_.bound(ftype, thread + 1);

                for (int i = nb0; i < nb1; i += nat1)
                {
                    for (int a = 1; a < nat1; a++)
                    {
                        forceBufferForThisThread.addAtomToMask(idef.il[ftype].iatoms[i + a]);
                    }
                }
            }
        }
    }

    forceBufferForThisThread.processMask();
}

void BondedThreading::updateAfterPartition(const int                     numAtomsForce,
                                           const bool                    useGpuForBondeds,
                                           const InteractionDefinitions& idef)
{
    GMX_ASSERT(numThreads_ >= 1, "Must have at least one thread");

    divideWork(useGpuForBondeds, idef);

    if (!haveBondeds_)
    {
        /* We don't have bondeds, so there is nothing to reduce */
        return;
    }

    /* Determine to which blocks each thread's bonded force calculation
     * contributes. Store this as a mask for each thread.
     */
#pragma omp parallel for num_threads(numThreads_) schedule(static)
    for (int t = 0; t < numThreads_; t++)
    {
        try
        {
            // TODO eliminate this arg, t is already a param
            calculateReductionMask(numAtomsForce, idef, t);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }

    threadedForceBuffers_.setupReduction();
}

BondedThreading::BondedThreading(const int       numThreads,
                                 const int       numEnergyGroups,
                                 int             numComGroups,
                                 const MDLogger& mdlog) :
    numThreads_(numThreads),
    threadedForceBuffers_(numThreads, true, numEnergyGroups),
    haveBondeds_(false),
    workDivision_(numThreads_),
    foreignLambdaWorkDivision_(1)
{
    /* Note that we also use threadedForceBuffer when running single-threaded.
     * This is because the bonded force buffer uses type rvec4, whereas
     * the normal force buffer is uses type rvec. This leads to a little
     * reduction overhead, but the speed gain in the bonded calculations
     * of doing transposeScatterIncr/DecrU with aligment 4 instead of 3
     * is much larger than the reduction overhead.
     */

    /* The optimal value after which to switch from uniform to localized
     * bonded interaction distribution is 3, 4 or 5 depending on the system
     * and hardware.
     */
    const int maxNthreadUniformDefault = 4;
    char*     ptr;

    if ((ptr = std::getenv("GMX_BONDED_NTHREAD_UNIFORM")) != nullptr)
    {
        sscanf(ptr, "%d", &maxNthreadUniform_);
        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        "Max threads for uniform bonded distribution set to %d by env.var.\n",
                        maxNthreadUniform_);
    }
    else
    {
        maxNthreadUniform_ = maxNthreadUniformDefault;
    }

    centersOfMassScaledBuffers_.resize(numThreads);
#pragma omp parallel for schedule(static) num_threads(numThreads_)
    for (int thread = 0; thread < numThreads_; thread++)
    {
        auto& bufs = centersOfMassScaledBuffers_[thread];
        bufs.comA_.resize(std::max(numComGroups, 1), { 0.0, 0.0, 0.0 });
        bufs.comB_.resize(std::max(numComGroups, 1), { 0.0, 0.0, 0.0 });
    }
}

} // namespace gmx
