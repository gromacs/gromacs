/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
/* \internal \file
 *
 * \brief Implements the UpdateGroupsCog class.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "updategroupscog.h"

#include <cstddef>

#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/updategroups.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxomp.h"

namespace gmx
{

UpdateGroupsCog::UpdateGroupsCog(const gmx_mtop_t& mtop,
                                 gmx::ArrayRef<const gmx::RangePartitioning> updateGroupingsPerMoleculeType,
                                 real temperature) :
    mtop_(mtop)
{
    int firstUpdateGroupInMolecule = 0;
    for (const auto& molblock : mtop.molblock)
    {
        const auto& updateGrouping = updateGroupingsPerMoleculeType[molblock.type];
        indicesPerMoleculeblock_.push_back({ firstUpdateGroupInMolecule, updateGrouping.numBlocks(), {} });
        auto& groupIndex = indicesPerMoleculeblock_.back().groupIndex_;

        for (int block = 0; block < updateGrouping.numBlocks(); block++)
        {
            groupIndex.insert(groupIndex.end(), updateGrouping.block(block).size(), block);
        }

        firstUpdateGroupInMolecule += molblock.nmol * updateGrouping.numBlocks();

        std::vector<bool> isFirstAtomInUpdateGroup(mtop.moltype[molblock.type].atoms.nr, false);
        for (int group = 0; group < updateGrouping.numBlocks(); group++)
        {
            isFirstAtomInUpdateGroup[*updateGrouping.block(group).begin()] = true;
        }
        isFirstAtomInUpdateGroup_.push_back(isFirstAtomInUpdateGroup);
    }

    maxUpdateGroupRadius_ = computeMaxUpdateGroupRadius(mtop, updateGroupingsPerMoleculeType, temperature);

    // Note that threadData_ is not allocated here as OpenMP might not be initialized yet
}

int UpdateGroupsCog::addCogsThread(ArrayRef<const int>  globalAtomIndices,
                                   ArrayRef<const RVec> coordinates,
                                   const int            cogBegin,
                                   const int            numThreads,
                                   const int            thread,
                                   const Range<int>&    threadAtomRange)
{
    // Prefetch the global update group indices, as we can OpenMP parallelize this
    int numUpdateGroups = 0;
    int moleculeBlock   = 0;
    for (int localAtom : threadAtomRange)
    {
        const int globalAtom = globalAtomIndices[localAtom];
        if (globalAtom < 0)
        {
            // This is a filler particle
            cogIndices_[localAtom] = -1;

            continue;
        }

        int moleculeIndex;
        int atomIndexInMolecule;
        mtopGetMolblockIndex(mtop_, globalAtom, &moleculeBlock, &moleculeIndex, &atomIndexInMolecule);
        const auto& indicesForBlock        = indicesPerMoleculeblock_[moleculeBlock];
        const int   globalUpdateGroupIndex = indicesForBlock.groupStart_
                                           + moleculeIndex * indicesForBlock.numGroupsPerMolecule_
                                           + indicesForBlock.groupIndex_[atomIndexInMolecule];

        // Temporarily store the global update group index in cogIndices_
        cogIndices_[localAtom] = globalUpdateGroupIndex;

        if (isFirstAtomInUpdateGroup_[moleculeBlock][atomIndexInMolecule])
        {
            numUpdateGroups++;
        }
    }

    threadData_[thread].numUpdateGroups_ = numUpdateGroups;

    if (numThreads > 1)
    {
        // OpenMP barrier so all threads can accumulate each others update group counts
#pragma omp barrier
    }

    int localCogIndexStart = cogBegin;
    for (int t = 0; t < thread; t++)
    {
        localCogIndexStart += threadData_[t].numUpdateGroups_;
    }

    auto& globalToLocalMap = threadData_[thread].globalToLocalMap_;

    int localCogIndex = localCogIndexStart;
    for (int localAtom : threadAtomRange)
    {
        const int globalUpdateGroupIndex = cogIndices_[localAtom];

        if (globalUpdateGroupIndex < 0)
        {
            continue;
        }

        if (const int* cogIndexPtr = globalToLocalMap.find(globalUpdateGroupIndex))
        {
            GMX_ASSERT(*cogIndexPtr >= cogBegin,
                       "Added atoms should not be part of previously present groups");

            cogIndices_[localAtom] = *cogIndexPtr;

            cogs_[*cogIndexPtr] += coordinates[localAtom];
            numAtomsPerCog_[*cogIndexPtr]++;
        }
        else
        {
            globalToLocalMap.insert(globalUpdateGroupIndex, localCogIndex);
            cogIndices_[localAtom]         = localCogIndex;
            cogs_[localCogIndex]           = coordinates[localAtom];
            numAtomsPerCog_[localCogIndex] = 1;

            localCogIndex++;
        }
    }

    GMX_ASSERT(
            localCogIndex - localCogIndexStart == numUpdateGroups,
            "We should have assigned as many local COGs as there are update groups on our thread");

    /* Divide sum of coordinates for each COG by the number of atoms */
    for (int i = localCogIndexStart; i < localCogIndex; i++)
    {
        const int numAtoms = numAtomsPerCog_[i];
        if (numAtoms > 1)
        {
            cogs_[i] /= numAtoms;
        }
    }

    return localCogIndex;
}

void UpdateGroupsCog::addCogs(ArrayRef<const int>  globalAtomIndices,
                              ArrayRef<const RVec> coordinates,
                              ArrayRef<const int>  numAtomsPerNbnxmGridColumn)
{
    const int numThreads = std::max(1, gmx_omp_nthreads_get(ModuleMultiThread::Domdec));
    if (threadData_.empty())
    {
        threadData_.resize(numThreads, { 0, HashedMap<int>(coordinates.size() / numThreads) });
    }
    else
    {
        GMX_RELEASE_ASSERT(gmx::ssize(threadData_) == numThreads,
                           "The number of threads should not change");
    }

    const int localAtomBegin = cogIndices_.size();
    const int localCogBegin  = cogs_.size();

    GMX_RELEASE_ASSERT(globalAtomIndices.ssize() >= localAtomBegin,
                       "addCogs should only be called to add COGs to the list that is already "
                       "present (which could be empty)");

    cogIndices_.resize(globalAtomIndices.size());
    // We overallocate the COG buffers with the number of atoms to avoid an extra OMP barrier
    cogs_.resize(globalAtomIndices.size());
    numAtomsPerCog_.resize(globalAtomIndices.size());

    int numCogs = -1;
    if (numAtomsPerNbnxmGridColumn.empty())
    {
        // When this list was passed as empty, we likely only have a few atoms to add,
        // do this single-threaded. Also we don't know how we can split the work.
        numCogs = addCogsThread(globalAtomIndices,
                                coordinates,
                                localCogBegin,
                                1,
                                0,
                                { localAtomBegin, int(globalAtomIndices.size()) });
    }
    else
    {
        const int totalNumAtomsToAdd = globalAtomIndices.size() - localAtomBegin;

        // Many atoms to process, use OpenMP threading
#pragma omp parallel num_threads(numThreads)
        {
            const int thread = gmx_omp_get_thread_num();

            // Divide the atom range over threads by assigning whole grid columns.
            // Grid columns can not broken up because then update groups can be broken up.
            // We try to assign roughly equal atoms count to the threads.
            int threadAtomBegin = 0;
            int threadAtomEnd   = 0;
            int numAtoms        = 0;
            for (int numAtomsInColumn : numAtomsPerNbnxmGridColumn)
            {
                numAtoms += numAtomsInColumn;
                if (numAtoms * numThreads <= totalNumAtomsToAdd * thread)
                {
                    threadAtomBegin = numAtoms;
                }
                if (numAtoms * numThreads <= totalNumAtomsToAdd * (thread + 1))
                {
                    threadAtomEnd = numAtoms;
                }
            }

            int numCogsThread =
                    addCogsThread(globalAtomIndices,
                                  coordinates,
                                  localCogBegin,
                                  numThreads,
                                  thread,
                                  { localAtomBegin + threadAtomBegin, localAtomBegin + threadAtomEnd });

            if (thread == numThreads - 1)
            {
                numCogs = numCogsThread;
            }
        }
    }

    GMX_ASSERT(numCogs >= 0, "numCogs should have been set in one of the two branches above");
    cogs_.resize(numCogs);
    numAtomsPerCog_.resize(numCogs);
}

void UpdateGroupsCog::clear()
{
    cogIndices_.clear();
    cogs_.clear();
    numAtomsPerCog_.clear();

    if (!threadData_.empty())
    {
#pragma omp parallel for schedule(static) num_threads(gmx::ssize(threadData_))
        for (int thread = 0; thread < gmx::ssize(threadData_); thread++)
        {
            threadData_[thread].globalToLocalMap_.clearAndResizeHashTable();
        }
    }
}

} // namespace gmx
