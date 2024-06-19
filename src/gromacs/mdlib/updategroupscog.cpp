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

#include "gromacs/mdlib/updategroups.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

UpdateGroupsCog::UpdateGroupsCog(const gmx_mtop_t&                           mtop,
                                 gmx::ArrayRef<const gmx::RangePartitioning> updateGroupingsPerMoleculeType,
                                 real                                        temperature,
                                 int                                         numHomeAtoms) :
    globalToLocalMap_(numHomeAtoms), mtop_(mtop)
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
    }

    maxUpdateGroupRadius_ = computeMaxUpdateGroupRadius(mtop, updateGroupingsPerMoleculeType, temperature);
}

void UpdateGroupsCog::addCogs(gmx::ArrayRef<const int>       globalAtomIndices,
                              gmx::ArrayRef<const gmx::RVec> coordinates)
{
    const int    localAtomBegin = cogIndices_.size();
    const size_t cogBegin       = cogs_.size();

    GMX_RELEASE_ASSERT(globalAtomIndices.ssize() >= localAtomBegin,
                       "addCogs should only be called to add COGs to the list that is already "
                       "present (which could be empty)");

    cogIndices_.reserve(globalAtomIndices.size());

    int moleculeBlock = 0;
    for (gmx::Index localAtom = localAtomBegin; localAtom < globalAtomIndices.ssize(); localAtom++)
    {
        const int globalAtom = globalAtomIndices[localAtom];
        int       moleculeIndex;
        int       atomIndexInMolecule;
        mtopGetMolblockIndex(mtop_, globalAtom, &moleculeBlock, &moleculeIndex, &atomIndexInMolecule);
        const auto& indicesForBlock        = indicesPerMoleculeblock_[moleculeBlock];
        int         globalUpdateGroupIndex = indicesForBlock.groupStart_
                                     + moleculeIndex * indicesForBlock.numGroupsPerMolecule_
                                     + indicesForBlock.groupIndex_[atomIndexInMolecule];

        if (const int* cogIndexPtr = globalToLocalMap_.find(globalUpdateGroupIndex))
        {
            GMX_ASSERT(static_cast<size_t>(*cogIndexPtr) >= cogBegin,
                       "Added atoms should not be part of previously present groups");

            cogIndices_.push_back(*cogIndexPtr);

            cogs_[*cogIndexPtr] += coordinates[localAtom];
            numAtomsPerCog_[*cogIndexPtr]++;
        }
        else
        {
            const int cogIndex = cogs_.size();

            globalToLocalMap_.insert(globalUpdateGroupIndex, cogIndex);
            cogIndices_.push_back(cogIndex);
            cogs_.push_back(coordinates[localAtom]);
            numAtomsPerCog_.push_back(1);
        }
    }

    /* Divide sum of coordinates for each COG by the number of atoms */
    for (size_t i = cogBegin; i < cogs_.size(); i++)
    {
        const int numAtoms = numAtomsPerCog_[i];
        if (numAtoms > 1)
        {
            cogs_[i] /= numAtoms;
        }
    }
}

void UpdateGroupsCog::clear()
{
    cogIndices_.clear();
    cogs_.clear();
    numAtomsPerCog_.clear();
    globalToLocalMap_.clearAndResizeHashTable();
}

} // namespace gmx
