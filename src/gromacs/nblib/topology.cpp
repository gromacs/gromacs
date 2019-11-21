//
// Created by sebkelle on 19.11.19.
//

#include <vector>
#include <array>

#include "gromacs/topology/block.h"
#include "gromacs/utility/smalloc.h"

#include "topology.h"

#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/topology/block.h"

namespace nblib {

t_blocka TopologyBuilder::fillExclusionsList(std::vector<std::tuple<MoleculeType, int>> moleculesList)
{
    /*!
     * 1. Create an empty global exclusions list
     * 2. Iterate through the molecules list and grab its exclusions list
     * 3. Convert that exclusions list into a gmx::ExclusionBlock
     * 4. Merge it with the number of molecules of that type
     * 5. Do the same for all molecule types and return the big ExclusionBlock
     */
    size_t numMolTypes = moleculesList.size();
    size_t numAtomsTotal = 0;
    for (size_t iMolTypes = 0; iMolTypes < numMolTypes; iMolTypes++)
    {
        int numAtomsInMolecule = std::get<0>(moleculesList[iMolTypes]).numAtomsInMolecule();
        numAtomsTotal += numAtomsInMolecule;
    }

    std::array<gmx::ExclusionBlock, numAtomsTotal> exclusionBlockGlobal;

    t_blocka tBlockGlobal;

    for (size_t iMolTypes = 0; iMolTypes < moleculesList.size(); iMolTypes++)
    {
        int numMols = std::get<1>(moleculesList[iMolTypes]);

        int numAtomsInMolecule = std::get<0>(moleculesList[iMolTypes]).numAtomsInMolecule();

        std::vector<std::tuple<int, int>> exclusionList = std::get<0>(moleculesList[iMolTypes]).exclusions;

        for (int iAtoms = 0; iAtoms < numAtomsInMolecule; iAtoms++)
        {
            /*!
             * 1. From the exclusion list, extract all exclusions for the atom marked
             *    by iAtom
             * 2. Create an ExclusionBlock object and put into array of exclusion blocks (repeat this
             *    process numAtomsInMolecule times
             * 3. Convert the array of exclusion blocks to a global t_blocka
             */

            for (size_t iPair = 0; iPair < exclusionList.size(); iPair++)
            {
                std::vector<int> atomsToExclude;
                if (std::get<0>(exclusionList[iPair]) == iAtoms)
                {
                    atomsToExclude.emplace_back(std::get<1>(exclusionList[iPair]));
                }
                gmx::ExclusionBlock exclusionBlock;
                exclusionBlock.atomNumber = atomsToExclude;


            }

        }
    }
    //! At the very end, outside of these 4 for loops, convert the global exclusion block into
    //! a massive t_blocka and return

    return tBlockGlobal;
}

Topology TopologyBuilder::buildTopology()
{
    return topology_;
}

TopologyBuilder& TopologyBuilder::add(const MoleculeType moleculeType, const int nMolecules)
{
    /*!
     * 1. Push-back a tuple of molecule type and nMolecules
     * 2. Append exclusion list into the data structure
     */

    molecules_.emplace_back(std::make_tuple(moleculeType, nMolecules));
    numAtoms_ += nMolecules*moleculeType.numAtomsInMolecule();
}


}
