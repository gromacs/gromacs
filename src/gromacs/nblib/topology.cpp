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

void TopologyBuilder::calculateTotalNumAtoms(std::vector<std::tuple<MoleculeType, int>> moleculesList)
{
    size_t numMolTypes = moleculesList.size();
    size_t numAtomsTotal = 0;
    for (size_t iMolTypes = 0; iMolTypes < numMolTypes; iMolTypes++)
    {
        int numAtomsInMolecule = std::get<0>(moleculesList[iMolTypes]).numAtomsInMolecule();
        int numMoleculesOfType = std::get<1>(moleculesList[iMolTypes]);
        numAtomsTotal += numAtomsInMolecule*numMoleculesOfType;
    }

    numAtoms_ = numAtomsTotal;
}

t_blocka TopologyBuilder::fillExclusionsList(std::vector<std::tuple<MoleculeType, int>> moleculesList)
{

    calculateTotalNumAtoms(moleculesList);
    //std::array<gmx::ExclusionBlock, numAtomsTotal> exclusionBlockGlobal;
    std::vector<gmx::ExclusionBlock> exclusionBlockGlobal(numAtoms_);

    t_blocka tBlockGlobal;
    size_t globalIndex = 0;

    for (size_t iMolTypes = 0; iMolTypes < moleculesList.size(); iMolTypes++)
    {
        MoleculeType& molecule = std::get<0>(moleculesList[iMolTypes]);
        size_t numMols = std::get<1>(moleculesList[iMolTypes]);

        int numAtomsInMolecule = molecule.numAtomsInMolecule();

        std::vector<std::tuple<int, int>> exclusionList = molecule.exclusions;
        std::vector<gmx::ExclusionBlock> exclusionsWithinMolecule(numAtomsInMolecule);

        for (int iAtoms = 0; iAtoms < numAtomsInMolecule; iAtoms++)
        {
            /*!
             * 1. From the exclusion list, extract all exclusions for the atom marked
             *    by iAtom
             * 2. Create an ExclusionBlock object and put into array of exclusion blocks (repeat this
             *    process numAtomsInMolecule)
             * 3. Duplicate this numMol times
             * 4. Concatonate this array with the global array
             * 5. Convert the array of exclusion blocks to a global t_blocka
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

                exclusionsWithinMolecule[iPair] = exclusionBlock;
            }
            for (size_t iMols = 0; iMols < numMols; iMols++ )
            {
                exclusionBlockGlobal[globalIndex] = exclusionsWithinMolecule[iAtoms];
                globalIndex++;
            }
        }
    }

    //! At the very end, outside of these 4 for loops, convert the exclusionBlockGlobal into
    //! a massive t_blocka and return

    gmx::exclusionBlocksToBlocka(exclusionBlockGlobal, &tBlockGlobal);

    return tBlockGlobal;
}

Topology TopologyBuilder::buildTopology()
{
    topology_.excls = fillExclusionsList(molecules_);
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

const std::vector<real>& Topology::getMasses() const
{
    return {};
}

const std::vector<real>& Topology::getCharges() const
{
    return {};
}

const std::vector<int>& Topology::getAtomTypes() const
{
    return {};
}

}
