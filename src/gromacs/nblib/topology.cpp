//
// Created by sebkelle on 19.11.19.
//

#include <vector>

#include "gromacs/topology/block.h"
#include "gromacs/utility/smalloc.h"

#include "topology.h"

namespace nblib {

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
