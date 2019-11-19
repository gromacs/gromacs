//
// Created by sebkelle on 19.11.19.
//

#include <vector>

#include "gromacs/topology/block.h"
#include "gromacs/utility/smalloc.h"

#include "topology.h"

namespace nblib {

Topology TopologyBuilder::buildTopology(int numAtoms)
{
    numAtoms_ = numAtoms;
    return topology_;
}

void TopologyBuilder::setExclusions(std::vector<int> indices, std::vector<int> exclusions)
{
    snew(topology_.excls.index, numAtoms_ + 1);
    snew(topology_.excls.a, numAtoms_*numAtomsInMolecule);
    for (int index = 0; index < indices; index++)
    {
        topology_.excls[index].index = indices[index];
    }
    for (int exclusion = 0; exclusion < indices; exclusion++)
    {
        topology_.excls[exclusion].a = exlcusions[exclusion];
    }


    /*
    topology_.excls.index[0] = 0;
    for (int atom = 0; atom < numAtoms_; atom++)
    {
        const int firstAtomInMolecule = atom - (atom % numAtomsInMolecule);
        for (int atomJ = 0; atomJ < numAtomsInMolecule; atomJ++)
        {
            topology_.excls.a[atom*numAtomsInMolecule + atomJ] = firstAtomInMolecule + atomJ;
        }
        topology_.excls.index[atom + 1] = (atom + 1)*numAtomsInMolecule;
    }
    */
}

}
