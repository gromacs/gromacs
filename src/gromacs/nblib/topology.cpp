//
// Created by sebkelle on 19.11.19.
//

#include <vector>

#include "gromacs/topology/block.h"

#include "topology.h"

/*
Topology buildTopology(int numAtoms)
{
    numAtoms_ = numAtoms;
}
TopologyBuilder::setExclusions(std::vector<int> indices, std::vector<int> exclusions, int numAtomsInMolecule)
{
    snew(topology_.excls.index, numAtoms_ + 1);
    snew(topology_.excls.a, numAtoms_*numAtomsInMolecule);
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

}
*/
