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
    return topology_;
}

Topology& TopologyBuilder::add(MoleculeType, int nMolecules)
{

}


}
