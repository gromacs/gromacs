//
// Created by sebkelle on 19.11.19.
//

#include <vector>

#include "gromacs/math/vec.h"

#include "molecules.h"

#ifndef GROMACS_TOPOLOGY_H
#define GROMACS_TOPOLOGY_H

struct t_blocka;

namespace nblib {

class Topology {
public:

    const std::vector<int>& getAtomTypes() const;

    const std::vector<real>& getCharges() const;

    const std::vector<real>& getMasses() const;

private:
    Topology();

    friend class TopologyBuilder;

    //! Storage for parameters for short range interactions.
    std::vector<real>      nonbondedParameters;
    //! Storage for atom type parameters.
    std::vector<int>       atomTypes;
    //! Storage for atom partial charges.
    std::vector<real>      charges;
    //! Atom masses
    std::vector<real>      masses;
    //! Atom info where all atoms are marked to have Van der Waals interactions
    std::vector<int>       atomInfoAllVdw;
    //! Information about exclusions.
    t_blocka               excls;
};

class TopologyBuilder {
public:
    TopologyBuilder();

    Topology buildTopology();

    TopologyBuilder& add(MoleculeType moleculeType, int nMolecules);

private:
    Topology topology_;

    int numAtoms_;
    std::vector<std::tuple<MoleculeType, int>> molecules_;

    std::vector<std::tuple<int, int>> exclusions_;
};

}
#endif //GROMACS_TOPOLOGY_H
