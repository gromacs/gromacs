//
// Created by sebkelle on 19.11.19.
//

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/topology/block.h"
#include "molecules.h"

#ifndef GROMACS_TOPOLOGY_H
#define GROMACS_TOPOLOGY_H

struct t_blocka;

namespace nblib {

class Topology {
public:

    const std::vector<int>& getAtoms() const;

    const std::vector<real>& getCharges() const;

    const std::vector<real>& getMasses() const;

    const std::vector<real>& getNonbondedParameters() const;

    const std::vector<int>& getAtomInfoAllVdw() const;

    // TODO: This function is only needed for testing. Need
    //       another way for testing exclusion correctness
    const t_blocka& getGMXexclusions() const
    {
        return excls;
    }

private:
    Topology() = default;

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
    TopologyBuilder() = default;

    Topology buildTopology();

    TopologyBuilder& addMolecule(Molecule moleculeType, int nMolecules);

private:
    Topology topology_;

    int numAtoms_;
    std::vector<std::tuple<Molecule, int>> molecules_;

    t_blocka createExclusionsList() const;

    template <class Extractor>
    std::vector<real> extractAtomTypeQuantity(Extractor extractor);

    std::vector<real> extractCharge();
};

} // namespace nblib

#endif //GROMACS_TOPOLOGY_H
