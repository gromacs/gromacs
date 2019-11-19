//
// Created by sebkelle on 19.11.19.
//

#ifndef GROMACS_TOPOLOGY_H
#define GROMACS_TOPOLOGY_H


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

    void setNonbondedParameters(std::vector<int> params);

    void setAtomTypes(std::vector<int> types);

    void setCharges(std::vector<int> charges);

    void setMasses(std::vector<int> masses);

    //! hardcoded version to converto to t_blocka
    void setExclusions(std::vector<int>, std::vector<int>);

    //! set exclusion rules (molecules or connectivity)
    //void setExclusions(exclusionRules);

    //! set exclusion rules based on a tuple
    //void setExclusions(someTupleWhatHaveYou);

};

#endif //GROMACS_TOPOLOGY_H
