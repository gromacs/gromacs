//
// Created by sebkelle on 20.11.19.
//

#ifndef GROMACS_MOLECULES_H
#define GROMACS_MOLECULES_H

#include <tuple>
#include <unordered_map>
#include <string>
#include <vector>

#include "atoms.h"
#include "interactions.h"

#include "gromacs/math/vectypes.h"

class TopologyBuilder;

namespace nblib
{

class MoleculeType {
public:
    MoleculeType(std::string moleculeName);

    MoleculeType& addAtom(const std::string &atomName, const std::string &residueName, AtomType const &atomType);

    MoleculeType& addAtom(const std::string &atomName, AtomType const &atomType);

    void addHarmonicBond(HarmonicType harmonicBond);

    void addExclusion(const int atomWithExclusion, const int atomToExclude);

    int numAtomsInMolecule() const;

    friend class TopologyBuilder;

private:
    std::string name_;

    //! one entry per atom in molecule
    std::vector<std::tuple<std::string, std::string>> atoms_;
    //! collection of distinct AtomTypes in molecule
    std::unordered_map<std::string, AtomType> atomTypes_;

    std::vector<std::tuple<int, int>> exclusions_;

    std::vector<HarmonicType> harmonicInteractions_;
};

} //namespace nblib
#endif //GROMACS_MOLECULES_H
