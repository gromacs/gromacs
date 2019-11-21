//
// Created by sebkelle on 20.11.19.
//

#ifndef GROMACS_MOLECULES_H
#define GROMACS_MOLECULES_H

#include <tuple>
#include <unordered_map>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "interactions.h"

class TopologyBuilder;

namespace nblib
{

class AtomType {
public:
    AtomType();

    AtomType(std::string atomName,
             real mass,
             real charge,
             real c6,
             real c12);

    std::string getName() const;

private:
    std::string atomName_;

    real mass_;
    real charge_;
    real c6_;
    real c12_;
};

class MoleculeType {
public:
    explicit MoleculeType(std::string moleculeName);

    MoleculeType &addAtom(const std::string &moleculeAtomName, const std::string &residueName, AtomType const &atom);

    MoleculeType &addAtom(const std::string &moleculeAtomName, AtomType const &atom);

    void addExclusion(int atomWithExclusion, int atomToExclude);

    void addHarmonicBond(HarmonicType harmonicBond);

    int numAtomsInMolecule() const;

    friend class TopologyBuilder;

private:
    std::string moleculeName_;

    //! one entry per atom in molecule
    std::vector<std::tuple<std::string, std::string>> atoms_;
    //! collection of distinct AtomTypes in molecule
    std::unordered_map<std::string, AtomType> atomTypes_;

    std::vector<std::tuple<int, int>> exclusions_;

    std::vector<HarmonicType> harmonicInteractions_;
};

} //namespace nblib
#endif //GROMACS_MOLECULES_H
