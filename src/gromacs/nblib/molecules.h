//
// Created by sebkelle on 20.11.19.
//

#ifndef GROMACS_MOLECULES_H
#define GROMACS_MOLECULES_H

#include <tuple>
#include <unordered_map>
#include <string>

#include "gromacs/math/vectypes.h"
#include "interactions.h"

class TopologyBuilder;

using MoleculeAtomName = std::string;
using ResidueName = std::string;

class AtomType {
public:
    AtomType(std::string name,
             real mass,
             real charge,
             real c6,
             real c12);

    std::string getName() const;

private:
    std::string name_;

    real mass_;
    real charge_;
    real c6_;
    real c12_;
};

class MoleculeType {
public:
    MoleculeType(std::string name);

    MoleculeType& addAtom(std::string moleculeAtomName, std::string residueName, AtomType const& atom);
    MoleculeType& addAtom(std::string moleculeAtomName, AtomType const& atom);

    MoleculeType &addAtom(std::string name, AtomType const &atomType);

    void addHarmonicBond(int atom1, int atom2, HarmonicType harmonicType);

    int numAtomsInMolecule() const;

    friend class TopologyBuilder;

private:
    std::string name_;

    //! one entry per atom in molecule
    std::vector<std::tuple<MoleculeAtomName, ResidueName>> atoms_;
    //! collection of distinct AtomTypes in molecule
    std::unordered_map<MoleculeAtomName, AtomType> atomTypes_;

    std::vector<std::tuple<int, int>> exclusions;

    std::vector<HarmonicType> harmonicInteractions;

    int numAtomsInMolecule();
};

} //namespace nblib
#endif //GROMACS_MOLECULES_H
