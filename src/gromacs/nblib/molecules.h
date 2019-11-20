//
// Created by sebkelle on 20.11.19.
//

#ifndef GROMACS_MOLECULES_H
#define GROMACS_MOLECULES_H

#include <tuple>
#include <unordered_map>

#include "gromacs/math/vectypes.h"
#include "interactions.h"

class TopologyBuilder;

using AtomName = std::string;
using ResidueName = std::string;

class AtomType
{
public:
    AtomType(std::string name,
             real mass,
             real charge,
             real c6,
             real c12);

private:
    std::string name_;

    real mass_;
    real charge_;
    real c6_;
    real c12_;
};

class MoleculeType
{
public:
    MoleculeType(int nAtoms);

    MoleculeType& addAtom(std::string name, std::string residueName, AtomType const& atomType);
    MoleculeType& addAtom(std::string name, AtomType const& atomType);

    void addHarmonicBond(int atom1, int atom2, HarmonicType harmonicType);

    friend class TopologyBuilder;

private:
    std::string name_;

    //! one entry per atom in molecule
    std::vector<std::tuple<AtomName, ResidueName>> atoms_;
    //! list of distinct AtomTypes in molecule
    std::unordered_map<AtomName, AtomType> atomsTypes_;

    std::vector<std::tuple<int, int>> exclusions;

    std::vector<HarmonicType>         harmonicInteractions;
};

#endif //GROMACS_MOLECULES_H
