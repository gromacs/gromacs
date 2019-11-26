//
// Created by sebkelle on 20.11.19.
//

#include <tuple>

#include "atomtype.h"
#include "molecules.h"

namespace nblib {


Molecule::Molecule(std::string moleculeName) : name_(std::move(moleculeName)) {}

void Molecule::addAtomSelfExclusion(std::string atomName, std::string resName)
{
    bool found = false;
    int atomIndex = atomNameAndResidueToIndex(std::make_tuple(atomName, resName));

    for(auto &tuple : exclusions_)
    {
        bool found = false;
        if(std::get<0>(tuple) == atomIndex) {
            if(std::get<1>(tuple) == atomIndex) {
                found = true;
            }
        }
    }
    if(!found) {
        exclusions_.emplace_back(std::make_tuple(atomIndex, atomIndex));
    }
}

Molecule& Molecule::addAtom(const AtomName& atomName, const ResidueName& residueName, const Charge& charge, AtomType const &atomType)
{
    // this check can only ensure that we don't use the same atomName twice, not that an (AtomType, charge) pair is not repeated
    if (!atomTypes_.count(atomName))
    {
        atomTypes_[atomName] = std::make_tuple(atomType, charge);
    }

    atoms_.emplace_back(std::make_tuple(atomName, residueName));
    addAtomSelfExclusion(atomName, residueName);

    return *this;
}

Molecule& Molecule::addAtom(const AtomName& atomName, const ResidueName& residueName, AtomType const &atomType)
{
    real charge = 0;
    addAtom(atomName, residueName, charge, atomType);

    return *this;
}

Molecule& Molecule::addAtom(const AtomName& atomName, const Charge& charge, AtomType const &atomType)
{
    addAtom(atomName, name_, charge, atomType);

    return *this;
}

Molecule& Molecule::addAtom(const AtomName& atomName, const AtomType& atomType)
{
    real charge = 0;
    addAtom(atomName, name_, charge, atomType);

    return *this;
}

int Molecule::numAtomsInMolecule() const
{
    return atoms_.size();
}

void Molecule::addHarmonicBond(HarmonicType harmonicBond)
{
    harmonicInteractions_.push_back(harmonicBond);
}


int Molecule::atomNameAndResidueToIndex(std::tuple<std::string, std::string> atomResNameTuple)
{
    auto equal = [](auto tup1, auto tup2) { return (std::get<0>(tup1) == std::get<0>(tup2) and std::get<1>(tup1) == std::get<1>(tup2)); };

    auto posIter = std::find_if(begin(atoms_), end(atoms_),
        [&](std::tuple<std::string, std::string> & atomAndResName)
        {
            return equal(atomAndResName, atomResNameTuple);
        });

    if(posIter != end(atoms_)) {
        return posIter - begin(atoms_);
    } else {
        // TODO throw exception
        return -1;
    }
}

void Molecule::addExclusion(const int atomIndex, const int atomIndexToExclude)
{
    // We do not need to add exclusion in case the atom indexes are the same
    // because self exclusion are added by addAtom
    if(atomIndex != atomIndexToExclude){
        exclusions_.emplace_back(std::make_tuple(atomIndex, atomIndexToExclude));
        exclusions_.emplace_back(std::make_tuple(atomIndexToExclude, atomIndex));
    }
}

void Molecule::addExclusion(std::tuple<std::string, std::string> atom, std::tuple<std::string, std::string> atomToExclude)
{
    auto atomNameIndex = atomNameAndResidueToIndex(atom);
    auto atomToExcludeIndex = atomNameAndResidueToIndex(atomToExclude);

    addExclusion(atomNameIndex, atomToExcludeIndex);
}

void Molecule::addExclusion(std::string atomName, std::string atomNameToExclude)
{
    auto atomNameIndex = atomNameAndResidueToIndex(std::make_tuple(atomName, name_));
    auto atomToExcludeIndex = atomNameAndResidueToIndex(std::make_tuple(atomNameToExclude, name_));

    addExclusion(atomNameIndex, atomToExcludeIndex);
}

} // namespace nblib
