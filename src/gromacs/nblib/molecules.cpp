//
// Created by sebkelle on 20.11.19.
//

#include <tuple>

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

Molecule& Molecule::addAtom(const std::string &atomName, const std::string &residueName, Atom const &atomType)
{
    // check whether we already have the atom type
    if (!atomTypes_.count(atomType.name()))
    {
        atomTypes_[atomName] = atomType;
    }

    atoms_.emplace_back(std::make_tuple(atomName, residueName));
    addAtomSelfExclusion(atomName, residueName);

    return *this;
}

Molecule& Molecule::addAtom(const std::string &atomName, Atom const &atomType)
{
    addAtom(atomName, name_, atomType);
    addAtomSelfExclusion(atomName, name_);

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
    // TODO
    // avoid duplicating the exclusions
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
