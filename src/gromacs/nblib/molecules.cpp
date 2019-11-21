//
// Created by sebkelle on 20.11.19.
//

#include "molecules.h"

namespace nblib {


MoleculeType::MoleculeType(std::string name) : name_(std::move(name)) {}

MoleculeType& MoleculeType::addAtom(const std::string &atomName, const std::string &residueName, AtomType const &atomType)
{
    // check whether we already have the atom type
    if (!atomTypes_.count(atomName))
    {
        atomTypes_[atomName] = atomType;
    }

    atoms_.emplace_back(std::make_tuple(atomName, residueName));

    return *this;
}

MoleculeType& MoleculeType::addAtom(const std::string &atomName, AtomType const &atomType)
{
    if (name_.length() > 0){
        return this->addAtom(atomName, atomName, atomType);
    }
    else {
        return this->addAtom(atomName, "", atomType);
    }
}

int MoleculeType::numAtomsInMolecule() const
{
    return atoms_.size();
}

void MoleculeType::addHarmonicBond(HarmonicType harmonicBond)
{
    harmonicInteractions_.push_back(harmonicBond);
}

void MoleculeType::addExclusion(const int atomWithExclusion, const int atomToExclude)
{
    exclusions_.emplace_back(std::make_tuple(atomWithExclusion, atomToExclude));
}

} // namespace nblib
