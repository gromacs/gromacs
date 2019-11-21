//
// Created by sebkelle on 20.11.19.
//

#include "molecules.h"

namespace nblib {


MoleculeType::MoleculeType(std::string moleculeName) : name_(std::move(moleculeName)) {}

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

void MoleculeType::addExclusion(std::string atomWithExclusion, std::string atomToExclude)
{
    int indexAtomWithExclusion, indexAtomToExclude;

    auto iterWithExclusion = std::find_if(atoms_.begin(), atoms_.end(),
            [&atomWithExclusion](std::tuple<std::string, std::string> & val)
            {
                if (std::get<0>(val) == atomWithExclusion)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            });

    auto iterToExclude = std::find_if(atoms_.begin(), atoms_.end(),
              [&atomToExclude](std::tuple<std::string, std::string> & val)
              {
                  if (std::get<0>(val) == atomToExclude)
                  {
                      return true;
                  }
                  else
                  {
                      return false;
                  }
              });

    indexAtomWithExclusion = atoms_.begin() - iterWithExclusion;
    indexAtomToExclude = atoms_.begin() - iterToExclude;

    addExclusion(indexAtomWithExclusion, indexAtomToExclude);
}

} // namespace nblib
