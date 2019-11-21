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

void MoleculeType::addExclusion(const int atomIndex, const int atomIndexToExclude)
{
    exclusions_.emplace_back(std::make_tuple(atomIndex, atomIndexToExclude));
}

void MoleculeType::addExclusion(std::string atomName, std::string atomNameToExclude)
{
    int indexAtomWithExclusion, indexAtomToExclude;

    auto iterWithExclusion = std::find_if(begin(atoms_), end(atoms_),
            [&atomName](std::tuple<std::string, std::string> & exclusionTuple)
            {
                if (std::get<0>(exclusionTuple) == atomName)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            });

    auto iterToExclude = std::find_if(begin(atoms_), end(atoms_),
              [&atomNameToExclude](std::tuple<std::string, std::string> & exclusionTuple)
              {
                  if (std::get<0>(exclusionTuple) == atomNameToExclude)
                  {
                      return true;
                  }
                  else
                  {
                      return false;
                  }
              });

    indexAtomWithExclusion = begin(atoms_) - iterWithExclusion;
    indexAtomToExclude = begin(atoms_) - iterToExclude;

    addExclusion(indexAtomWithExclusion, indexAtomToExclude);
}

} // namespace nblib
