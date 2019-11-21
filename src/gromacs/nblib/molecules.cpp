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

    // bool found = false;
    // for(auto tuple : exclusions_)
    // {
    //     using firstAtomIndex = std::get<0>(tuple);
    //     using secondAtomIndex = std::get<1>(tuple);

    //     bool found = false;
    //     if(firstAtomIndex == atomIndex) {
    //         if(secondAtomIndex == atomIndexToExclude) {
    //             found = true;
    //         }
    //     }
    // }
    // if(!found) {
    //     exclusions_.emplace_back(std::make_tuple(atomIndex, atomIndex));
    // }

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

int MoleculeType::atomNameToIndex(std::string atomName)
{
    auto iterAtom = std::find_if(begin(atoms_), end(atoms_),
        [&atomName](std::tuple<std::string, std::string> & atomAndResNameList)
        {
            if (std::get<0>(atomAndResNameList) == atomName)
            {
                return true;
            }
            else
            {
                return false;
            }
        });

    if(iterAtom - begin(atoms_) >= 0) {
        return iterAtom - begin(atoms_);
    } else {
        return 0;
        // TODO throw exception
    }
}

void MoleculeType::addExclusion(const int atomIndex, const int atomIndexToExclude)
{
    exclusions_.emplace_back(std::make_tuple(atomIndex, atomIndexToExclude));
}

void MoleculeType::addExclusion(std::string atomName, std::string atomNameToExclude)
{
    auto atomNameIndex = atomNameToIndex(atomName);
    auto atomToExcludeIndex = atomNameToIndex(atomNameToExclude);

    addExclusion(atomNameIndex, atomToExcludeIndex);
}

} // namespace nblib
