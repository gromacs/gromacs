//
// Created by sebkelle on 20.11.19.
//

#include "molecules.h"

namespace nblib {

AtomType::AtomType(std::string name, real mass, real charge, real c6, real c12)
: name_(name),
  mass_(mass),
  charge_(charge),
  c6_(c6),
  c12_(c12)
{}

std::string AtomType::getName() const { return name_; }


MoleculeType::MoleculeType(std::string name) name_(name) {}

MoleculeType& MoleculeType::addAtom(std::string moleculeAtomName, std::string residueName, AtomType const &atom)
{
    // check whether we already have the atom type
    if (!atomTypes_.count(moleculeAtomName))
    {
        atomTypes_[moleculeAtomName] = atom;
    }

    atoms_.push_back(std::make_tuple(moleculeAtomName, residueName));
}

MoleculeType& MoleculeType::addAtom(std::string name, AtomType const &atom)
{
    this->addAtom(name, "", atom);
}

int MoleculeType::numAtomsInMolecule()
{
    return atoms_.size();
}

} // namespace nblib
