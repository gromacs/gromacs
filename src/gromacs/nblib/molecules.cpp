//
// Created by sebkelle on 20.11.19.
//

#include "molecules.h"

namespace nblib {

AtomType::AtomType()
: atomName_(""),
  mass_(0),
  charge_(0),
  c6_(0),
  c12_(0)
{}

AtomType::AtomType(std::string atomName, real mass, real charge, real c6, real c12)
: atomName_(std::move(atomName)),
  mass_(mass),
  charge_(charge),
  c6_(c6),
  c12_(c12)
{}

std::string AtomType::getName() const { return atomName_; }


MoleculeType::MoleculeType(std::string moleculeName) : moleculeName_(std::move(moleculeName)) {}

MoleculeType& MoleculeType::addAtom(const std::string &moleculeAtomName, const std::string &residueName, AtomType const &atom)
{
    // check whether we already have the atom type
    if (!atomTypes_.count(moleculeAtomName))
    {
        atomTypes_[moleculeAtomName] = atom;
    }

    atoms_.emplace_back(std::make_tuple(moleculeAtomName, residueName));

    return *this;
}

MoleculeType& MoleculeType::addAtom(const std::string &name, AtomType const &atom)
{
    return this->addAtom(name, "", atom);
}

int MoleculeType::numAtomsInMolecule() const
{
    return atoms_.size();
}

void MoleculeType::addHarmonicBond(HarmonicType harmonicBond)
{
    harmonicInteractions_.push_back(harmonicBond);
}

} // namespace nblib
