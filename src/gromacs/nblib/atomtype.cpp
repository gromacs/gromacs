//
// Created by sebkelle on 20.11.19.
//

#include "atomtype.h"

namespace nblib {

AtomType::AtomType() noexcept :
  name_(""),
  mass_(0),
  c6_(0),
  c12_(0)
{}

AtomType::AtomType(std::string atomName, real mass, real c6, real c12)
: name_(std::move(atomName)),
  mass_(mass),
  c6_(c6),
  c12_(c12)
{}

std::string AtomType::name() const { return name_; }

real AtomType::mass() const { return mass_; }

real AtomType::c6() const { return c6_; }

real AtomType::c12() const { return c12_; }

} // namespace nblib
