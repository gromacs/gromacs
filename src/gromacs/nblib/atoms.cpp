//
// Created by sebkelle on 20.11.19.
//

#include "atoms.h"

namespace nblib {

Atom::Atom() noexcept :
  name_(""),
  mass_(0),
  charge_(0),
  c6_(0),
  c12_(0)
{}

Atom::Atom(std::string atomName, real mass, real charge, real c6, real c12)
: name_(std::move(atomName)),
  mass_(mass),
  charge_(charge),
  c6_(c6),
  c12_(c12)
{}

std::string Atom::name() const { return name_; }

real Atom::mass() const { return mass_; }

real Atom::charge() const { return charge_; }

real Atom::c6() const { return c6_; }

real Atom::c12() const { return c12_; }

} // namespace nblib
