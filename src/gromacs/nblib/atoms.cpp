//
// Created by sebkelle on 20.11.19.
//

#include "atoms.h"

namespace nblib {

AtomType::AtomType() noexcept :
  name_(""),
  mass_(0),
  charge_(0),
  c6_(0),
  c12_(0)
{}

AtomType::AtomType(std::string atomName, real mass, real charge, real c6, real c12)
: name_(std::move(atomName)),
  mass_(mass),
  charge_(charge),
  c6_(c6),
  c12_(c12)
{}

std::string AtomType::getName() const
{
    return name_;
}

} // namespace nblib
