//
// Created by sebkelle on 20.11.19.
//

#ifndef GROMACS_ATOMS_H
#define GROMACS_ATOMS_H

#include <tuple>
#include <unordered_map>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "interactions.h"

class TopologyBuilder;

namespace nblib
{

class Atom {
public:
    Atom() noexcept;

    Atom(std::string atomName,
             real mass,
             real charge,
             real c6,
             real c12);

    std::string name() const;

    real mass() const;

    real charge() const;

    real c6() const;

    real c12() const;

private:
    std::string name_;
    real mass_;
    real charge_;
    real c6_;
    real c12_;
};

} //namespace nblib
#endif //GROMACS_MOLECULES_H
