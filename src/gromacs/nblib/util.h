//
// Created by sebkelle on 19.11.19.
//

#ifndef GROMACS_UTIL_H
#define GROMACS_UTIL_H

#include <vector>

#include "gromacs/math/vectypes.h"

namespace nblib {

//! generate Velocites from a Maxwell Boltzmann distro, masses should be the
//! same as the ones specified for the Topology object
std::vector<gmx::RVec> generateVelocity(real Temperature, unsigned int seed, std::vector<real> const& masses);

bool checkNumericValues(const std::vector<gmx::RVec> &values);

} // namespace nblib

#endif //GROMACS_UTIL_H
