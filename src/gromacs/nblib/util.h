//
// Created by sebkelle on 19.11.19.
//

#ifndef GROMACS_UTIL_H
#define GROMACS_UTIL_H

//! generate Velocites from a Maxwell Boltzmann distro, masses should be the
//! same as the ones specified for the Topology object
void generateVelocity(real Temperature, std::vector<real> const& masses);

void checkNumericValues(const std::vector<gmx::RVec> &values);

#endif //GROMACS_UTIL_H
