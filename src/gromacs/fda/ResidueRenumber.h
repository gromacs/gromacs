/*
 * ResiduesRenumber.h
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_RESIDUESRENUMBER_H_
#define SRC_GROMACS_FDA_RESIDUESRENUMBER_H_

#include <cstdint>
#include <iostream>
#include <string>

namespace fda {

enum class ResiduesRenumber : std::int8_t
{
    AUTO,
    DO,
    DONT
};

/// Output stream for ResiduesRenumber
std::ostream& operator << (std::ostream& os, ResiduesRenumber r);

/// Input stream for ResiduesRenumber
std::istream& operator >> (std::istream& is, ResiduesRenumber& r);

} // namespace fda

#endif /* SRC_GROMACS_FDA_RESIDUESRENUMBER_H_ */
