/*
 * OnePair.h
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_ONEPAIR_H_
#define SRC_GROMACS_FDA_ONEPAIR_H_

#include <cstdint>
#include <iostream>
#include <string>

namespace fda {

/// OnePair defines the way the interactions between the same pair of atoms are stored
enum class OnePair : std::int8_t
{
    DETAILED, ///< each interaction is stored separately
              ///< it's possible to have the same pair in several interaction lists (default)
    SUMMED    ///< each interaction is stored once
};

/// Output stream for OnePair
std::ostream& operator << (std::ostream& os, OnePair r);

/// Input stream for OnePair
std::istream& operator >> (std::istream& is, OnePair& r);

} // namespace fda

#endif /* SRC_GROMACS_FDA_ONEPAIR_H_ */
