/*
 * ResultType.h
 *
 *  Created on: Nov 3, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_RESULTTYPE_H_
#define SRC_GROMACS_FDA_RESULTTYPE_H_

#include <cstdint>
#include <iostream>
#include <string>

namespace fda {

enum class ResultType : std::int8_t
{
    NO,                       // no storing (default)
    PAIRWISE_FORCES_VECTOR,
    PAIRWISE_FORCES_SCALAR,
    PUNCTUAL_STRESS,
    VIRIAL_STRESS,
    VIRIAL_STRESS_VON_MISES,
    COMPAT_BIN,               // DEPRICATED! compatibility mode (signed scalars) in binary
    COMPAT_ASCII              // DEPRICATED! compatibility mode (signed scalars) in ascii
};

/// Output stream for ResultType
std::ostream& operator << (std::ostream& os, ResultType r);

/// Input stream for ResultType
std::istream& operator >> (std::istream& is, ResultType& r);

} // namespace fda

#endif /* SRC_GROMACS_FDA_RESULTTYPE_H_ */
