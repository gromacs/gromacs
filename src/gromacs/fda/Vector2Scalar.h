/*
 * Vector2Scalar.h
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_VECTOR2SCALAR_H_
#define SRC_GROMACS_FDA_VECTOR2SCALAR_H_

#include <cstdint>
#include <iostream>
#include <string>

namespace fda {

/// Vector2Scalar defines the way the force vector to scalar conversion is done
/// As the conversion affects all scalar writing modes (PF_FILE_OUT*), this is kept as a separate
/// setting rather than creating separates modes for norm and projection.
enum class Vector2Scalar : std::int8_t
{
    NORM,       ///< Takes the norm of the vector
    PROJECTION  ///< Takes the projection of the force on the direction of the 2 atoms
};

/// Output stream for Vector2Scalar
std::ostream& operator << (std::ostream& os, Vector2Scalar r);

/// Input stream for Vector2Scalar
std::istream& operator >> (std::istream& is, Vector2Scalar& r);

} // namespace fda

#endif /* SRC_GROMACS_FDA_VECTOR2SCALAR_H_ */
