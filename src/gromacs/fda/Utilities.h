/*
 * Utilities.h
 *
 *  Created on: Nov 22, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_UTILITIES_H_
#define SRC_GROMACS_FDA_UTILITIES_H_

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"
#include "Vector2Scalar.h"

namespace fda {

/**
 * Takes a force vector v and returns its norm (magnitude) along with a sign
 * or the projection on the position vector formed by atoms i and j;
 * the sign will be negative (=attractive force) when the vector is in the same direction
 * as the position vector formed by atoms i and j and positive otherwise
 */
real vector2signedscalar(const rvec v, const rvec xi, const rvec xj, Vector2Scalar v2s);

real vector2unsignedscalar(const rvec v, int i, int j, PaddedRVecVector const& x);

} // namespace fda

#endif /* SRC_GROMACS_FDA_UTILITIES_H_ */
