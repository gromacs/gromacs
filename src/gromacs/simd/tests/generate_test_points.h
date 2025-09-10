/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

#ifndef GMX_SIMD_TESTS_GENERATE_TEST_POINTS_H
#define GMX_SIMD_TESTS_GENERATE_TEST_POINTS_H

#include <vector>

#include "gromacs/utility/real.h"

#include "simd.h"

namespace gmx
{
namespace test
{

/*! \brief Generate test point vector
 *
 *  \param inputRange  The test interval, half open. Upper limit is not included.
 *                     Pass by value, since we need to modify in method anyway.
 *  \param inputPoints Number of points to generate. This might be increased
 *                     slightly to account both for extra special values like 0.0
 *                     and the SIMD width.
 *
 * This routine generates a vector with test points separated by constant
 * multiplicative factors, based on the range and number of points in the
 * class. If the range includes both negative and positive values, points
 * will be generated separately for the negative/positive intervals down
 * to the smallest real number that can be represented, and we also include
 * 0.0 explicitly.
 *
 * This is highly useful for large test ranges. For example, with a linear
 * 1000-point division of the range (1,1e10) the first three values to test
 * would be 1, 10000000.999, and 20000000.998, etc. For large values we would
 * commonly hit the point where adding the small delta has no effect due to
 * limited numerical precision.
 * When we instead use this routine, the values will be 1, 1.0239, 1.0471, etc.
 * This will spread the entropy over all bits in the IEEE754 representation,
 * and be a much better test of all potential input values.
 */
std::vector<real> generateTestPoints(RealRange inputRange, std::size_t inputPoints);

} // namespace test
} // namespace gmx

#endif // GMX_SIMD_TESTS_GENERATE_TEST_POINTS_H
