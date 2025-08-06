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

#include "generate_test_points.h"

#include <algorithm>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#include "simd.h"

namespace gmx
{
namespace test
{

std::vector<real> generateTestPoints(RealRange inputRange, std::size_t inputPoints)
{
    std::vector<real> testPoints;
    testPoints.reserve(inputPoints);

    GMX_RELEASE_ASSERT(inputRange.first < inputRange.second,
                       "The start of the interval must come before the end");

    std::vector<RealRange> testRanges;

    if (inputRange.first < 0 && inputRange.second > 0)
    {
        testRanges.emplace_back(inputRange.first, -std::numeric_limits<real>::min());
        testRanges.emplace_back(0.0, inputRange.second);
    }
    else
    {
        if (inputRange.second == 0)
        {
            inputRange.second = -std::numeric_limits<real>::min();
            inputRange.first  = std::min(inputRange.first, inputRange.second);
        }
        testRanges.push_back(inputRange);
    }

    for (RealRange& range : testRanges)
    {
        std::size_t points = inputPoints / testRanges.size();

        // The value 0 is special, and can only occur at the start of
        // the interval after the corrections outside this loop.
        // Add it explicitly, and adjust the interval to continue
        // at the first valid non-zero positive number.
        if (range.first == 0)
        {
            testPoints.push_back(0.0);
            range.first = std::numeric_limits<real>::min();
            points--; // Used one point
        }

        union
        {
            real                                                                               r;
            std::conditional<sizeof(real) == sizeof(double), std::int64_t, std::int32_t>::type i;
        } low, high, x;

        low.r  = range.first;
        high.r = range.second;

        // IEEE754 floating-point numbers have the cool property that for any range of
        // constant sign, for all non-zero numbers a constant (i.e., linear) difference
        // in the bitwise representation corresponds to a constant multiplicative factor.
        //
        // Divide the ulp difference evenly
        std::int64_t ulpDiff = high.i - low.i;
        // dividend and divisor must both be signed types
        std::int64_t ulpDelta    = ulpDiff / static_cast<std::int64_t>(points);
        std::int64_t minUlpDelta = (ulpDiff > 0) ? 1 : -1;

        if (ulpDelta == 0)
        {
            // Very short interval or very many points caused round-to-zero.
            // Select the smallest possible change, which is one ulp (with correct sign)
            ulpDelta = minUlpDelta;
            points   = std::abs(ulpDiff);
        }

        x.r = low.r;
        // Use an index-based loop to avoid floating-point comparisons with
        // values that might have overflowed. Save one point for the very last
        // bitwise value that is part of the interval
        for (std::size_t i = 0; i < points - 1; i++)
        {
            testPoints.push_back(x.r);
            x.i += ulpDelta;
        }

        // Make sure we test the very last point that is inside the interval
        x.r = high.r;
        x.i -= minUlpDelta;
        testPoints.push_back(x.r);
    }
    return testPoints;
}

} // namespace test
} // namespace gmx
