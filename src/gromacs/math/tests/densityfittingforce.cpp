/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Tests for the density fitting force.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/densityfittingforce.h"

#include <cmath>

#include <algorithm>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/math/multidimarray.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"
namespace gmx
{

namespace test
{

namespace
{

TEST(DensityFittingForce, isZeroWhenMatchingDensity)
{
    const float sigma      = 1;
    const float nSigma     = 5;
    const RVec  gridCenter = { 1, 1, 1 };

    DensityFittingForce forceEvaluator({ { sigma, sigma, sigma }, nSigma });
    MultiDimArray<std::vector<float>, dynamicExtents3D> densityDerivative(3, 3, 3);
    std::fill(begin(densityDerivative), end(densityDerivative), 0);

    const RVec result =
            forceEvaluator.evaluateForce({ gridCenter, 1.0 }, densityDerivative.asConstView());
    const RVec             expected = { 0, 0, 0 };
    FloatingPointTolerance tolerance(defaultFloatTolerance());

    EXPECT_FLOAT_EQ_TOL(expected[XX], result[XX], tolerance);
    EXPECT_FLOAT_EQ_TOL(expected[YY], result[YY], tolerance);
    EXPECT_FLOAT_EQ_TOL(expected[ZZ], result[ZZ], tolerance);
}

TEST(DensityFittingForce, isZeroWhenMismatchingSameAllDirections)
{
    const BasicVector<double> sigma      = { 1., 1., 1. };
    const float               nSigma     = 5;
    const RVec                gridCenter = { 1, 1, 1 };

    DensityFittingForce                                 forceEvaluator({ { sigma }, nSigma });
    MultiDimArray<std::vector<float>, dynamicExtents3D> densityDerivative(3, 3, 3);
    std::fill(begin(densityDerivative), end(densityDerivative), 1);

    const RVec result =
            forceEvaluator.evaluateForce({ gridCenter, 1. }, densityDerivative.asConstView());
    const RVec expected = { 0, 0, 0 };

    // need to increase the tolerance here, because the checked value is the
    // result of a complex calculation (the sum of 27 3d Gaussians weighted with the distance to
    // center) \todo increase tightness of this bound with better implementations
    FloatingPointTolerance tolerance(absoluteTolerance(1e-8));

    EXPECT_FLOAT_EQ_TOL(expected[XX], result[XX], tolerance);
    EXPECT_FLOAT_EQ_TOL(expected[YY], result[YY], tolerance);
    EXPECT_FLOAT_EQ_TOL(expected[ZZ], result[ZZ], tolerance);
}

TEST(DensityFittingForce, pullsTowardsDerivative)
{
    const float sigma      = 1;
    const float nSigma     = 5;
    const RVec  gridCenter = { 1, 1, 1 };

    DensityFittingForce forceEvaluator({ { sigma, sigma, sigma }, nSigma });
    MultiDimArray<std::vector<float>, dynamicExtents3D> densityDerivative(3, 3, 3);
    std::fill(begin(densityDerivative), end(densityDerivative), 0);
    densityDerivative(0, 0, 0) = 1;

    const RVec result =
            forceEvaluator.evaluateForce({ gridCenter, 1. }, densityDerivative.asConstView());
    real expected = -1 / std::sqrt(2 * 2 * 2 * M_PI * M_PI * M_PI) * std::exp(-0.5) * std::exp(-0.5)
                    * std::exp(-0.5);

    EXPECT_FLOAT_EQ(expected, result[XX]);
    EXPECT_FLOAT_EQ(expected, result[YY]);
    EXPECT_FLOAT_EQ(expected, result[ZZ]);
}

} // namespace

} // namespace test

} // namespace gmx
