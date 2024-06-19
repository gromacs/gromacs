/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 *
 * \brief Tests routines in neldermead.h .
 *
 * \author Christian Blau <blau@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/math/neldermead.h"

#include <functional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{


/* \brief Test the NelderMeadSimplex class.
 *
 * Use a one-dimensional function to create a simplex with two vertices.
 * The simplex is created at [1, 1.05] with function values [2, 2.1]
 * respectively.
 */
class NelderMeadSimplexTest : public ::testing::Test
{
public:
    //! Set up a one-d, two-vertices Nelder-Mead simplex
    NelderMeadSimplexTest() :
        initialGuess_{ 1 }, simplex_{ NelderMeadSimplexTest::doubleFirstCoordinateValue, initialGuess_ }
    {
    }

    //! Function that the simplex is meant help optimise
    static real doubleFirstCoordinateValue(ArrayRef<const real> x) { return 2 * x[0]; }

protected:
    std::vector<real> initialGuess_;
    NelderMeadSimplex simplex_;
};

TEST_F(NelderMeadSimplexTest, BestVertex)
{
    EXPECT_REAL_EQ(simplex_.bestVertex().coordinate_[0], 1);
    EXPECT_REAL_EQ(simplex_.bestVertex().value_, 2);
}

TEST_F(NelderMeadSimplexTest, WorstVertex)
{
    EXPECT_REAL_EQ(simplex_.worstVertex().coordinate_[0], 1.05);
    EXPECT_REAL_EQ(simplex_.worstVertex().value_, 2.1);
}

TEST_F(NelderMeadSimplexTest, SecondWorstValue)
{
    EXPECT_REAL_EQ(simplex_.secondWorstValue(), 2.00);
}

TEST_F(NelderMeadSimplexTest, ReflectionPoint)
{
    EXPECT_REAL_EQ(
            simplex_.evaluateReflectionPoint(NelderMeadSimplexTest::doubleFirstCoordinateValue).coordinate_[0],
            0.95);
    EXPECT_REAL_EQ(
            simplex_.evaluateReflectionPoint(NelderMeadSimplexTest::doubleFirstCoordinateValue).value_, 1.9);
}

TEST_F(NelderMeadSimplexTest, EvaluateExpansionPoint)
{
    EXPECT_REAL_EQ(
            simplex_.evaluateExpansionPoint(NelderMeadSimplexTest::doubleFirstCoordinateValue).coordinate_[0],
            0.9);
    EXPECT_REAL_EQ(
            simplex_.evaluateExpansionPoint(NelderMeadSimplexTest::doubleFirstCoordinateValue).value_, 1.8);
}

TEST_F(NelderMeadSimplexTest, EvaluateContractionPoint)
{
    EXPECT_REAL_EQ(
            simplex_.evaluateContractionPoint(NelderMeadSimplexTest::doubleFirstCoordinateValue).coordinate_[0],
            1.025);
    EXPECT_REAL_EQ(
            simplex_.evaluateContractionPoint(NelderMeadSimplexTest::doubleFirstCoordinateValue).value_,
            2.05);
}


TEST_F(NelderMeadSimplexTest, SwapOutWorst)
{
    // introduce a new vertex that we know is better than any in the
    // initial simplex
    RealFunctionvalueAtCoordinate newVertex = { { 0 }, 0 };
    simplex_.swapOutWorst(newVertex);
    EXPECT_REAL_EQ(simplex_.bestVertex().coordinate_[0], 0);
    EXPECT_REAL_EQ(simplex_.bestVertex().value_, 0);
    // introduce a new vertex that we know is worse than any in the
    // initial simplex
    newVertex = { { 3 }, 6 };
    simplex_.swapOutWorst(newVertex);
    EXPECT_REAL_EQ(simplex_.worstVertex().coordinate_[0], 3);
    EXPECT_REAL_EQ(simplex_.worstVertex().value_, 6);
    // check that also the reflection point has changed now
    EXPECT_REAL_EQ(
            simplex_.evaluateReflectionPoint(NelderMeadSimplexTest::doubleFirstCoordinateValue).coordinate_[0],
            -3);
    EXPECT_REAL_EQ(
            simplex_.evaluateReflectionPoint(NelderMeadSimplexTest::doubleFirstCoordinateValue).value_, -6);
}

TEST_F(NelderMeadSimplexTest, ShrinkSimplexPointsExceptBest)
{
    simplex_.shrinkSimplexPointsExceptBest(NelderMeadSimplexTest::doubleFirstCoordinateValue);
    EXPECT_REAL_EQ(simplex_.worstVertex().coordinate_[0], 1.025);
    EXPECT_REAL_EQ(simplex_.worstVertex().value_, 2.05);
    // check that also the reflection point has changed now
    EXPECT_REAL_EQ(
            simplex_.evaluateReflectionPoint(NelderMeadSimplexTest::doubleFirstCoordinateValue).coordinate_[0],
            0.975);
    EXPECT_REAL_EQ(
            simplex_.evaluateReflectionPoint(NelderMeadSimplexTest::doubleFirstCoordinateValue).value_, 1.95);
}

TEST_F(NelderMeadSimplexTest, OrientedLength)
{
    // here, the distance between the two vertices
    EXPECT_REAL_EQ(simplex_.orientedLength(), 0.05);
}

} // namespace
} // namespace test
} // namespace gmx
