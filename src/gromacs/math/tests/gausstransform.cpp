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
 * Tests for Gaussian spreading
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/gausstransform.h"

#include <cstddef>

#include <algorithm>
#include <array>
#include <functional>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/mdspan.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace gmx
{

namespace test
{

namespace
{

TEST(GaussianOn1DLattice, sumsCloseToOne)
{
    const int           spreadWidth = 7;
    const real          sigma       = 1.9;
    const real          shift       = 0.1;
    GaussianOn1DLattice gauss1d(spreadWidth, sigma);
    const real          amplitude = 1;
    gauss1d.spread(amplitude, shift);
    auto sumOverLattice = std::accumulate(std::begin(gauss1d.view()), std::end(gauss1d.view()), 0.);

    FloatingPointTolerance tolerance(defaultFloatTolerance());
    // The sum over the lattice should be roughly one
    EXPECT_FLOAT_EQ_TOL(0.99993300437927246, sumOverLattice, tolerance);
}

TEST(GaussianOn1DLattice, isCorrect)
{
    const int              spreadWidth = 2;
    const real             sigma       = 0.57;
    const real             shift       = -0.2;
    GaussianOn1DLattice    gauss1d(spreadWidth, sigma);
    const auto             viewOnResult = gauss1d.view();
    FloatingPointTolerance tolerance(defaultFloatTolerance());
    const real             amplitude = 1.0;
    gauss1d.spread(amplitude, shift);
    std::array<float, 2 * spreadWidth + 1> expected = { 0.0047816522419452667236328125,
                                                        0.2613909542560577392578125,
                                                        0.65811407566070556640625,
                                                        0.07631497085094451904296875,
                                                        0.000407583254855126142501831054688 };
    EXPECT_THAT(expected, Pointwise(FloatEq(tolerance), viewOnResult));
}

TEST(GaussianOn1DLattice, complementaryAmplitudesSumToZero)
{
    const int              spreadWidth = 2;
    const real             sigma       = 0.57;
    const real             shift       = -0.2;
    GaussianOn1DLattice    gauss1d(spreadWidth, sigma);
    const auto             viewOnResult = gauss1d.view();
    FloatingPointTolerance tolerance(defaultFloatTolerance());
    const real             amplitude = 2.0;
    gauss1d.spread(amplitude, shift);
    std::vector<float> sumOfComplementaryGaussians;
    // keep a copy of the first Gaussian
    std::copy(std::begin(viewOnResult), std::end(viewOnResult), std::back_inserter(sumOfComplementaryGaussians));

    gauss1d.spread(-amplitude, shift);
    // add the two spread Gaussians
    std::transform(std::begin(viewOnResult),
                   std::end(viewOnResult),
                   std::begin(sumOfComplementaryGaussians),
                   std::begin(sumOfComplementaryGaussians),
                   std::plus<>());
    // Expect all zeros
    std::array<float, 2 * spreadWidth + 1> expected = {};
    EXPECT_THAT(expected, Pointwise(FloatEq(tolerance), sumOfComplementaryGaussians));
}

TEST(GaussianOn1DLattice, doesNotOverflowForLargeRange)
{
    const int           spreadWidth = 200;
    const real          sigma       = 1;
    const real          shift       = -0.5;
    GaussianOn1DLattice gauss1d(spreadWidth, sigma);
    const auto          viewOnResult = gauss1d.view();
    const real          weightFirst  = 1;
    gauss1d.spread(weightFirst, shift);

    for (size_t i = 0; i < 187; i++)
    {
        EXPECT_FLOAT_EQ(0, viewOnResult[i]);
    }

    std::array<float, 12> expectedResult = {
        7.64165518045829568433117268116e-30, 4.57537550091150099862240391475e-25,
        1.00779356059942059652096780012e-20, 8.1662360418003390174698785664e-17,
        2.43432064870457987026952650922e-13, 2.66955679784075528004905208945e-10,
        1.07697594842193211661651730537e-07, 1.59837418323149904608726501465e-05,
        0.000872682663612067699432373046875, 0.01752830110490322113037109375,
        0.12951759994029998779296875,        0.3520653247833251953125
    };

    for (size_t i = 188; i < 200; i++)
    {
        EXPECT_FLOAT_EQ(expectedResult[i - 188], viewOnResult[i]);
    }

    for (size_t i = 200; i < 212; i++)
    {
        EXPECT_FLOAT_EQ(expectedResult[211 - i], viewOnResult[i]);
    }

    EXPECT_FLOAT_EQ(4.69519506491195688717869977423e-35, viewOnResult[212]);

    for (size_t i = 213; i < 2 * spreadWidth + 1; i++)
    {
        EXPECT_FLOAT_EQ(0, viewOnResult[i]);
    }
}

class GaussTransformTest : public ::testing::Test
{
public:
    void isZeroWithinFloatTolerance()
    {
        for (const auto& x : gaussTransform_.constView())
        {
            EXPECT_FLOAT_EQ_TOL(0, x, tolerance_);
        }
    }

protected:
    extents<dynamic_extent, dynamic_extent, dynamic_extent> latticeExtent_ = { 3, 3, 3 };
    DVec                                                    sigma_         = { 1., 1., 1. };
    double                                                  nSigma_        = 5;
    GaussTransform3D             gaussTransform_ = { latticeExtent_, { sigma_, nSigma_ } };
    test::FloatingPointTolerance tolerance_      = test::defaultFloatTolerance();
    const RVec                   latticeCenter_  = { 1, 1, 1 };
};

TEST_F(GaussTransformTest, isZeroUponConstruction)
{
    isZeroWithinFloatTolerance();
}

TEST_F(GaussTransformTest, isZeroAddingZeroAmplitudeGauss)
{
    gaussTransform_.add({ latticeCenter_, 0. });
    isZeroWithinFloatTolerance();
}

TEST_F(GaussTransformTest, isZeroAfterSettingZero)
{
    gaussTransform_.add({ latticeCenter_, 1. });
    for (const auto value : gaussTransform_.constView())
    {
        EXPECT_GT(value, 0);
    }
    gaussTransform_.setZero();
    isZeroWithinFloatTolerance();
}

TEST_F(GaussTransformTest, isZeroWhenOutsideRangeinX)
{
    DVec coordinateOutsideX(-nSigma_ * sigma_[XX], 0, 0);
    gaussTransform_.add({ coordinateOutsideX.toRVec(), 1. });
    isZeroWithinFloatTolerance();
}

TEST_F(GaussTransformTest, isZeroWhenOutsideRangeinY)
{
    DVec coordinateOutsideY(0, -nSigma_ * sigma_[YY], 0);
    gaussTransform_.add({ coordinateOutsideY.toRVec(), 1. });
    isZeroWithinFloatTolerance();
}

TEST_F(GaussTransformTest, isZeroWhenOutsideRangeinZ)
{
    RVec coordinateOutsideZ(0, 0, -nSigma_ * sigma_[ZZ]);
    gaussTransform_.add({ coordinateOutsideZ.toRVec(), 1. });
    isZeroWithinFloatTolerance();
}

TEST_F(GaussTransformTest, complementaryGaussAddToZero)
{
    gaussTransform_.add({ latticeCenter_, -2.0 });
    gaussTransform_.add({ latticeCenter_, 0.8 });
    gaussTransform_.add({ latticeCenter_, 1.2 });
    isZeroWithinFloatTolerance();
}

TEST_F(GaussTransformTest, centerGaussianInCubeHasExpectedValues)
{
    gaussTransform_.add({ latticeCenter_, 1. });
    const float        center         = 0.0634936392307281494140625;      // center
    const float        f              = 0.0385108403861522674560546875;   // face
    const float        e              = 0.0233580060303211212158203125;   // edge
    const float        c              = 0.014167346991598606109619140625; // corner
    std::vector<float> expectedValues = { c, e, c, e, f,      e, c, e, c,

                                          e, f, e, f, center, f, e, f, e,

                                          c, e, c, e, f,      e, c, e, c };
    // This assignment to std::vector is needed because the view() on the GaussTransform (aka basic_mdspan) does not provide size() needed by Pointwise
    std::vector<float> gaussTransformVector;
    gaussTransformVector.assign(gaussTransform_.constView().data(),
                                gaussTransform_.constView().data()
                                        + gaussTransform_.constView().mapping().required_span_size());
    EXPECT_THAT(expectedValues, testing::Pointwise(FloatEq(tolerance_), gaussTransformVector));
}

TEST_F(GaussTransformTest, view)
{
    gaussTransform_.add({ latticeCenter_, 1. });
    const float        f              = 0.0385108403861522674560546875;   // face
    const float        e              = 0.0233580060303211212158203125;   // edge
    const float        c              = 0.014167346991598606109619140625; // corner
    std::vector<float> expectedValues = { c, e, c, e, f,  e, c, e, c,

                                          e, f, e, f, 0., f, e, f, e,

                                          c, e, c, e, f,  e, c, e, c };

    gaussTransform_.view()[1][1][1] = 0;
    // This assignment to std::vector is needed because the view() on the
    // GaussTransform (aka basic_mdspan) does not provide size() that is
    // needed by Pointwise
    // \todo Update when mdspan functionality is extended
    std::vector<float> gaussTransformVector;
    gaussTransformVector.assign(
            gaussTransform_.view().data(),
            gaussTransform_.view().data() + gaussTransform_.view().mapping().required_span_size());
    EXPECT_THAT(expectedValues, testing::Pointwise(FloatEq(tolerance_), gaussTransformVector));
}

} // namespace

} // namespace test

} // namespace gmx
