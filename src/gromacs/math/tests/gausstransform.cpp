/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Tests canonical vector basis
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/gausstransform.h"

#include <array>
#include <numeric>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

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
    const int              spreadWidth = 7;
    const real             sigma       = 1.9;
    const real             shift       = 0.1;
    GaussianOn1DLattice    gauss1d(spreadWidth, sigma);
    const real             amplitude = 1;
    gauss1d.spread(amplitude, shift);
    auto                   sumOverLattice = std::accumulate(std::begin(gauss1d.view()), std::end(gauss1d.view()), 0.);

    FloatingPointTolerance tolerance(defaultFloatTolerance());
    // The sum over the lattice should be roughly one
    EXPECT_FLOAT_EQ_TOL(0.99993300437927246, sumOverLattice, tolerance);
}

TEST(GaussianOn1DLattice, isCorrect)
{
    const int                             spreadWidth = 2;
    const real                            sigma       = 0.57;
    const real                            shift       = -0.2;
    GaussianOn1DLattice                   gauss1d(spreadWidth, sigma);
    const auto                            viewOnResult = gauss1d.view();
    FloatingPointTolerance                tolerance(defaultFloatTolerance());
    const real                            amplitude = 1.0;
    gauss1d.spread(amplitude, shift);
    std::array<float, 2 *spreadWidth + 1> expected =  {
        0.0047816522419452667236328125, 0.2613909542560577392578125,
        0.65811407566070556640625, 0.07631497085094451904296875, 0.000407583254855126142501831054688
    };
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
    std::vector<float>     sumOfComplementaryGaussians;
    // keep a copy of the first Gaussian
    std::copy(std::begin(viewOnResult), std::end(viewOnResult),
              std::back_inserter(sumOfComplementaryGaussians));

    gauss1d.spread(-amplitude, shift);
    // add the two spread Gaussians
    std::transform(std::begin(viewOnResult), std::end(viewOnResult),
                   std::begin(sumOfComplementaryGaussians), std::begin(sumOfComplementaryGaussians), std::plus<>());
    // Expect all zeros
    std::array<float, 2 *spreadWidth + 1> expected = {};
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
        0.12951759994029998779296875, 0.3520653247833251953125
    };

    for (size_t i = 188; i < 200; i++)
    {
        EXPECT_FLOAT_EQ(expectedResult[i-188], viewOnResult[i]);
    }

    for (size_t i = 200; i < 212; i++)
    {
        EXPECT_FLOAT_EQ(expectedResult[211-i], viewOnResult[i]);
    }

    EXPECT_FLOAT_EQ(4.69519506491195688717869977423e-35, viewOnResult[212]);

    for (size_t i = 213; i < 2*spreadWidth+1; i++)
    {
        EXPECT_FLOAT_EQ(0, viewOnResult[i]);
    }
}

} // namespace

} // namespace test

} // namespace gmx
