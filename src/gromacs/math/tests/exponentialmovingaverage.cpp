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
 * Tests for the exponential moving average.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/exponentialmovingaverage.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace
{

TEST(ExponentialMovingAverage, ThrowsWhenLagTimeIsZero)
{
    EXPECT_THROW_GMX(ExponentialMovingAverage(0), InconsistentInputError);
}

TEST(ExponentialMovingAverage, ThrowsWhenLagTimeIsNegative)
{
    EXPECT_THROW_GMX(ExponentialMovingAverage(-10), InconsistentInputError);
}

TEST(ExponentialMovingAverage, LagTimeOneYieldsInstantaneousValue)
{
    const real               lagTime = 1;
    ExponentialMovingAverage exponentialMovingAverage(lagTime);

    exponentialMovingAverage.updateWithDataPoint(10);
    EXPECT_REAL_EQ(10, exponentialMovingAverage.biasCorrectedAverage());

    exponentialMovingAverage.updateWithDataPoint(-10);
    EXPECT_REAL_EQ(-10, exponentialMovingAverage.biasCorrectedAverage());
}

TEST(ExponentialMovingAverage, YieldsCorrectValue)
{
    const real               lagTime = 100;
    ExponentialMovingAverage exponentialMovingAverage(lagTime);

    exponentialMovingAverage.updateWithDataPoint(10);
    EXPECT_REAL_EQ(10, exponentialMovingAverage.biasCorrectedAverage());

    exponentialMovingAverage.updateWithDataPoint(-10);
    EXPECT_REAL_EQ(-0.050251256281406857, exponentialMovingAverage.biasCorrectedAverage());

    exponentialMovingAverage.updateWithDataPoint(0);
    EXPECT_REAL_EQ(-0.03333221103666531, exponentialMovingAverage.biasCorrectedAverage());
}

TEST(ExponentialMovingAverage, SetAverageCorrectly)
{
    const real               lagTime = 100;
    ExponentialMovingAverage exponentialMovingAverage(lagTime);

    exponentialMovingAverage.updateWithDataPoint(10);
    EXPECT_REAL_EQ(10, exponentialMovingAverage.biasCorrectedAverage());

    ExponentialMovingAverageState thisState = exponentialMovingAverage.state();

    ExponentialMovingAverage other(lagTime, thisState);

    other.updateWithDataPoint(-10);
    EXPECT_REAL_EQ(-0.050251256281406857, other.biasCorrectedAverage());

    other.updateWithDataPoint(0);
    EXPECT_REAL_EQ(-0.03333221103666531, other.biasCorrectedAverage());
}

TEST(ExponentialMovingAverage, DeterminesCorrectlyIfIncreasing)
{
    const real               lagTime = 100;
    ExponentialMovingAverage exponentialMovingAverage(lagTime);

    exponentialMovingAverage.updateWithDataPoint(10);
    exponentialMovingAverage.updateWithDataPoint(9.99);

    EXPECT_FALSE(exponentialMovingAverage.increasing());

    exponentialMovingAverage.updateWithDataPoint(-10);

    EXPECT_FALSE(exponentialMovingAverage.increasing());

    exponentialMovingAverage.updateWithDataPoint(100);
    EXPECT_TRUE(exponentialMovingAverage.increasing());
}


TEST(ExponentialMovingAverage, InverseLagTimeCorrect)
{
    const real               lagTime = 2.;
    ExponentialMovingAverage exponentialMovingAverage(lagTime);
    EXPECT_REAL_EQ(0.5, exponentialMovingAverage.inverseTimeConstant());
}

TEST(ExponentialMovingAverage, RoundTripAsKeyValueTree)
{
    KeyValueTreeBuilder           builder;
    const real                    weightedSum   = 9;
    const real                    weightedCount = 1;
    const bool                    increasing    = true;
    ExponentialMovingAverageState state         = { weightedSum, weightedCount, increasing };
    exponentialMovingAverageStateAsKeyValueTree(builder.rootObject(), state);
    state                     = {};
    KeyValueTreeObject result = builder.build();
    state                     = exponentialMovingAverageStateFromKeyValueTree(result);
    EXPECT_EQ(weightedSum, state.weightedSum_);
    EXPECT_EQ(weightedCount, state.weightedCount_);
    EXPECT_EQ(increasing, state.increasing_);
}

} // namespace

} // namespace test

} // namespace gmx
