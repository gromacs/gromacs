/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * Tests for energyanalysis energy term
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/energyanalysis/energyterm.h"

#include <cstring>

#include <optional>
#include <string>

#include <gtest/gtest.h>

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(EnergyTermTest, ConstructWorks)
{
    EnergyTerm term(0, true, "test", "test");
    EXPECT_FALSE(term.slopeOfLinearFit().has_value());
    EXPECT_FALSE(term.errorEstimate(0).has_value());
}

TEST(EnergyTermTest, AddFrameWorks)
{
    EnergyTerm term(0, true, "test", "test");
    term.addFrame(2, 1000, 10, 50, 5, 255);
    term.addFrame(4, 2000, 10, 100, 10, 155);
    EXPECT_REAL_EQ(term.average(), 7.5);
    EXPECT_REAL_EQ(term.standardDeviation(), 2.6457513110645907);
    EXPECT_EQ(term.numFrames(), 2);
    EXPECT_EQ(term.numSteps(), 1000);
    EXPECT_REAL_EQ(term.timeSpan(), 2);
    auto errorEstimate = term.errorEstimate(1);
    ASSERT_TRUE(errorEstimate.has_value());
    EXPECT_REAL_EQ(errorEstimate.value(), 0);
    term.addFrame(6, 3000, 10, 75, 7, 175);
    auto slope = term.slopeOfLinearFit();
    ASSERT_TRUE(slope.has_value());
    EXPECT_REAL_EQ(slope.has_value(), 1);
}

} // namespace
} // namespace test
} // namespace gmx
