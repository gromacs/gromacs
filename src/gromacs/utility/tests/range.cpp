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
 * \brief
 * Tests for the Range class.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/range.h"

#include <cstddef>

#include <algorithm>
#include <limits>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(Range, EmptyRangeWorks)
{
    Range<int> range;

    EXPECT_EQ(range.empty(), true);
    EXPECT_EQ(range.size(), 0);
}

TEST(Range, NonEmptyRangeWorks)
{
    const Range<char> range(3, 5);

    EXPECT_EQ(range.empty(), false);
    EXPECT_EQ(range.size(), 2);
}

TEST(Range, BeginEnd)
{
    const Range<long> range(-2, 9);

    EXPECT_EQ(range.begin(), -2);
    EXPECT_EQ(*range.begin(), -2);
    EXPECT_EQ(range.end(), 9);
    EXPECT_EQ(*range.end(), 9);
}

TEST(Range, IsInRangeWorks)
{
    const Range<size_t> range(5, 8);

    EXPECT_EQ(range.isInRange(4), false);
    EXPECT_EQ(range.isInRange(5), true);
    EXPECT_EQ(range.isInRange(6), true);
    EXPECT_EQ(range.isInRange(7), true);
    EXPECT_EQ(range.isInRange(8), false);
}

TEST(Range, IteratorWorks)
{
    const Range<Index> range(-1, 3);

    int minValue = std::numeric_limits<int>::max();
    int maxValue = std::numeric_limits<int>::min();
    for (int i : range)
    {
        minValue = std::min(minValue, i);
        maxValue = std::max(maxValue, i);
    }
    EXPECT_EQ(minValue, -1);
    EXPECT_EQ(maxValue, 2);
}

} // namespace
} // namespace test
} // namespace gmx
