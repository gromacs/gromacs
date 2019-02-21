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
 * Tests multidimensional indexes
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/multidimindexrange.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testutils/testasserts.h"

namespace gmx
{

TEST(MultiDimIndexRange, MultiDimIteratorConstructs)
{
    MultiDimIndexRange<2>::iterator::value_type beginPosition({0, 0});
    MultiDimIndexRange<2>::iterator::value_type endPosition({2, 2});

    MultiDimIndexRange<2>::iterator             multiDimIterator(beginPosition, endPosition, beginPosition);
    MultiDimIndexRange<2>::iterator::value_type expectedPosition = {0, 0};
    EXPECT_THAT(expectedPosition, testing::Pointwise(testing::Eq(), *multiDimIterator));
}

TEST(MultiDimIndexRange, MultiDimIndexRangeIsAllZero)
{
    MultiDimIndexRange<2>::iterator::value_type beginPosition({0, 0});
    MultiDimIndexRange<2>::iterator::value_type endPosition({2, 2});
    MultiDimIndexRange<2>::iterator             extentsIterator(beginPosition, endPosition, beginPosition);

    MultiDimIndexRange<2>                       range({2, 2});

    EXPECT_EQ(extentsIterator, range.begin());
}

TEST(MultiDimIndexRange, MultiDimIndexRangeEndIsPastLast)
{
    MultiDimIndexRange<2>           range({2, 2});
    MultiDimIndexRange<2>::iterator extentsIterator(*range.begin(), *range.end(), {1, 1});
    ++extentsIterator;
    EXPECT_EQ(extentsIterator, range.end());
}


TEST(MultiDimIndexRange, IteratingOverMultiDim)
{
    MultiDimIndexRange<2> range({2, 2});
    std::vector<typename MultiDimIndexRange<2>::iterator::value_type> expectedPositions = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    auto                  expectedPositionIt = expectedPositions.cbegin();
    for (const auto &positionWithinMultiDim : range)
    {
        EXPECT_THAT(*expectedPositionIt, testing::Pointwise(testing::Eq(), positionWithinMultiDim));
        ++expectedPositionIt;
    }
}

TEST(MultiDimIndexRange, IteratingFromOrigin)
{
    std::vector<typename MultiDimIndexRange<2>::iterator::value_type> expectedPositions = {{1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};
    auto          expectedPositionIt = expectedPositions.cbegin();
    for (const auto &positionWithinMultiDim : MultiDimIndexRange<2>({1, 0}, {3, 3}))
    {
        EXPECT_THAT(*expectedPositionIt, testing::Pointwise(testing::Eq(), positionWithinMultiDim));
        ++expectedPositionIt;
    }
}

TEST(MultiDimIndexRange, IsEqual)
{
    MultiDimIndexRange<2> range({0, 0}, {2, 2});
    MultiDimIndexRange<2> otherRange({0, 0}, {2, 2});
    EXPECT_EQ(range, otherRange);
}

TEST(MultiDimIndexRange, IsNotEqualDifferentBegin)
{
    MultiDimIndexRange<2> range({0, 1}, {2, 2});
    MultiDimIndexRange<2> otherRange({0, 0}, {2, 2});
    EXPECT_NE(range, otherRange);
}

TEST(MultiDimIndexRange, IsNotEqualDifferentEnd)
{
    MultiDimIndexRange<2> range({0, 0}, {2, 3});
    MultiDimIndexRange<2> otherRange({0, 0}, {2, 2});
    EXPECT_NE(range, otherRange);
}

namespace
{
//! Test wether the intersection of two two-dimensional input ranges matches the expected range
void test2DIntersection(MultiDimIndexRange<2> firstRange,
                        MultiDimIndexRange<2> otherRange, MultiDimIndexRange<2> expectedRange)
{
    auto intersectionRange = multiDimIndexRangeIntersection(firstRange, otherRange);
    EXPECT_THAT(*expectedRange.begin(), testing::Pointwise(testing::Eq(), *intersectionRange.begin()));
    EXPECT_THAT(*expectedRange.end(), testing::Pointwise(testing::Eq(), *intersectionRange.end()));
}
}   // namespace

TEST(MultiDimIndexRange, IntersectionWithEncompassingRange)
{
    test2DIntersection({{0, 0}, {2, 2}}, {{-1, -1}, {3, 3}}, {{0, 0}, {2, 2}});
}

TEST(MultiDimIndexRange, IntersectionSameBeginDifferentEnds)
{
    test2DIntersection({{0, 0}, {3, 2}}, {{0, 0}, {2, 3}}, {{0, 0}, {2, 2}});
}

TEST(MultiDimIndexRange, IntersectionSameEndDifferentBegins)
{
    test2DIntersection({{1, 0}, {2, 2}}, {{0, 1}, {2, 2}}, {{1, 1}, {2, 2}});
}

TEST(MultiDimIndexRange, IntersectionDifferentEndDifferentBegins)
{
    test2DIntersection({{1, 0}, {3, 2}}, {{0, 1}, {2, 3}}, {{1, 1}, {2, 2}});
}

TEST(MultiDimIndexRange, IntersectionEmpty)
{
    MultiDimIndexRange<2> firstRange({0, 0}, {2, 2});
    MultiDimIndexRange<2> otherRange({-1, 2}, {4, 4});
    auto                  intersectionRange = multiDimIndexRangeIntersection(firstRange, otherRange);
    // how an empty range is implemented should not be tested
    EXPECT_EQ(intersectionRange.begin(), intersectionRange.end());
}

TEST(MultiDimIndexRange, NullShiftYieldsSameRange)
{
    MultiDimIndexRange<2> range({0, 0}, {2, 2});
    auto                  shiftedRange = multiDimIndexRangeShift(range, {{0, 0}});
    EXPECT_EQ(range, shiftedRange);
}

TEST(MultiDimIndexRange, ShiftOkay)
{
    MultiDimIndexRange<2> range({0, 0}, {2, 2});
    auto                  shiftedRange = multiDimIndexRangeShift(range, {{-10, 10}});

    MultiDimIndexRange<2> expectedRange({-10, 10}, {-8, 12});
    EXPECT_THAT(*expectedRange.begin(), testing::Pointwise(testing::Eq(), *shiftedRange.begin()));
    EXPECT_THAT(*expectedRange.end(), testing::Pointwise(testing::Eq(), *shiftedRange.end()));
}

} // namespace gmx
