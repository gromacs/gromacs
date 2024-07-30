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
 * \brief Tests for gmx::FixedCapacityVector
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/fixedcapacityvector.h"

#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(FixedCapacityVectorTest, IsEmpty)
{
    FixedCapacityVector<int, 2> empty;

    EXPECT_EQ(0U, empty.size());
    EXPECT_EQ(2U, empty.capacity());
    EXPECT_EQ(2U, empty.max_size());
    EXPECT_TRUE(empty.empty());
}

TEST(FixedCapacityVectorTest, PushWorks)
{
    FixedCapacityVector<int, 2> v;

    v.push_back(8);
    EXPECT_EQ(1U, v.size());
    EXPECT_EQ(8, v[0]);

    v.push_back(5);
    EXPECT_EQ(2U, v.size());
    EXPECT_EQ(8, v[0]);
    EXPECT_EQ(5, v[1]);
}

TEST(FixedCapacityVectorTest, PopWorks)
{
    FixedCapacityVector<int, 3> v;

    v.push_back(8);
    EXPECT_EQ(1U, v.size());

    v.pop_back();
    EXPECT_EQ(0U, v.size());
    EXPECT_TRUE(v.empty());
}

TEST(FixedCapacityVectorTest, ClearWorks)
{
    FixedCapacityVector<int, 3> v;

    v.push_back(8);
    EXPECT_EQ(1U, v.size());

    v.clear();
    EXPECT_EQ(0U, v.size());
    EXPECT_TRUE(v.empty());
}

TEST(FixedCapacityVectorTest, EmplaceBackWorks)
{
    FixedCapacityVector<std::vector<int>, 2> v;

    const auto& elem = v.emplace_back(5);
    EXPECT_EQ(1U, v.size());
    EXPECT_EQ(5U, elem.size());
}

TEST(FixedCapacityVectorTest, AtThrows)
{
    FixedCapacityVector<int, 3> v;

    v.push_back(8);
    EXPECT_EQ(8, v.at(0));

    EXPECT_THROW(v.at(1), std::out_of_range);
}

TEST(FixedCapacityVectorTest, IteratorWorks)
{
    FixedCapacityVector<int, 4> v;

    std::vector<int> ref = { 7, 4, 5 };

    for (auto elem : ref)
    {
        v.push_back(elem);
    }
    EXPECT_EQ(3U, v.size());

    std::vector<int> loopResult;
    for (auto elem : v)
    {
        loopResult.push_back(elem);
    }
    EXPECT_EQ(3U, loopResult.size());
    EXPECT_EQ(ref[0], loopResult[0]);
    EXPECT_EQ(ref[1], loopResult[1]);
    EXPECT_EQ(ref[2], loopResult[2]);
}

TEST(FixedCapacityVectorTest, ReverseIteratorWorks)
{
    FixedCapacityVector<int, 4> v;

    std::vector<int> ref = { 7, 4, 5 };

    for (auto elem : ref)
    {
        v.push_back(elem);
    }
    EXPECT_EQ(3U, v.size());

    std::vector<int> loopResult;
    for (auto it = v.rbegin(); it != v.rend(); ++it)
    {
        loopResult.push_back(*it);
    }
    EXPECT_EQ(3U, loopResult.size());
    EXPECT_EQ(ref[0], loopResult[2]);
    EXPECT_EQ(ref[1], loopResult[1]);
    EXPECT_EQ(ref[2], loopResult[0]);
}

TEST(FixedCapacityVectorTest, ZeroCapacityWorks)
{
    FixedCapacityVector<int, 0> v;

    EXPECT_EQ(0U, v.size());
    EXPECT_EQ(0U, v.capacity());
    EXPECT_EQ(0U, v.max_size());
    EXPECT_TRUE(v.empty());
}

TEST(FixedCapacityVectorTest, CopyConstructorWorks)
{
    FixedCapacityVector<int, 4> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    FixedCapacityVector<int, 4> copy(v);
    EXPECT_EQ(1, copy[0]);
    EXPECT_EQ(2, copy[1]);
    EXPECT_EQ(3, copy[2]);
    EXPECT_EQ(3, copy.size());
}

TEST(FixedCapacityVectorTest, CopyAssignmentWorks)
{
    FixedCapacityVector<int, 4> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    FixedCapacityVector<int, 4> copy;
    copy = v;
    EXPECT_EQ(1, copy[0]);
    EXPECT_EQ(2, copy[1]);
    EXPECT_EQ(3, copy[2]);
    EXPECT_EQ(3, copy.size());
}

TEST(FixedCapacityVectorTest, MoveConstructorWorks)
{
    FixedCapacityVector<int, 3> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    // We are testing for correctness and don't care about performance
    // NOLINTNEXTLINE(performance-move-const-arg)
    FixedCapacityVector<int, 3> copy(std::move(v));
    EXPECT_EQ(1, copy[0]);
    EXPECT_EQ(2, copy[1]);
    EXPECT_EQ(3, copy[2]);
    EXPECT_EQ(3, copy.size());
}

TEST(FixedCapacityVectorTest, MoveAssignmentWorks)
{
    FixedCapacityVector<int, 3> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    FixedCapacityVector<int, 3> copy;
    // We are testing for correctness and don't care about performance
    // NOLINTNEXTLINE(performance-move-const-arg)
    copy = std::move(v);
    EXPECT_EQ(1, copy[0]);
    EXPECT_EQ(2, copy[1]);
    EXPECT_EQ(3, copy[2]);
    EXPECT_EQ(3, copy.size());
}

TEST(FixedCapacityVectorTest, ElementAssignmentWorks)
{
    FixedCapacityVector<int, 4> v;
    v.push_back(1);
    v.push_back(2);
    v[0]    = 11;
    v.at(1) = 12;
    EXPECT_EQ(11, v[0]);
    EXPECT_EQ(12, v[1]);
}

TEST(FixedCapacityVectorTest, DataWorks)
{
    FixedCapacityVector<int, 4> v;
    v.push_back(1);
    v.push_back(2);
    EXPECT_EQ(1, v.data()[0]);
    EXPECT_EQ(2, v.data()[1]);
}

TEST(FixedCapacityVectorTest, ConstMethodsWork)
{
    FixedCapacityVector<int, 4> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    const FixedCapacityVector<int, 4> copy = v;
    EXPECT_EQ(3, copy.size());
    EXPECT_EQ(3, copy.ssize());
    EXPECT_FALSE(copy.empty());
    static_assert(copy.capacity() == 4);
    static_assert(copy.max_size() == 4);
    EXPECT_EQ(1, copy[0]);
    EXPECT_EQ(2, copy.at(1));
    EXPECT_EQ(3, copy.data()[2]);
    for (const auto* it = copy.begin(); it != copy.end(); it++)
    {
        EXPECT_EQ(*it, (it - copy.begin() + 1));
    }
}

} // namespace
} // namespace test
} // namespace gmx
