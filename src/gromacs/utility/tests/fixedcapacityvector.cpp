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
 * \brief Tests for gmx::FixedCapacityVector
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/fixedcapacityvector.h"

#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

namespace
{

TEST(FixedCapacityVectorTest, IsEmpty)
{
    FixedCapacityVector<int, 2> empty;

    EXPECT_EQ(0U, empty.size());
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
    EXPECT_TRUE(v.empty());
}

} // namespace

} // namespace gmx
