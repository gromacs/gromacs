/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * Tests for the ListOfLists class.
 *
 * \author berk Hess <hess@kth.se>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/listoflists.h"

#include <gtest/gtest.h>

#include "testutils/testasserts.h"

namespace gmx
{

namespace
{

//! Compares two lists element by element
void compareLists(ArrayRef<const char> a, const std::vector<char>& b)
{
    EXPECT_EQ(a.size(), b.size());
    if (a.size() == b.size())
    {
        for (size_t i = 0; i < b.size(); i++)
        {
            EXPECT_EQ(a[i], b[i]);
        }
    }
}

TEST(ListOfLists, EmptyListOfListsWorks)
{
    ListOfLists<char> list;

    EXPECT_EQ(list.size(), 0);
    EXPECT_EQ(list.empty(), true);
    EXPECT_EQ(list.numElements(), 0);
}

TEST(ListOfLists, AppendWorks)
{
    ListOfLists<char> list;

    std::vector<char> v1 = { 5, 3 };
    std::vector<char> v2 = { -1, 7, 4 };
    list.pushBack(v1);
    list.pushBack(v2);
    EXPECT_EQ(list.size(), 2);
    compareLists(list[0], v1);
    compareLists(list[1], v2);
}

TEST(ListOfLists, EmptyListWorks)
{
    ListOfLists<char> list;

    std::vector<char> v = { 5, 3 };
    list.pushBack(v);
    list.pushBack({});
    EXPECT_EQ(list.size(), 2);
    auto a = list[1];
    EXPECT_EQ(a.empty(), true);
}

TEST(ListOfLists, ClearWorks)
{
    ListOfLists<char> list;

    std::vector<char> v = { 5, 3 };
    list.pushBack(v);
    list.pushBack({});
    list.clear();
    EXPECT_EQ(list.empty(), true);
    EXPECT_EQ(list.numElements(), 0);
}

TEST(ListOfLists, OutOfRangeAccessThrows)
{
    ListOfLists<char> list;

    std::vector<char> v = { 5, 3 };
    EXPECT_THROW(list.at(1), std::out_of_range);
}

TEST(ListOfLists, ExtractsAndRestores)
{
    ListOfLists<char> list1;

    std::vector<char> v1 = { 5, 3 };
    std::vector<char> v3 = { -1, 4 };
    list1.pushBack(v1);
    list1.pushBack({});
    list1.pushBack(v3);
    auto             listRanges = list1.listRangesView();
    auto             elements   = list1.elementsView();
    std::vector<int> listRangesVector;
    listRangesVector.insert(listRangesVector.begin(), listRanges.begin(), listRanges.end());
    std::vector<char> elementsVector;
    elementsVector.insert(elementsVector.begin(), elements.begin(), elements.end());
    ListOfLists<char> list2(std::move(listRangesVector), std::move(elementsVector));
    compareLists(list2[0], v1);
    EXPECT_EQ(list2[1].empty(), true);
    compareLists(list2[2], v3);
}

// Test that we can extract raw data from one object and use it correctly generate a new object
TEST(ListOfLists, AppendsListOfLists)
{
    ListOfLists<char> list1;
    ListOfLists<char> list2;

    std::vector<char> v1 = { 5, 3 };
    list1.pushBack(v1);
    std::vector<char> v2 = { 2, -1 };
    std::vector<char> v3 = { 4 };
    list2.pushBack(v2);
    list2.pushBack(v3);
    const char offset = 2;
    list1.appendListOfLists(list2, offset);
    EXPECT_EQ(list1.size(), 3);
    auto a = list1[1];
    EXPECT_EQ(a.size(), 2);
    EXPECT_EQ(a[0], v2[0] + offset);
    EXPECT_EQ(a[1], v2[1] + offset);
}

} // namespace

} // namespace gmx
