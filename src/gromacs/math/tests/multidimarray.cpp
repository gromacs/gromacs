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
 * Tests multidimensional arrays
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/multidimarray.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace
{

class MultiDimArrayTest : public ::testing::Test
{
    public:
        MultiDimArrayTest()
        {
            std::fill(begin(twoDArray_), end(twoDArray_), testNumber_ - 1);
            std::fill(begin(twoDDyanmicArray_), end(twoDDyanmicArray_), testNumber_ - 1);
        }
    protected:
        using container_type = std::vector<float>;

        using extents_type   = extents<3, 3>;
        using array_type     = MultiDimArray < container_type, extents_type>;

        using dynamic_extents_type = extents<dynamic_extent, dynamic_extent>;
        using dynamic_array_type   = MultiDimArray<container_type, dynamic_extents_type>;

        array_type         twoDArray_        = array_type::allocate();
        dynamic_array_type twoDDyanmicArray_ = dynamic_array_type::allocate(3, 3);

        float              testNumber_      = 42;
};


TEST_F(MultiDimArrayTest, canConstructAndFill)
{
    for (const auto &x : twoDArray_)
    {
        EXPECT_EQ(testNumber_ - 1, x);
    }
}

TEST_F(MultiDimArrayTest, canSetValues)
{
    twoDArray_(1, 1) = testNumber_;
    EXPECT_EQ(testNumber_, twoDArray_(1, 1) );

}

TEST_F(MultiDimArrayTest, canMoveConstruct)
{
    auto other(std::move(twoDArray_));
    for (const auto &x : other)
    {
        EXPECT_EQ(testNumber_ - 1, x);
    }
}

TEST_F(MultiDimArrayTest, canMoveAssign)
{
    auto other = array_type::allocate();
    other = std::move(twoDArray_);
    for (const auto &x : other)
    {
        EXPECT_EQ(testNumber_ - 1, x);
    }
}

TEST_F(MultiDimArrayTest, canCopyConstruct)
{
    auto other       = twoDArray_;
    auto twoDArrayIt = begin(twoDArray_);
    for (const auto &x : other)
    {
        EXPECT_EQ(*twoDArrayIt, x);
        ++twoDArrayIt;
    }
}

TEST_F(MultiDimArrayTest, canCopyAssign)
{
    auto other = array_type::allocate();
    other = twoDArray_;
    auto twoDArrayIt = begin(twoDArray_);
    for (const auto &x : other)
    {
        EXPECT_EQ(*twoDArrayIt, x);
        ++twoDArrayIt;
    }
}

TEST_F(MultiDimArrayTest, canSwap)
{
    auto other = array_type::allocate();
    other.swap(twoDArray_);
    for (const auto &x : other)
    {
        EXPECT_EQ(testNumber_ - 1, x);
    }
}

TEST_F(MultiDimArrayTest, staticMultiDimArrayExtent)
{
    EXPECT_EQ(twoDArray_.extent(0), 3);
    EXPECT_EQ(twoDArray_.extent(1), 3);
}

TEST_F(MultiDimArrayTest, dynamicMultiDimArrayResize)
{
    twoDDyanmicArray_.resize(5, 4);
    EXPECT_EQ(twoDDyanmicArray_.extent(0), 5);
    EXPECT_EQ(twoDDyanmicArray_.extent(1), 4);
}

TEST_F(MultiDimArrayTest, dynamicMultiDimArrayResizeAndSetValue)
{
    twoDDyanmicArray_.resize(5, 4);
    twoDDyanmicArray_(4, 3) = testNumber_;
    EXPECT_EQ(twoDDyanmicArray_(4, 3), testNumber_);
}

TEST_F(MultiDimArrayTest, dynamicMultiDimArrayExtent)
{
    EXPECT_EQ(twoDDyanmicArray_.extent(0), 3);
    EXPECT_EQ(twoDDyanmicArray_.extent(1), 3);
}

TEST_F(MultiDimArrayTest, dynamicMultidDimArrayWorks)
{
    twoDDyanmicArray_(1, 2) = testNumber_;
    EXPECT_EQ(testNumber_, twoDDyanmicArray_(1, 2));
}

} // namespace

} // namespace test

} // namespace gmx
