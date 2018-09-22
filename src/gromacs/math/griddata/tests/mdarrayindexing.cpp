/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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
 * Test multidimensional array indexing routines
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/griddata/mdarrayindexing.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace internal
{

namespace
{

TEST(MDArrayOffset, canConstruct)
{
    offset<1> oneD {
        1
    };
    EXPECT_EQ(oneD[0], 1);
    offset<4> fourD {
        4, 2, 78, 1
    };
    EXPECT_EQ(fourD[0], 4);
    EXPECT_EQ(fourD[1], 2);
    EXPECT_EQ(fourD[2], 78);
    EXPECT_EQ(fourD[3], 1);

    std::array<ptrdiff_t, 4> fourDArray = { 4, 2, 78, 1 };
    offset<4>                fourDFromArray {
        fourDArray
    };
    EXPECT_TRUE(fourD == fourDFromArray);

    offset<2> twoDDefault;
    EXPECT_EQ(twoDDefault[0], 0);
    EXPECT_EQ(twoDDefault[1], 0);

}

TEST(MDArrayOffset, defaultConstructorComparesToZero)
{
    offset<2> defaultOffset;
    offset<2> zeroOffset = {0, 0};
    ASSERT_TRUE(defaultOffset == zeroOffset);
}

TEST(MDArrayOffset, canAccessAndCopy)
{
    offset<2> o;
    o[0] = 2;
    o[1] = 3;
    const offset<2> constOffset = o;
    ASSERT_TRUE(constOffset[0] == o[0]);
}

TEST(MDArrayBounds, defaultSizeZero)
{
    bounds<1> oneDBounds;
    bounds<4> fourDBounds;

    ASSERT_EQ(0, oneDBounds.size());
    ASSERT_EQ(0, fourDBounds.size());
}

TEST(MDArrayBounds, containsOffsets)
{
    bounds<1> oneDBounds;
    ASSERT_FALSE(oneDBounds.contains({0}));

    oneDBounds[0] = 5;
    ASSERT_TRUE(oneDBounds.contains({0}));
    ASSERT_FALSE(oneDBounds.contains({5}));

    bounds<4> fourDBounds = { 4, 2, 78, 1 };
    ASSERT_TRUE(fourDBounds.contains({0, 0, 0, 0}));
    ASSERT_TRUE(fourDBounds.contains({3, 1, 77, 0}));
    ASSERT_FALSE(fourDBounds.contains({4, 1, 77, 0}));
    ASSERT_FALSE(fourDBounds.contains({3, 2, 77, 0}));
    ASSERT_FALSE(fourDBounds.contains({3, 1, 78, 0}));
    ASSERT_FALSE(fourDBounds.contains({3, 1, 77, 1}));
}

TEST(MDArrayBounds, sizeCorrect)
{
    ASSERT_EQ(bounds<4>({ 4, 2, 78, 5}).size(), 3120);
}

TEST(MDArrayBounds, beginIteratorDereferencesToZeroOffset)
{
    ASSERT_EQ(*(bounds<4>({ 4, 2, 78, 5}).begin()), offset<4>({0, 0, 0, 0}));
}

TEST(MDArrayBounds, iteratorIteratesCorrectly)
{
    bounds<3>          mdBounds = {2, 5, 3};

    bounds_iterator<3> boundsIterator = {mdBounds, {0, 3, 2}};
    // test iteration order
    // (0,3,2) -> (1,3,2) -> (0,4,2) -> (1,4,2) -> end

    // test pre-fix increment
    ASSERT_EQ(*(++boundsIterator), offset<3>({1, 3, 2}));

    // test post-fix increment
    boundsIterator++;
    ASSERT_EQ(*boundsIterator++, offset<3>({0, 4, 2}));
    ASSERT_EQ(*boundsIterator, offset<3>({1, 4, 2}));

    // test end
    ASSERT_TRUE(++boundsIterator == mdBounds.end());

}

} // namespace

} // namespace internal

} // namespace test

} // namespace gmx
