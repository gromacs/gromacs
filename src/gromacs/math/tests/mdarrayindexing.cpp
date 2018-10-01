/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*! \internal
 * \file
 * \brief
 * Test multidimensional array indexing routines
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/mdarrayindexing.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"

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

    std::array<offset<4>::value_type, 4> fourDArray = { 4, 2, 78, 1 };
    offset<4>                            fourDFromArray {
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

    EXPECT_TRUE(defaultOffset == zeroOffset);
}

TEST(MDArrayOffset, canAccessAndCopy)
{
    offset<2> o;
    o[0] = 2;
    o[1] = 3;
    const offset<2> constOffset = o;

    EXPECT_TRUE(constOffset[0] == o[0]);
}

TEST(MDArrayBounds, defaultSizeZero)
{
    bounds<1> oneDBounds;
    bounds<4> fourDBounds;

    EXPECT_EQ(0, oneDBounds.size());
    EXPECT_EQ(0, fourDBounds.size());
}

TEST(MDArrayBounds, containsOffsets)
{
    bounds<1> oneDBounds;
    EXPECT_FALSE(oneDBounds.contains({0}));

    oneDBounds[0] = 5;
    EXPECT_TRUE(oneDBounds.contains({0}));

    EXPECT_FALSE(oneDBounds.contains({5}));
    bounds<4> fourDBounds = { 4, 2, 78, 1 };

    // Extreme points still in grid
    EXPECT_TRUE(fourDBounds.contains({0, 0, 0, 0}));
    EXPECT_TRUE(fourDBounds.contains({3, 1, 77, 0}));
    // one step outside the bounds in each dimension
    EXPECT_FALSE(fourDBounds.contains({4, 1, 77, 0}));
    EXPECT_FALSE(fourDBounds.contains({3, 2, 77, 0}));
    EXPECT_FALSE(fourDBounds.contains({3, 1, 78, 0}));
    EXPECT_FALSE(fourDBounds.contains({3, 1, 77, 1}));
}

TEST(MDArrayBounds, sizeCorrect)
{
    EXPECT_EQ(bounds<4>({ 4, 2, 78, 5}).size(), 3120);
}

TEST(MDArrayBounds, beginIteratorDereferencesToZeroOffset)
{
    EXPECT_EQ(*(bounds<4>({ 4, 2, 78, 5}).begin()), offset<4>({0, 0, 0, 0}));
}

TEST(MDArrayBounds, creationFromIVec)
{
    gmx::IVec offsetIVec {
        10, 2, 4
    };
    offset<DIM> offset(offsetIVec);
    // calls like "offset<2> offset(offsetIVec);" shall not compile

    gmx::IVec boundsIVec {
        10, 2, 4
    };
    bounds<DIM> bounds(boundsIVec);

    EXPECT_FALSE(bounds.contains(offset));
}

TEST(MDArrayBounds, iteratorIteratesCorrectly)
{
    bounds<3>          mdBounds       = {3, 5, 2};
    bounds_iterator<3> boundsIterator = {mdBounds, {2, 3, 0}};
    // test iteration order
    // (2,3,0) -> (2,3,1) -> (2,4,0) -> (2,4,1) -> end

    // test pre-fix increment
    EXPECT_EQ(*(++boundsIterator), offset<3>({2, 3, 1}));

    // test post-fix increment
    boundsIterator++;

    EXPECT_EQ(*boundsIterator++, offset<3>({2, 4, 0}));
    EXPECT_EQ(*boundsIterator, offset<3>({2, 4, 1}));

    // test if end iterator is recognised
    EXPECT_TRUE(++boundsIterator == mdBounds.end());

}

TEST(MDArrayBounds, coversArea)
{
    // 10 point long cube
    const int r        = 10;
    bounds<3> mdBounds = {r, r, r};
    //count grid points in cube that have radius smaller r
    int       pointsWithinSphere = 0;
    for (const auto &gridPoint : mdBounds)
    {
        if (gmx::square(gridPoint[XX])+gmx::square(gridPoint[YY])+gmx::square(gridPoint[ZZ]) < gmx::square(r))
        {
            pointsWithinSphere++;
        }
    }
    // should be ca r^2 * pi/6
    EXPECT_EQ(639, pointsWithinSphere);

}


} // namespace

} // namespace internal

} // namespace test

} // namespace gmx
