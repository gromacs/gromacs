/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests canonical vector basis
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/griddata/grid.h"

#include <vector>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testutils/testasserts.h"

// #include "testutils/testfilemanager.h"

namespace gmx
{

namespace test
{

namespace internal
{

namespace
{

TEST(GridTest, canConstructandCopy)
{
    auto g1 = Grid<1>(CanonicalVectorBasis<1>({1}), ColumnMajorLattice<1>({1}));
    auto g2 = g1;
    auto g3(g1);
}

TEST(GridTest, coordinateToFloorMultiIndex)
{
    auto grid = Grid<3>(CanonicalVectorBasis<3>({1, 0.5, 0.25}), ColumnMajorLattice<3>({10, 40, 80}));
    std::vector<CanonicalVectorBasis<3>::NdVector> toIndex        = {{0, 0, 0}, {0, 0.49, 0}, {0, 10, 0}, {0, 0, -2}, {-0.12, 0.3, 2.1}};
    std::vector<ColumnMajorLattice<3>::MultiIndex> expectedResult = {{0, 0, 0}, {0, 39, 0}, {0, 800, 0}, {0, 0, -640}, {-2, 24, 672}};
    auto currExpected = std::begin(expectedResult);
    for (const auto x : toIndex)
    {
        ASSERT_EQ((*currExpected)[0], grid.coordinateToFloorMultiIndex(x)[0]);
        ASSERT_EQ((*currExpected)[1], grid.coordinateToFloorMultiIndex(x)[1]);
        ASSERT_EQ((*currExpected)[2], grid.coordinateToFloorMultiIndex(x)[2]);
        ++currExpected;
    }
}

TEST(GridTest, multiIndexToCoordinate)
{
    auto grid = Grid<3>(CanonicalVectorBasis<3>({-1, 0.25, 0.5}), ColumnMajorLattice<3>({10, 40, 80}));
    std::vector<ColumnMajorLattice<3>::MultiIndex> toCoordinate   = {{0, 0, 0}, {0, 39, 0}, {0, 800, 0}, {0, 0, -640}, {-2, 24, 672}};
    std::vector<CanonicalVectorBasis<3>::NdVector> expectedResult = {{0, 0, 0}, {0, 0.24375001, 0}, {0, 5, 0}, {0, 0, -4}, {0.2, 0.15, 4.2}};
    auto currExpected = std::begin(expectedResult);
    for (const auto i : toCoordinate)
    {
        ASSERT_FLOAT_EQ((*currExpected)[0], grid.multiIndexToCoordinate(i)[0]);
        ASSERT_FLOAT_EQ((*currExpected)[1], grid.multiIndexToCoordinate(i)[1]);
        ASSERT_FLOAT_EQ((*currExpected)[2], grid.multiIndexToCoordinate(i)[2]);
        ++currExpected;
    }
}

TEST(GridTest, gridVectorFromGridPointToCoordinate)
{
    auto grid1DUnitSpacing = Grid<1>(CanonicalVectorBasis<1>({1}), ColumnMajorLattice<1>({1}));
    ASSERT_FLOAT_EQ(0.1, grid1DUnitSpacing.gridVectorFromGridPointToCoordinate({0.1}, {0})[0]);
    ASSERT_FLOAT_EQ(1.1, grid1DUnitSpacing.gridVectorFromGridPointToCoordinate({1.1}, {0})[0]);
    ASSERT_FLOAT_EQ(-1.1, grid1DUnitSpacing.gridVectorFromGridPointToCoordinate({-1.1}, {0})[0]);
    ASSERT_FLOAT_EQ(-2.1, grid1DUnitSpacing.gridVectorFromGridPointToCoordinate({-1.1}, {1})[0]);

    auto grid1DTenthSpacing = Grid<1>(CanonicalVectorBasis<1>({1}), ColumnMajorLattice<1>({10}));
    ASSERT_FLOAT_EQ(1., grid1DTenthSpacing.gridVectorFromGridPointToCoordinate({0.1}, {0})[0]);
    ASSERT_FLOAT_EQ(10., grid1DTenthSpacing.gridVectorFromGridPointToCoordinate({1.1}, {1})[0]);
    ASSERT_FLOAT_EQ(-11., grid1DTenthSpacing.gridVectorFromGridPointToCoordinate({-1.1}, {0})[0]);
    ASSERT_FLOAT_EQ(-12., grid1DTenthSpacing.gridVectorFromGridPointToCoordinate({-1.1}, {1})[0]);

    auto grid = Grid<3>(CanonicalVectorBasis<3>({-1, 0.25, 0.5}), ColumnMajorLattice<3>({10, 40, 80}));
    std::vector<ColumnMajorLattice<3>::MultiIndex> gridIndex       = {{0, 0, 0}, {0, 39, 0}, {0, 800, 0}, {0, 0, -640}, {-2, 24, 672}};
    std::vector<CanonicalVectorBasis<3>::NdVector> spaceCoordinate = {{0, 0, 0}, {0, 0.24375001, 0}, {0, 5, 0}, {0, 0, -4}, {0.2, 0.15, 4.2}};
    auto currCoordinate = std::begin(spaceCoordinate);
    for (const auto i : gridIndex)
    {
        ASSERT_FLOAT_EQ(0, grid.gridVectorFromGridPointToCoordinate(*currCoordinate, i)[0]);
        ASSERT_FLOAT_EQ(0, grid.gridVectorFromGridPointToCoordinate(*currCoordinate, i)[1]);
        ASSERT_FLOAT_EQ(0, grid.gridVectorFromGridPointToCoordinate(*currCoordinate, i)[2]);
        ++currCoordinate;
    }
}

TEST(GridTest, duplicateBehavesAsExpected)
{
    auto grid          = Grid<3>(CanonicalVectorBasis<3>({-1, 0.25, 0.5}), ColumnMajorLattice<3>({10, 40, 80}));
    auto gridDuplicate = grid.duplicate();
    std::vector<ColumnMajorLattice<3>::MultiIndex> gridIndex       = {{0, 0, 0}, {0, 39, 0}, {0, 800, 0}, {0, 0, -640}, {-2, 24, 672}};
    std::vector<CanonicalVectorBasis<3>::NdVector> spaceCoordinate = {{0, 0, 0}, {0, 0.24375001, 0}, {0, 5, 0}, {0, 0, -4}, {0.2, 0.15, 4.2}};
    auto currCoordinate = std::begin(spaceCoordinate);
    for (const auto i : gridIndex)
    {
        ASSERT_FLOAT_EQ(0, gridDuplicate->gridVectorFromGridPointToCoordinate(*currCoordinate, i)[0]);
        ASSERT_FLOAT_EQ(0, gridDuplicate->gridVectorFromGridPointToCoordinate(*currCoordinate, i)[1]);
        ASSERT_FLOAT_EQ(0, gridDuplicate->gridVectorFromGridPointToCoordinate(*currCoordinate, i)[2]);
        ++currCoordinate;
    }
}

} // namespace

} // internal

} // test

} // gmx
