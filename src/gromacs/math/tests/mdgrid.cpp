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
/*! \internal \file
 * \brief
 * Tests canonical vector basis
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/mdgrid.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
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

TEST(GridTest, coordinateToFloorOffset)
{
    FloatingPointTolerance    tolerance(defaultRealTolerance());

    MdGrid<3>                 grid({{1, 0.5, 0.25}}, {10, 40, 80});
    std::vector < MdFloatVector < 3>> toIndex = {{{0, 0, 0}}, {{0, 0.49, 0}}, {{0, 10, 0}}, {{0, 0, -2}}, {{-0.12, 0.3, 2.1}}};
    std::vector < offset < 3>> expectedResult = {{0, 0, 0}, {0, 39, 0}, {0, 800, 0}, {0, 0, -640}, {-2, 24, 672}};
    auto currExpected = std::begin(expectedResult);

    for (const auto &x : toIndex)
    {
        EXPECT_FLOAT_EQ_TOL((*currExpected)[0], grid.coordinateToFloorOffset(x)[0], tolerance);
        EXPECT_FLOAT_EQ_TOL((*currExpected)[1], grid.coordinateToFloorOffset(x)[1], tolerance);
        EXPECT_FLOAT_EQ_TOL((*currExpected)[2], grid.coordinateToFloorOffset(x)[2], tolerance);
        ++currExpected;
    }
}

TEST(GridTest, multiIndexToCoordinate)
{
    FloatingPointTolerance    tolerance(defaultRealTolerance());

    MdGrid<3>                 grid({{-1, 0.25, 0.5}}, {10, 40, 80});
    std::vector < offset < 3>> toCoordinate          = {{0, 0, 0}, {0, 39, 0}, {0, 800, 0}, {0, 0, -640}, {-2, 24, 672}};
    std::vector < MdFloatVector < 3>> expectedResult = {{{0, 0, 0}}, {{0, 0.24375001, 0}}, {{0, 5, 0}}, {{0, 0, -4}}, {{0.2, 0.15, 4.2}}};
    auto expected = std::begin(expectedResult);
    for (const auto &i : toCoordinate)
    {
        EXPECT_FLOAT_EQ_TOL((*expected)[0], grid.multiIndexToCoordinate(i)[0], tolerance);
        EXPECT_FLOAT_EQ_TOL((*expected)[1], grid.multiIndexToCoordinate(i)[1], tolerance);
        EXPECT_FLOAT_EQ_TOL((*expected)[2], grid.multiIndexToCoordinate(i)[2], tolerance);
        ++expected;
    }
}

TEST(GridTest, internalVectorFromGridPointToCoordinate)
{
    FloatingPointTolerance    tolerance(defaultRealTolerance());

    MdGrid<1>                 grid1DUnitSpacing({{1}}, {1});
    EXPECT_FLOAT_EQ_TOL(0.1, grid1DUnitSpacing.internalVectorFromGridPointToCoordinate({{0.1}}, {0})[0], tolerance);
    EXPECT_FLOAT_EQ_TOL(1.1, grid1DUnitSpacing.internalVectorFromGridPointToCoordinate({{1.1}}, {0})[0], tolerance);
    EXPECT_FLOAT_EQ_TOL(-1.1, grid1DUnitSpacing.internalVectorFromGridPointToCoordinate({{-1.1}}, {0})[0], tolerance);
    EXPECT_FLOAT_EQ_TOL(-2.1, grid1DUnitSpacing.internalVectorFromGridPointToCoordinate({{-1.1}}, {1})[0], tolerance);

    MdGrid<1> grid1DTenthSpacing({{1}}, {10});
    EXPECT_FLOAT_EQ_TOL(1., grid1DTenthSpacing.internalVectorFromGridPointToCoordinate({{0.1}}, {0})[0], tolerance);
    EXPECT_FLOAT_EQ_TOL(10., grid1DTenthSpacing.internalVectorFromGridPointToCoordinate({{1.1}}, {1})[0], tolerance);
    EXPECT_FLOAT_EQ_TOL(-11., grid1DTenthSpacing.internalVectorFromGridPointToCoordinate({{-1.1}}, {0})[0], tolerance);
    EXPECT_FLOAT_EQ_TOL(-12., grid1DTenthSpacing.internalVectorFromGridPointToCoordinate({{-1.1}}, {1})[0], tolerance);

    MdGrid<3> grid({{-1, 0.25, 0.5}}, {10, 40, 80});
    std::vector < offset < DIM>> gridIndex            = {{0, 0, 0}, {0, 39, 0}, {0, 800, 0}, {0, 0, -640}, {-2, 24, 672}};
    std::vector < MdFloatVector < 3>> spaceCoordinate = {{{0, 0, 0}}, {{0, 0.24375001, 0}}, {{0, 5, 0}}, {{0, 0, -4}}, {{0.2, 0.15, 4.2}}};
    auto currCoordinate = std::begin(spaceCoordinate);
    for (const auto &i : gridIndex)
    {
        EXPECT_FLOAT_EQ_TOL(0, grid.internalVectorFromGridPointToCoordinate(*currCoordinate, i)[0], tolerance);
        EXPECT_FLOAT_EQ_TOL(0, grid.internalVectorFromGridPointToCoordinate(*currCoordinate, i)[1], tolerance);
        EXPECT_FLOAT_EQ_TOL(0, grid.internalVectorFromGridPointToCoordinate(*currCoordinate, i)[2], tolerance);
        ++currCoordinate;
    }
}

TEST(GridTest, copyBehavesAsExpected)
{
    FloatingPointTolerance    tolerance(defaultRealTolerance());
    MdGrid<3>                 grid({{-1, 0.25, 0.5}}, {10, 40, 80});
    const auto                gridDuplicate = grid;
    std::vector < offset < DIM>> gridIndex            = {{0, 0, 0}, {0, 39, 0}, {0, 800, 0}, {0, 0, -640}, {-2, 24, 672}};
    std::vector < MdFloatVector < 3>> spaceCoordinate = {{{0, 0, 0}}, {{0, 0.24375001, 0}}, {{0, 5, 0}}, {{0, 0, -4}}, {{0.2, 0.15, 4.2}}};
    auto currCoordinate = std::begin(spaceCoordinate);
    for (const auto &i : gridIndex)
    {
        EXPECT_FLOAT_EQ_TOL(0, gridDuplicate.internalVectorFromGridPointToCoordinate(*currCoordinate, i)[0], tolerance);
        EXPECT_FLOAT_EQ_TOL(0, gridDuplicate.internalVectorFromGridPointToCoordinate(*currCoordinate, i)[1], tolerance);
        EXPECT_FLOAT_EQ_TOL(0, gridDuplicate.internalVectorFromGridPointToCoordinate(*currCoordinate, i)[2], tolerance);
        ++currCoordinate;
    }
}

} // namespace

} // internal

} // test

} // gmx
