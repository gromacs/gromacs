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

#include "gromacs/math/griddata/griddata.h"

#include <string>
#include <vector>

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

TEST(GridDataTest, memoryOffsetIsCorrect)
{
    // auto grid = Grid<3>(CanonicalVectorBasis<3>({1, 0.5, 0.25}), {10, 40, 80});
    // GridDataFloat3D griddata(grid);
    // griddata.memoryOffset();
    // std::vector<CanonicalVectorBasis<3>::NdVector> toIndex        = {{{0, 0, 0}}, {{0, 0.49, 0}}, {{0, 10, 0}}, {{0, 0, -2}}, {{-0.12, 0.3, 2.1}}};
    // std::vector<offset<3>> expectedResult = {{0, 0, 0}, {0, 39, 0}, {0, 800, 0}, {0, 0, -640}, {-2, 24, 672}};
    // auto currExpected = std::begin(expectedResult);
    // for (const auto x : toIndex)
    // {
    //     ASSERT_EQ((*currExpected)[0], grid.coordinateToFloorMultiIndex(x)[0]);
    //     ASSERT_EQ((*currExpected)[1], grid.coordinateToFloorMultiIndex(x)[1]);
    //     ASSERT_EQ((*currExpected)[2], grid.coordinateToFloorMultiIndex(x)[2]);
    //     ++currExpected;
    // }
}

} // namespace

} // internal

} // test

} // gmx
