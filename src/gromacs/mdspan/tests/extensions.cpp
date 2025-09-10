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
 * \brief
 * Tests gromacs extensions to mdspan proposal
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_mdspan
 */
#include "gmxpre.h"

#include "gromacs/mdspan/extensions.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/mdspan/mdspan.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(MdSpanExtension, SlicingAllStatic)
{
    std::array<int, 2 * 3>           data = { 1, 2, 3, 4, 5, 6 };
    basic_mdspan<int, extents<3, 2>> span{ data.data() };
    EXPECT_EQ(span[0][0], 1);
    EXPECT_EQ(span[0][1], 2);

    EXPECT_EQ(span[1][0], 3);
    EXPECT_EQ(span[1][1], 4);

    EXPECT_EQ(span[2][0], 5);
    EXPECT_EQ(span[2][1], 6);
}

TEST(MdSpanExtension, SlicingDynamic)
{
    std::array<int, 2 * 3>                                     data = { 1, 2, 3, 4, 5, 6 };
    basic_mdspan<int, extents<dynamic_extent, dynamic_extent>> span{ data.data(), 2, 2 };

    EXPECT_EQ(span[0][0], 1);
    EXPECT_EQ(span[0][1], 2);

    EXPECT_EQ(span[1][0], 3);
    EXPECT_EQ(span[1][1], 4);

    EXPECT_EQ(span[2][0], 5);
    EXPECT_EQ(span[2][1], 6);
}

TEST(MdSpanExtension, SlicingAllStatic3D)
{
    std::array<int, 2 * 2 * 2>          data = { 1, 2, 3, 4, 5, 6, 7, 8 };
    basic_mdspan<int, extents<2, 2, 2>> span{ data.data() };
    EXPECT_EQ(span[0][0][0], 1);
    EXPECT_EQ(span[0][0][1], 2);

    EXPECT_EQ(span[0][1][0], 3);
    EXPECT_EQ(span[0][1][1], 4);

    EXPECT_EQ(span[1][0][0], 5);
    EXPECT_EQ(span[1][0][1], 6);

    EXPECT_EQ(span[1][1][0], 7);
    EXPECT_EQ(span[1][1][1], 8);
}

TEST(MdSpanExtension, SlicingEqualsView3D)
{
    std::array<int, 2 * 2 * 2>          data = { 1, 2, 3, 4, 5, 6, 7, 8 };
    basic_mdspan<int, extents<2, 2, 2>> span{ data.data() };
    EXPECT_EQ(span[0][0][0], span(0, 0, 0));
    EXPECT_EQ(span[0][0][1], span(0, 0, 1));
    EXPECT_EQ(span[0][1][0], span(0, 1, 0));
    EXPECT_EQ(span[0][1][1], span(0, 1, 1));
    EXPECT_EQ(span[1][0][0], span(1, 0, 0));
    EXPECT_EQ(span[1][0][1], span(1, 0, 1));
    EXPECT_EQ(span[1][1][0], span(1, 1, 0));
    EXPECT_EQ(span[1][1][1], span(1, 1, 1));
}

TEST(MdSpanExtension, additionWorks)
{
    std::array<int, 2 * 2 * 2>          arr1 = { { -4, -3, -2, -1, 0, 1, 2, 3 } };
    std::array<int, 2 * 2 * 2>          arr2 = { { 1, 1, 1, 1, 1, 1, 1, 1 } };
    basic_mdspan<int, extents<2, 2, 2>> span1{ arr1.data() };
    basic_mdspan<int, extents<2, 2, 2>> span2{ arr2.data() };

    auto result = addElementwise(span1, span2);
    EXPECT_EQ(result[0][0][1], -2);
    EXPECT_EQ(result[0][1][0], -1);
    EXPECT_EQ(result[1][0][0], 1);
}

TEST(MdSpanExtension, subtractionWorks)
{
    std::array<int, 2 * 2 * 2>          arr1 = { { -4, -3, -2, -1, 0, 1, 2, 3 } };
    std::array<int, 2 * 2 * 2>          arr2 = { { 1, 1, 1, 1, 1, 1, 1, 1 } };
    basic_mdspan<int, extents<2, 2, 2>> span1{ arr1.data() };
    basic_mdspan<int, extents<2, 2, 2>> span2{ arr2.data() };

    auto result = subtractElementwise(span1, span2);
    EXPECT_EQ(result[0][0][1], -4);
    EXPECT_EQ(result[0][1][0], -3);
    EXPECT_EQ(result[1][0][0], -1);
}

TEST(MdSpanExtension, multiplicationWorks)
{
    std::array<int, 2 * 2 * 2>          arr1 = { { -4, -3, -2, -1, 0, 1, 2, 3 } };
    std::array<int, 2 * 2 * 2>          arr2 = { {
            2,
            2,
            2,
            2,
            2,
            2,
            2,
    } };
    basic_mdspan<int, extents<2, 2, 2>> span1{ arr1.data() };
    basic_mdspan<int, extents<2, 2, 2>> span2{ arr2.data() };

    auto result = multiplyElementwise(span1, span2);
    EXPECT_EQ(result[0][0][1], -6);
    EXPECT_EQ(result[0][1][0], -4);
    EXPECT_EQ(result[1][0][0], 0);
}

TEST(MdSpanExtension, divisionWorks)
{
    std::array<float, 2 * 2 * 2>          arr1 = { { -4, -3, -2, -1, 0, 1, 2, 3 } };
    std::array<float, 2 * 2 * 2>          arr2 = { {
            2,
            2,
            2,
            2,
            2,
            2,
            2,
            2,
    } };
    basic_mdspan<float, extents<2, 2, 2>> span1{ arr1.data() };
    basic_mdspan<float, extents<2, 2, 2>> span2{ arr2.data() };

    auto result = divideElementwise(span1, span2);
    EXPECT_EQ(result[0][0][1], -1.5);
    EXPECT_EQ(result[0][1][0], -1);
    EXPECT_EQ(result[1][0][0], 0);
}

} // namespace
} // namespace test
} // namespace gmx
