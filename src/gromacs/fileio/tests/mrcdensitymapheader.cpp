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
 * Tests for mrc file header structure.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/mrcdensitymapheader.h"

#include <cstddef>

#include <array>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/math/coordinatetransformation.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(MrcDensityMapHeaderTest, DataSizeIsZeroForDefaultHeader)
{
    MrcDensityMapHeader header;
    EXPECT_EQ(0, numberOfExpectedDataItems(header));
}

TEST(MrcDensityMapHeaderTest, DataSizeIsCorrect)
{
    MrcDensityMapHeader header;
    header.numColumnRowSection_ = { 1, 2, 3 };
    EXPECT_EQ(6, numberOfExpectedDataItems(header));
}

TEST(MrcDensityMapHeaderTest, DataSizeThrowsWhenInvalid)
{
    MrcDensityMapHeader header;
    header.numColumnRowSection_ = { -1, 2, 3 };
    EXPECT_THROW(numberOfExpectedDataItems(header), InternalError);
}
TEST(MrcDensityMapHeaderTest, GetsCorrectCoordinateTransformNoOriginGiven)
{
    MrcDensityMapHeader header;

    header.extent_                = { 100, 200, 300 };
    header.columnRowSectionStart_ = { 50, 200, 0 };
    header.cellLength_            = { 10, 20, 15 };

    std::array<RVec, 2> testVectors = { RVec{ 0., 0., 0. }, RVec{ 1., 1., 1. } };
    getCoordinateTransformationToLattice(header)(testVectors);

    std::vector<RVec> expectedVectors = { { -50, -200, 0 }, { 50, -100, 200 } };
    EXPECT_THAT(expectedVectors,
                testing::Pointwise(test::RVecEq(test::defaultFloatTolerance()), testVectors));
}

TEST(MrcDensityMapHeaderTest, GetsCorrectCoordinateTransformWithOriginDefined)
{
    MrcDensityMapHeader header;
    header.userDefinedFloat_[12] = 1.;
    header.userDefinedFloat_[13] = 2.;
    header.userDefinedFloat_[14] = 3.;
    header.extent_               = { 100, 200, 300 };
    // setting the columnRowSectionStart values that are to be ignored if userDefinedFloat_ is not zero
    header.columnRowSectionStart_ = { 50, 200, 0 };
    header.cellLength_            = { 10, 20, 15 };

    std::array<RVec, 2> testVectors = { RVec{ 0., 0., 0. }, RVec{ 1., 1., 1. } };
    getCoordinateTransformationToLattice(header)(testVectors);

    std::vector<RVec> expectedVectors = { { -10, -20, -60 }, { 90, 80, 140 } };
    EXPECT_THAT(expectedVectors,
                testing::Pointwise(test::RVecEq(test::defaultFloatTolerance()), testVectors));
}

TEST(MrcDensityMapHeaderTest, GetsCorrectCoordinateTransformWithStartValues)
{
    MrcDensityMapHeader header;
    header.userDefinedFloat_[12]  = 0;
    header.userDefinedFloat_[13]  = 0;
    header.userDefinedFloat_[14]  = 0;
    header.extent_                = { 100, 200, 300 };
    header.columnRowSectionStart_ = { 50, 200, 0 };
    header.cellLength_            = { 10, 20, 15 };

    std::array<RVec, 2> testVectors = { RVec{ 0., 0., 0. }, RVec{ 1., 1., 1. } };
    getCoordinateTransformationToLattice(header)(testVectors);

    std::vector<RVec> expectedVectors = { { -50, -200, 0 }, { 50, -100, 200 } };
    EXPECT_THAT(expectedVectors,
                testing::Pointwise(test::RVecEq(test::defaultFloatTolerance()), testVectors));
}

TEST(MrcDensityMapHeaderTest, GetsCorrectExtents)
{
    MrcDensityMapHeader header;
    header.numColumnRowSection_ = { 100, 200, 300 };

    const auto                      extents         = getDynamicExtents3D(header);
    std::array<std::ptrdiff_t, DIM> expectedExtents = { 300, 200, 100 };
    EXPECT_EQ(expectedExtents[XX], extents.extent(XX));
    EXPECT_EQ(expectedExtents[YY], extents.extent(YY));
    EXPECT_EQ(expectedExtents[ZZ], extents.extent(ZZ));
}

TEST(MrcDensityMapHeaderTest, IsSane)
{
    MrcDensityMapHeader header;
    EXPECT_TRUE(mrcHeaderIsSane(header));

    header.numColumnRowSection_[YY] = -1;
    EXPECT_FALSE(mrcHeaderIsSane(header));

    header.numColumnRowSection_[YY] = 10'000'000;
    EXPECT_FALSE(mrcHeaderIsSane(header));

    header = {};
    EXPECT_TRUE(mrcHeaderIsSane(header));

    header.cellAngles_[XX] = -20;
    EXPECT_FALSE(mrcHeaderIsSane(header));
}

} // namespace
} // namespace test
} // namespace gmx
