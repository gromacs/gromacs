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
 * Tests for NonbondedOnGpuFromUser
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 */
#include "gmxpre.h"

#include "gromacs/taskassignment/usergpuids.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"

namespace gmx
{

namespace test
{
namespace
{

TEST(GpuIdStringHandlingTest, ParsingAndReconstructionWork)
{
    using ::testing::UnorderedElementsAreArray;

    // TODO It would be nicer to use EXPECT_THAT(assignment,
    // UnorderedElementsAreArray({0,1}) but MSVC 2015 does
    // not deduce the template arguments for <int, 2>.

    // Test simple assignments and back mappings
    {
        const char *strings[] = { "01", "0,1", "0,1," };
        for (const auto &s : strings)
        {
            auto assignment = parseUserGpuIds(s);
            auto matcher    = UnorderedElementsAreArray<int, 2>({0, 1});
            EXPECT_THAT(assignment, matcher) << "for string " << s;
            EXPECT_EQ("0",     makeGpuIdString(assignment, 1));
            EXPECT_EQ("0,1",   makeGpuIdString(assignment, 2));
            EXPECT_EQ("0,0,1", makeGpuIdString(assignment, 3));
        }
    }
    // Test an input that could be a single large index, or two small indices; and back mappings
    {
        auto assignment = parseUserGpuIds("11");
        auto matcher    = UnorderedElementsAreArray<int, 2>({1, 1});
        EXPECT_THAT(assignment, matcher);
        EXPECT_EQ("1",     makeGpuIdString(assignment, 1));
        EXPECT_EQ("1,1",   makeGpuIdString(assignment, 2));
        EXPECT_EQ("1,1,1", makeGpuIdString(assignment, 3));
    }
    // Test an input that must be a single large index; and back mappings
    {
        auto assignment = parseUserGpuIds("11,");
        auto matcher    = UnorderedElementsAreArray<int, 1>({11});
        EXPECT_THAT(assignment, matcher);
        EXPECT_EQ("11",       makeGpuIdString(assignment, 1));
        EXPECT_EQ("11,11",    makeGpuIdString(assignment, 2));
        EXPECT_EQ("11,11,11", makeGpuIdString(assignment, 3));
    }
    // Test multiple large indices; and back mappings
    {
        const char *strings[] = { "11,12", "11,12," };
        for (const auto &s : strings)
        {
            auto assignment = parseUserGpuIds(s);
            auto matcher    = UnorderedElementsAreArray<int, 2>({11, 12});
            EXPECT_THAT(assignment, matcher) << "for string " << s;
            EXPECT_EQ("11",       makeGpuIdString(assignment, 1));
            EXPECT_EQ("11,12",    makeGpuIdString(assignment, 2));
            EXPECT_EQ("11,11,12", makeGpuIdString(assignment, 3));
        }
    }
}

TEST(GpuIdStringHandlingTest, EmptyStringCanBeValid)
{
    using ::testing::IsEmpty;

    auto assignment = parseUserGpuIds("");
    EXPECT_THAT(assignment, IsEmpty());
    EXPECT_EQ("", makeGpuIdString(assignment, 0));
}

TEST(GpuIdStringHandlingTest, InvalidInputsThrow)
{
    {
        const char *strings[] = {
            "a", "0a", ",01", ",0,1", ",0,1,",
            ":0", "0a:1b", "0:1:2",
            ",", ";", ":", "-", "=",
        };
        for (const auto &s : strings)
        {
            EXPECT_THROW(parseUserGpuIds(s), InvalidInputError) << "for string " << s;
        }
    }
}

} // namespace
} // namespace
} // namespace
