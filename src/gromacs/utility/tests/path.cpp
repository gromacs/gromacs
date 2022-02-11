/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * Tests for (some) functions in path.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/path.h"

#include <utility>

#include <gtest/gtest.h>

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(PathTest, StripSourcePrefixWorks)
{
    EXPECT_STREQ("", Path::stripSourcePrefix(""));
    EXPECT_STREQ("foo.cpp", Path::stripSourcePrefix("foo.cpp"));
    EXPECT_STREQ("foo.cpp", Path::stripSourcePrefix("some/dir/foo.cpp"));
    EXPECT_STREQ("foo.cpp", Path::stripSourcePrefix("src/some/dir/foo.cpp"));
    EXPECT_STREQ("foo.cpp", Path::stripSourcePrefix("srcx/gromacs/foo.cpp"));
    EXPECT_STREQ("src/gromacs/foo.cpp", Path::stripSourcePrefix("src/gromacs/foo.cpp"));
    EXPECT_STREQ("src/gromacs/foo.cpp", Path::stripSourcePrefix("some/dir/src/gromacs/foo.cpp"));
    // TODO: For in-source builds, this might not work.
    EXPECT_EQ(Path::normalize("src/gromacs/utility/tests/path.cpp"), Path::stripSourcePrefix(__FILE__))
            << "stripSourcePrefix() does not work with compiler-produced file names. "
            << "This only affects source paths reported in fatal error messages.";
}

class PathSearchTest : public testing::TestWithParam<std::string>
{
};

TEST_P(PathSearchTest, SearchOperationsWork)
{
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker rootChecker(data.rootChecker());
    std::string                     input = GetParam();

    auto checker = rootChecker.checkCompound("PathToTest", input);
    {
        std::string result;
        ASSERT_NO_THROW_GMX(result = Path::getParentPath(input));
        checker.checkString(result, "getParentPath");
    }
    {
        std::string result;
        ASSERT_NO_THROW_GMX(result = Path::getFilename(input));
        checker.checkString(result, "getFilename");
    }
    {
        bool result = false;
        ASSERT_NO_THROW_GMX(result = Path::hasExtension(input));
        checker.checkBoolean(result, "hasExtension");
    }
    {
        bool result = false;
        ASSERT_NO_THROW_GMX(result = Path::extensionMatches(input, "pdb"));
        checker.checkBoolean(result, "extensionMatchesPdb");
        // The match is exclusive of the dot separator, so no
        // input string can match.
        ASSERT_FALSE(Path::extensionMatches(input, ".pdb"));
    }
    {
        std::string result;
        ASSERT_NO_THROW_GMX(result = Path::stripExtension(input));
        checker.checkString(result, "stripExtension");
    }
    {
        std::string result;
        ASSERT_NO_THROW_GMX(result = Path::concatenateBeforeExtension(input, "_34"));
        checker.checkString(result, "concatenateBeforeExtension");
    }
}

INSTANTIATE_TEST_SUITE_P(WithInputPaths,
                         PathSearchTest,
                         testing::Values("",
                                         "md.log",
                                         "md",
                                         "/tmp/absolute.txt",
                                         "simpledir/traj.tng",
                                         "simpledir/traj",
                                         "windowsdir\\traj.tng",
                                         "complex.dir/traj.tng",
                                         "complex.dir/traj",
                                         "nested/dir/conf.pdb",
                                         "/tmp/absolutedir/conf.pdb"));

} // namespace
} // namespace test
} // namespace gmx
