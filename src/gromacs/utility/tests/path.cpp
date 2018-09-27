/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2018, by the GROMACS development team, led by
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
 * Tests for (some) functions in path.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/path.h"

#include <utility>

#include <gtest/gtest.h>

namespace
{

TEST(PathTest, StripSourcePrefixWorks)
{
    EXPECT_STREQ("", gmx::Path::stripSourcePrefix(""));
    EXPECT_STREQ("foo.cpp", gmx::Path::stripSourcePrefix("foo.cpp"));
    EXPECT_STREQ("foo.cpp", gmx::Path::stripSourcePrefix("some/dir/foo.cpp"));
    EXPECT_STREQ("foo.cpp", gmx::Path::stripSourcePrefix("src/some/dir/foo.cpp"));
    EXPECT_STREQ("foo.cpp",
                 gmx::Path::stripSourcePrefix("srcx/gromacs/foo.cpp"));
    EXPECT_STREQ("src/gromacs/foo.cpp",
                 gmx::Path::stripSourcePrefix("src/gromacs/foo.cpp"));
    EXPECT_STREQ("src/gromacs/foo.cpp",
                 gmx::Path::stripSourcePrefix("some/dir/src/gromacs/foo.cpp"));
    // TODO: For in-source builds, this might not work.
    EXPECT_EQ(gmx::Path::normalize("src/gromacs/utility/tests/path.cpp"),
              gmx::Path::stripSourcePrefix(__FILE__))
    << "stripSourcePrefix() does not work with compiler-produced file names. "
    << "This only affects source paths reported in fatal error messages.";
}

TEST(PathTest, ConcatenateBeforeExtensionWorks)
{
    EXPECT_STREQ("md1.log", gmx::Path::concatenateBeforeExtension("md.log", "1").c_str());
    EXPECT_STREQ("ener0", gmx::Path::concatenateBeforeExtension("ener", "0").c_str());
    EXPECT_STREQ("simpledir/traj34.tng", gmx::Path::concatenateBeforeExtension("simpledir/traj.tng", "34").c_str());
    EXPECT_STREQ("complex.dir/traj34.tng", gmx::Path::concatenateBeforeExtension("complex.dir/traj.tng", "34").c_str());
    EXPECT_STREQ("complex.dir/traj34", gmx::Path::concatenateBeforeExtension("complex.dir/traj", "34").c_str());
}

TEST(PathTest, GetParentPathWorks)
{
    {
        auto result = gmx::Path::getParentPath("");
        EXPECT_EQ("", result);
    }
    {
        auto result = gmx::Path::getParentPath("file.txt");
        EXPECT_EQ("", result);
    }
    {
        auto result = gmx::Path::getParentPath("dir/file.txt");
        EXPECT_EQ("dir", result);
    }
    {
        auto result = gmx::Path::getParentPath("windowsdir\\file.txt");
        EXPECT_EQ("windowsdir", result);
    }
    {
        auto result = gmx::Path::getParentPath("dir/anotherdir/file.txt");
        EXPECT_EQ("dir/anotherdir", result);
    }
    {
        auto result = gmx::Path::getParentPath("dir");
        EXPECT_EQ("", result);
    }
    {
        auto result = gmx::Path::getParentPath("dir/anotherdir");
        EXPECT_EQ("dir", result);
    }
}

TEST(PathTest, GetParentPathAndBasenameWorks)
{
    {
        auto result = gmx::Path::getParentPathAndBasename("");
        EXPECT_EQ("", std::get<0>(result));
        EXPECT_EQ("", std::get<1>(result));
    }
    {
        auto result = gmx::Path::getParentPathAndBasename("file.txt");
        EXPECT_EQ("", std::get<0>(result));
        EXPECT_EQ("file.txt", std::get<1>(result));
    }
    {
        auto result = gmx::Path::getParentPathAndBasename("dir/file.txt");
        EXPECT_EQ("dir", std::get<0>(result));
        EXPECT_EQ("file.txt", std::get<1>(result));
    }
    {
        auto result = gmx::Path::getParentPathAndBasename("windowsdir\\file.txt");
        EXPECT_EQ("windowsdir", std::get<0>(result));
        EXPECT_EQ("file.txt", std::get<1>(result));
    }
    {
        auto result = gmx::Path::getParentPathAndBasename("dir/anotherdir/file.txt");
        EXPECT_EQ("dir/anotherdir", std::get<0>(result));
        EXPECT_EQ("file.txt", std::get<1>(result));
    }
    {
        auto result = gmx::Path::getParentPathAndBasename("dir");
        EXPECT_EQ("", std::get<0>(result));
        EXPECT_EQ("dir", std::get<1>(result));
    }
    {
        auto result = gmx::Path::getParentPathAndBasename("dir/anotherdir");
        EXPECT_EQ("dir", std::get<0>(result));
        EXPECT_EQ("anotherdir", std::get<1>(result));
    }
}

} // namespace
