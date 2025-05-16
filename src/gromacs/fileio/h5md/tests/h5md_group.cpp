/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Tests for H5MD file group creation, opening and manipulation routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_group.h"

#include <hdf5.h>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(H5mdGroupTest, CreateGroupWorks)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    constexpr char groupName[] = "testGroup";

    {
        SCOPED_TRACE("Create group");
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(file.fileid(), groupName));
    }
    {
        SCOPED_TRACE("Assert that group was created");
        const auto [group, groupGuard] =
                makeH5mdGroupGuard(H5Gopen(file.fileid(), groupName, H5P_DEFAULT));
        ASSERT_GT(H5Iis_valid(group), 0);
    }
}

TEST(H5mdGroupTest, CreateGroupWithEmptyNameThrows)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    EXPECT_THROW(createGroup(file.fileid(), ""), gmx::FileIOError);
}

TEST(H5mdGroupTest, CreateGroupWithDuplicateNameThrows)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    constexpr char groupName[] = "testGroup";

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(file.fileid(), groupName));
    EXPECT_THROW(createGroup(file.fileid(), groupName), gmx::FileIOError);
}

TEST(H5mdGroupTest, CreateGroupForNonexistingContainerThrows)
{
    constexpr char groupName[] = "testGroup";
    EXPECT_THROW(createGroup(H5I_INVALID_HID, groupName), gmx::FileIOError);
}

TEST(H5mdGroupTest, CreateGroupForReadOnlyFileThrows)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    constexpr char groupName[] = "testGroup";

    {
        SCOPED_TRACE("Create file.");
        H5md file(fileName, H5mdFileMode::Write);
    }
    {
        SCOPED_TRACE("Create group in read-only file");
        H5md file(fileName, H5mdFileMode::Read);
        EXPECT_THROW(createGroup(file.fileid(), groupName), gmx::FileIOError);
    }
}

TEST(H5mdGroupTest, OpenGroupWorksForWriteModeFile)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    constexpr char groupName[] = "testGroup";

    {
        SCOPED_TRACE("Create group");
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(file.fileid(), groupName));
    }
    {
        SCOPED_TRACE("Open group");
        const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(file.fileid(), groupName));
        ASSERT_GT(H5Iis_valid(group), 0);
    }
}

TEST(H5mdGroupTest, OpenGroupWorksForReadOnlyModeFile)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    constexpr char groupName[] = "testGroup";

    {
        SCOPED_TRACE("Create group");
        H5md file(fileName, H5mdFileMode::Write);
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(file.fileid(), groupName));
    }
    {
        SCOPED_TRACE("Open group");
        H5md file(fileName, H5mdFileMode::Read);
        const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(file.fileid(), groupName));
        ASSERT_GT(H5Iis_valid(group), 0);
    }
}

TEST(H5mdGroupTest, OpenNonExistingGroupThrows)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    constexpr char groupName[] = "testGroup";

    {
        SCOPED_TRACE("Open group with non-existing name");
        EXPECT_THROW_GMX(openGroup(file.fileid(), groupName), gmx::FileIOError);
    }
    {
        SCOPED_TRACE("Open group with empty name");
        EXPECT_THROW_GMX(openGroup(file.fileid(), ""), gmx::FileIOError);
    }
}

TEST(H5mdGroupTest, OpenGroupInNonexistingContainerThrows)
{
    constexpr char groupName[] = "testGroup";
    EXPECT_THROW_GMX(openGroup(H5I_INVALID_HID, groupName), gmx::FileIOError);
}

} // namespace
} // namespace test
} // namespace gmx
