/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Tests for gmx::ExtendableForceField file searching.
 *
 * Uses the static data files in tests/extendable-force-field/, which
 * emulate a base force field directory and two force field addon
 * directories, and checks that files are found in, and reported from,
 * the right places with the base force field taking precedence over the
 * addons, and the addons searched in the order they were provided.
 */
#include "gmxpre.h"

#include "gromacs/gmxpreprocess/extendableforcefield.h"

#include <filesystem>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

using testing::ElementsAre;
using testing::Optional;

/*! \brief Base + two addons with a small, overlapping set of data files:
 *
 * ```
 * base.ff/    common.dat  a.hdb  sub/deep.hdb
 * addon1.ffa/ common.dat  b.hdb
 * addon2.ffa/ c.hdb
 * ```
 */
class ExtendableForceFieldTest : public ::testing::Test
{
public:
    ExtendableForceFieldTest()
    {
        // TestFileManager resolves data file paths relative to this
        // executable's source directory, which contains the test data.
        root_ = TestFileManager::getInputFilePath("extendable-force-field/base.ff/common.dat")
                        .parent_path()  // base.ff
                        .parent_path(); // extendable-force-field
    }

    std::filesystem::path pathIn(const std::string& relativePath) const
    {
        return root_ / relativePath;
    }

    ExtendableForceField forceField() const
    {
        const std::vector<std::filesystem::path> addonDirs{ pathIn("addon1.ffa"), pathIn("addon2.ffa") };
        return ExtendableForceField(pathIn("base.ff"), addonDirs);
    }

    std::filesystem::path root_;
};

TEST_F(ExtendableForceFieldTest, SearchesBaseFirstThenAddonsInOrder)
{
    const ExtendableForceField ff = forceField();

    auto directories = ff.searchDirectories();

    ASSERT_EQ(directories.size(), 3);
    EXPECT_EQ(directories[0], pathIn("base.ff"));
    EXPECT_EQ(directories[1], pathIn("addon1.ffa"));
    EXPECT_EQ(directories[2], pathIn("addon2.ffa"));
}

TEST_F(ExtendableForceFieldTest, SinglePathPrefersTheBaseForceField)
{
    EXPECT_THAT(forceField().findSinglePath("common.dat"), Optional(pathIn("base.ff/common.dat")));
}

TEST_F(ExtendableForceFieldTest, SinglePathFindsAddonFiles)
{
    EXPECT_THAT(forceField().findSinglePath("a.hdb"), Optional(pathIn("base.ff/a.hdb")));
    EXPECT_THAT(forceField().findSinglePath("b.hdb"), Optional(pathIn("addon1.ffa/b.hdb")));
    EXPECT_THAT(forceField().findSinglePath("c.hdb"), Optional(pathIn("addon2.ffa/c.hdb")));
}

TEST_F(ExtendableForceFieldTest, SinglePathReturnsNulloptWhenMissing)
{
    EXPECT_FALSE(forceField().findSinglePath("nothing.dat").has_value());
}

TEST_F(ExtendableForceFieldTest, RequiredPathReturnsTheFirstMatch)
{
    EXPECT_EQ(forceField().findRequiredPath("common.dat"), pathIn("base.ff/common.dat"));
}

TEST_F(ExtendableForceFieldTest, RequiredPathThrowsWhenMissing)
{
    EXPECT_THROW_GMX(forceField().findRequiredPath("nothing.dat"), InvalidInputError);
}

TEST_F(ExtendableForceFieldTest, AllPathsReturnsEveryMatchBaseFirst)
{
    EXPECT_EQ(forceField().findAllPaths("common.dat"),
              (std::vector<std::filesystem::path>{ pathIn("base.ff/common.dat"),
                                                   pathIn("addon1.ffa/common.dat") }));
}

TEST_F(ExtendableForceFieldTest, AllPathsEndingWithFollowsTheSearchOrder)
{
    // The contract is that files are reported with the base force field
    // first, then the addons in the order they were provided, and sorted
    // by name within each directory. Files in subdirectories must not be
    // reported, since the search covers the top level of each directory.
    EXPECT_THAT(forceField().findAllPathsEndingWith(".hdb"),
                ElementsAre(pathIn("base.ff/a.hdb"), pathIn("addon1.ffa/b.hdb"), pathIn("addon2.ffa/c.hdb")));
}

TEST_F(ExtendableForceFieldTest, SwappingAddonOrderSwapsTheirResults)
{
    const std::vector<std::filesystem::path> addonDirs{ pathIn("addon2.ffa"), pathIn("addon1.ffa") };
    const ExtendableForceField ff(pathIn("base.ff"), addonDirs);

    EXPECT_THAT(ff.findAllPathsEndingWith(".hdb"),
                ElementsAre(pathIn("base.ff/a.hdb"), pathIn("addon2.ffa/c.hdb"), pathIn("addon1.ffa/b.hdb")));
    EXPECT_THAT(ff.findSinglePath("common.dat"), Optional(pathIn("base.ff/common.dat")));
}

TEST_F(ExtendableForceFieldTest, AllPathsEndingWithReturnsEmptyWhenNothingMatches)
{
    EXPECT_TRUE(forceField().findAllPathsEndingWith(".zzz").empty());
}

TEST_F(ExtendableForceFieldTest, AllPathsEndingWithThrowsWhenPolicyIsThrowError)
{
    EXPECT_THROW_GMX(forceField().findAllPathsEndingWith(
                             ".zzz", ExtendableForceField::OnEmptySearch::ThrowError),
                     InvalidInputError);
}

TEST_F(ExtendableForceFieldTest, EmptyAddonDirectoriesNeverMatch)
{
    const std::vector<std::filesystem::path> addonDirs{ std::filesystem::path(), pathIn("addon1.ffa") };
    const ExtendableForceField ff(pathIn("base.ff"), addonDirs);

    EXPECT_THAT(ff.findSinglePath("b.hdb"), Optional(pathIn("addon1.ffa/b.hdb")));
    EXPECT_FALSE(ff.findSinglePath("c.hdb").has_value());
}

TEST_F(ExtendableForceFieldTest, ThrowsOnEmptyForceFieldDirectory)
{
    const std::filesystem::path emptyDirectory;

    EXPECT_THROW_GMX(ExtendableForceField{ emptyDirectory }, InvalidInputError);
}

TEST_F(ExtendableForceFieldTest, ThrowsOnMissingForceFieldDirectory)
{
    EXPECT_THROW_GMX(ExtendableForceField{ pathIn("missing.ff") }, InvalidInputError);
}

TEST_F(ExtendableForceFieldTest, ForceFieldNameIsTheBaseDirectory)
{
    const std::vector<std::filesystem::path> addonDirs{ pathIn("addon1.ffa") };
    const ExtendableForceField               ff(pathIn("base.ff"), addonDirs);

    EXPECT_EQ(ff.forceFieldName(), pathIn("base.ff"));
}

} // namespace
} // namespace test
} // namespace gmx
