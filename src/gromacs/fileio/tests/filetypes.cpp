/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * Tests for routines for file type handling.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/filetypes.h"

#include <filesystem>
#include <string>
#include <tuple>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace test
{
namespace
{

using TypeAndName = std::tuple<int, std::string>;

using FileTypeTestParams = std::tuple<TypeAndName, std::string>;

class FileTypeTest : public ::testing::Test, public ::testing::WithParamInterface<FileTypeTestParams>
{
public:
    static void runTest(const TypeAndName& params);
};

void FileTypeTest::runTest(const TypeAndName& params)
{
    const int   type = std::get<0>(params);
    const auto& path = std::get<1>(params);
    EXPECT_EQ(type, fn2ftp(path));
    // also test
}

TEST_P(FileTypeTest, CorrectValueForNullptr)
{
    // A lot of places in the code still call fn2ftp(nullptr).
    ASSERT_EQ(fn2ftp(nullptr), efNR);
}

TEST_F(FileTypeTest, CorrectValueForEmptyString)
{
    runTest({ efNR, "" });
}

TEST_F(FileTypeTest, CorrectValueForNoExtension)
{
    runTest({ efNR, "test" });
}

TEST_F(FileTypeTest, CorrectValueForEmptyExtension)
{
    runTest({ efNR, "test." });
}

TEST_F(FileTypeTest, CorrectValueForLongExtensionWithStrangeCharacters)
{
    runTest({ efNR, "test.1234\\abc" });
}

TEST_P(FileTypeTest, CorrectValueForExtension)
{
    auto param       = GetParam();
    auto typeAndName = std::get<0>(param);
    auto prefix      = std::get<1>(param);
    auto fullName    = prefix + std::get<1>(typeAndName);
    runTest({ std::get<0>(typeAndName), fullName });
}

const std::vector<TypeAndName> testParams = {
    { 0, ".mdp" },  { 4, ".trr" },  { 6, ".xtc" },  { 7, ".tng" },  { 8, ".h5md" }, { 9, ".edr" },
    { 12, ".gro" }, { 13, ".g96" }, { 14, ".pdb" }, { 15, ".brk" }, { 16, ".ent" }, { 17, ".esp" },
    { 18, ".pqr" }, { 19, ".cpt" }, { 20, ".log" }, { 21, ".xvg" }, { 22, ".out" }, { 23, ".ndx" },
    { 24, ".top" }, { 25, ".itp" }, { 27, ".tpr" }, { 28, ".tex" }, { 29, ".rtp" }, { 30, ".atp" },
    { 31, ".hdb" }, { 32, ".dat" }, { 33, ".dlg" }, { 34, ".map" }, { 35, ".eps" }, { 36, ".mat" },
    { 37, ".m2p" }, { 38, ".mtx" }, { 39, ".edi" }, { 40, ".cub" }, { 41, ".xpm" }, { 43, ".csv" },
    { 44, ".inp" }
};

const std::vector<std::string> prefixes = { "",
                                            "test",
                                            "test.pdb",
                                            "a/../b/test",
                                            "james.gro/system.mdp/test.pdb" };

INSTANTIATE_TEST_SUITE_P(FileTypeMatch,
                         FileTypeTest,
                         ::testing::Combine(::testing::ValuesIn(testParams), ::testing::ValuesIn(prefixes)));


} // namespace
} // namespace test
} // namespace gmx
