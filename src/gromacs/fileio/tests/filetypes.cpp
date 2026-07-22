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
    { 12, ".gro" }, { 13, ".g96" }, { 14, ".pdb" }, { 15, ".cif" }, { 16, ".brk" }, { 17, ".ent" },
    { 18, ".esp" }, { 19, ".pqr" }, { 20, ".cpt" }, { 21, ".log" }, { 22, ".xvg" }, { 23, ".out" },
    { 24, ".ndx" }, { 25, ".top" }, { 26, ".itp" }, { 28, ".tpr" }, { 29, ".tex" }, { 30, ".rtp" },
    { 31, ".atp" }, { 32, ".hdb" }, { 33, ".dat" }, { 34, ".dlg" }, { 35, ".map" }, { 36, ".eps" },
    { 37, ".mat" }, { 38, ".m2p" }, { 39, ".mtx" }, { 40, ".edi" }, { 41, ".cub" }, { 42, ".xpm" },
    { 44, ".csv" }, { 45, ".inp" }
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
