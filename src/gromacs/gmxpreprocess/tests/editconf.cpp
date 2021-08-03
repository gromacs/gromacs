/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2021, by the GROMACS development team, led by
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
 * Tests for editconf
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/editconf.h"

#include <tuple>

#include "gromacs/fileio/filetypes.h"

#include "testutils/cmdlinetest.h"
#include "testutils/filematchers.h"
#include "testutils/refdata.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

using test::CommandLine;
using test::TextFileMatch;

//! Test parameter struct.
using CommandLineOptionParams = std::tuple<std::string, int>;

class EditconfTest :
    public test::CommandLineTestBase,
    public ::testing::WithParamInterface<CommandLineOptionParams>
{
public:
    EditconfTest()
    {
        const auto& params = GetParam();
        setInputFile("-f", std::get<0>(params));

        std::string outputfile = "output.";
        outputfile += ftp2ext(std::get<1>(params));
        ExactTextMatch settings;
        setOutputFile("-o", outputfile.c_str(), TextFileMatch(settings));
    }

    void runTest(const char* testName)
    {
        // Get the command line flags that were set up in the constructor
        CommandLine& cmdline = commandLine();

        // Provide the name of the module to call
        std::string module[] = { "editconf" };
        cmdline.merge(CommandLine(module));

        // Call the module
        ASSERT_EQ(0, gmx_editconf(cmdline.argc(), cmdline.argv()));

        // Check the output
        const auto*          extension = ftp2ext(std::get<1>(GetParam()));
        TestReferenceChecker rootChecker(this->rootChecker());
        rootChecker.checkString(extension, testName);
        checkOutputFiles();
    }
};

TEST_P(EditconfTest, ProducesMatchingOutputStructureFile)
{
    runTest("Output file type");
}

TEST_P(EditconfTest, ProducesMatchingOutputStructureFileUsingIndexGroup)
{
    setInputFile("-n", "fragment1.ndx");
    runTest("Output file type using index group");
}

// TODO These reproduce slightly differently in double precision, and
// we don't yet have a precision-agnostic way to check on the output
// coordinates. It's better to run the tests only in single than not
// have the tests.
#if !GMX_DOUBLE
INSTANTIATE_TEST_SUITE_P(
        SinglePeptideFragments,
        EditconfTest,
        ::testing::Combine(::testing::Values("fragment1.pdb", "fragment1.gro", "fragment1.g96"),
                           ::testing::Values(efPDB, efGRO, efG96)));
#else
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(EditconfTest);
#endif

} // namespace
} // namespace test
} // namespace gmx
