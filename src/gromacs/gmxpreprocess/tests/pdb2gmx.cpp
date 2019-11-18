/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * Tests for pdb2gmx
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/pdb2gmx.h"

#include "gromacs/fileio/filetypes.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/textreader.h"

#include "testutils/cmdlinetest.h"
#include "testutils/conftest.h"
#include "testutils/filematchers.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

using test::CommandLine;

//! Test parameter struct.
using CommandLineOptionParams =
        std::tuple<std::string, std::string, std::string, std::string, std::string, std::string, int>;

/*! \brief Strings containing regular expressions for lines to skip
 * when matching.
 *
 * \todo It would be preferable to just scrub the content that actually
 * varies, but we don't use enough regular expression support for that
 * yet.
 *
 * Note that the "\n" are needed so these regular expressions match
 * Windows line endings. */
std::vector<std::string> c_regexStringsToSkip = { "^;[[:blank:]] *File '.*' was generated.*\n",
                                                  "^;[[:blank:]]*By user:.*\n",
                                                  "^;[[:blank:]]*On host:.*\n",
                                                  "^;[[:blank:]]*At date:.*\n",
                                                  "^;[[:blank:]]*:-\\).*\\(-:.*\n",
                                                  "^;[[:blank:]]*Executable:.*\n",
                                                  "^;[[:blank:]]*Data prefix:.*\n",
                                                  "^;[[:blank:]]*Working dir:.*\n",
                                                  "^;[[:blank:]]*pdb2gmx.*-test.*\n" };
//! Compiled regular expressions for lines to skip when matching.
FilteringExactTextMatch c_textMatcher(c_regexStringsToSkip);

class Pdb2gmxTest : public test::CommandLineTestBase, public ::testing::WithParamInterface<CommandLineOptionParams>
{
public:
    Pdb2gmxTest()
    {
        int outputFileType = std::get<6>(GetParam());
        if (outputFileType == efPDB)
        {
            // If we're writing PDB output, we are interested in
            // testing things like TER records and chain IDs.
            std::string outputfile = "conf.";
            outputfile += ftp2ext(outputFileType);
            ExactTextMatch settings;
            setOutputFile("-o", outputfile.c_str(), TextFileMatch(settings));
        }
        else
        {
            setOutputFile("-o", "conf.gro", ConfMatch());
        }
        setOutputFile("-p", "topol.top", TextFileMatch(c_textMatcher));
    }

    void runTest(const CommandLine& args)
    {
        CommandLine& cmdline = commandLine();
        cmdline.merge(args);

        TestReferenceChecker rootChecker(this->rootChecker());

        ASSERT_EQ(0, CommandLineTestHelper::runModuleFactory(&pdb2gmxInfo::create, &cmdline));

        checkOutputFiles();
    }
};

TEST_P(Pdb2gmxTest, ProducesMatchingTopology)
{
    const auto& params    = GetParam();
    std::string cmdline[] = { "pdb2gmx",   "-ignh",
                              "-ff",       std::get<0>(params),
                              "-water",    std::get<1>(params),
                              "-vsite",    std::get<2>(params),
                              "-chainsep", std::get<3>(params),
                              "-merge",    std::get<4>(params) };
    setInputFile("-f", std::get<5>(params));
    runTest(CommandLine(cmdline));
}

// These tests are still rather slow when run with TSAN, so in the
// CMakeLists.txt file we split them into separtae test binaries.

#if OPLSAA
INSTANTIATE_TEST_CASE_P(ForOplsaa,
                        Pdb2gmxTest,
                        ::testing::Combine(::testing::Values("oplsaa"),
                                           ::testing::Values("tip3p", "tip4p", "tip5p"),
                                           ::testing::Values("none", "h"),
                                           ::testing::Values("id_or_ter"),
                                           ::testing::Values("no"),
                                           ::testing::Values("fragment1.pdb",
                                                             "fragment2.pdb",
                                                             "fragment3.pdb",
                                                             "fragment4.pdb"),
                                           ::testing::Values(efGRO)));
#endif

#if GROMOS
INSTANTIATE_TEST_CASE_P(ForGromos43a1,
                        Pdb2gmxTest,
                        ::testing::Combine(::testing::Values("gromos43a1"),
                                           ::testing::Values("spc", "spce"),
                                           ::testing::Values("none", "h"),
                                           ::testing::Values("id_or_ter"),
                                           ::testing::Values("no"),
                                           ::testing::Values("fragment1.pdb",
                                                             "fragment2.pdb",
                                                             "fragment3.pdb",
                                                             "fragment4.pdb"),
                                           ::testing::Values(efGRO)));

INSTANTIATE_TEST_CASE_P(ForGromos53a6,
                        Pdb2gmxTest,
                        ::testing::Combine(::testing::Values("gromos53a6"),
                                           ::testing::Values("spc", "spce"),
                                           ::testing::Values("none", "h"),
                                           ::testing::Values("id_or_ter"),
                                           ::testing::Values("no"),
                                           ::testing::Values("fragment1.pdb",
                                                             "fragment2.pdb",
                                                             "fragment3.pdb",
                                                             "fragment4.pdb"),
                                           ::testing::Values(efGRO)));
#endif

#if AMBER
INSTANTIATE_TEST_CASE_P(ForAmber99sb_ildn,
                        Pdb2gmxTest,
                        ::testing::Combine(::testing::Values("amber99sb-ildn"),
                                           ::testing::Values("tip3p"),
                                           ::testing::Values("none", "h"),
                                           ::testing::Values("id_or_ter"),
                                           ::testing::Values("no"),
                                           ::testing::Values("fragment1.pdb",
                                                             "fragment2.pdb",
                                                             "fragment3.pdb",
                                                             "fragment4.pdb"),
                                           ::testing::Values(efGRO)));
#endif

#if CHARMM
INSTANTIATE_TEST_CASE_P(ForCharmm27,
                        Pdb2gmxTest,
                        ::testing::Combine(::testing::Values("charmm27"),
                                           ::testing::Values("tip3p"),
                                           ::testing::Values("none", "h"),
                                           ::testing::Values("id_or_ter"),
                                           ::testing::Values("no"),
                                           ::testing::Values("fragment1.pdb",
                                                             "fragment2.pdb",
                                                             "fragment3.pdb",
                                                             "fragment4.pdb"),
                                           ::testing::Values(efGRO)));


INSTANTIATE_TEST_CASE_P(ChainSep,
                        Pdb2gmxTest,
                        ::testing::Combine(::testing::Values("charmm27"),
                                           ::testing::Values("tip3p"),
                                           ::testing::Values("none"),
                                           ::testing::Values("id", "ter", "id_or_ter", "id_and_ter"),
                                           ::testing::Values("all", "no"),
                                           ::testing::Values("chainTer.pdb"),
                                           ::testing::Values(efGRO)));

INSTANTIATE_TEST_CASE_P(ChainChanges,
                        Pdb2gmxTest,
                        ::testing::Combine(::testing::Values("charmm27"),
                                           ::testing::Values("tip3p"),
                                           ::testing::Values("none"),
                                           ::testing::Values("id", "ter", "id_or_ter", "id_and_ter"),
                                           ::testing::Values("no"),
                                           ::testing::Values("two-fragments.pdb"),
                                           ::testing::Values(efPDB)));
#endif

} // namespace
} // namespace test
} // namespace gmx
