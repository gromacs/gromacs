/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * Tests for pdb2gmx
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/pdb2gmx.h"

#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "testutils/cmdlinetest.h"
#include "testutils/conftest.h"
#include "testutils/filematchers.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
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
        std::tuple<std::string, std::string, std::string, std::string, std::string, std::string, int, bool>;

/*! \brief Strings containing regular expressions for lines to skip
 * when matching.
 *
 * \todo It would be preferable to just scrub the content that actually
 * varies, but we don't use enough regular expression support for that
 * yet. */
std::vector<std::string> c_regexStringsToSkip = { "^;[[:blank:]] *File '.*' was generated.*",
                                                  "^;[[:blank:]]*By user:.*",
                                                  "^;[[:blank:]]*On host:.*",
                                                  "^;[[:blank:]]*At date:.*",
                                                  "^;[[:blank:]]*:-\\).*\\(-:.*",
                                                  "^;[[:blank:]]*Executable:.*",
                                                  "^;[[:blank:]]*Data prefix:.*",
                                                  "^;[[:blank:]]*Working dir:.*",
                                                  "^;[[:blank:]]*pdb2gmx.*-test.*" };
//! Compiled regular expressions for lines to skip when matching.
FilteringExactTextMatch c_textMatcher(c_regexStringsToSkip, false, true);

class Pdb2gmxTest : public test::CommandLineTestBase, public ::testing::WithParamInterface<CommandLineOptionParams>
{
public:
    Pdb2gmxTest()
    {
        int outputFileType = std::get<6>(GetParam());
        // If the file type of the output configuration is one
        // commonly used (ie. pdb, gro), then check its content,
        // otherwise just check the output file exists.
        if (outputFileType == efPDB)
        {
            // If we're writing PDB output, we are interested in
            // testing things like TER records and chain IDs.
            std::string outputfile = "conf.";
            outputfile += ftp2ext(outputFileType);
            ExactTextMatch settings;
            setOutputFile("-o", outputfile.c_str(), TextFileMatch(settings));
        }
        else if (outputFileType == efGRO)
        {
            setOutputFile("-o", "conf.gro", ConfMatch().matchFullConfiguration(std::get<7>(GetParam())));
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

TEST_P(Pdb2gmxTest, Runs)
{
    checkTestNameLength();
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

//! Help GoogleTest name our test cases
std::string namesOfTests(const testing::TestParamInfo<Pdb2gmxTest::ParamType>& info)
{
    const auto& param = info.param;

    std::string testName = formatString(
            "ff_%s_"
            "%s_"
            "vsite_%s_"
            "%s_"
            "merge_%s_"
            "%s_"
            "format_%s_"
            "match_%s",
            std::get<0>(param).c_str(),
            std::get<1>(param).c_str(),
            std::get<2>(param).c_str(),
            std::get<3>(param).c_str(),
            std::get<4>(param).c_str(),
            std::get<5>(param).c_str(),
            ftp2ext(std::get<6>(param)),
            std::get<7>(param) ? "full" : "file");

    // Note that the returned names must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = replaceAll(testName, "-", "");
    testName = replaceAll(testName, ".", "");
    return testName;
}

// These tests are still rather slow when run with TSAN, so in the
// CMakeLists.txt file we split them into separtae test binaries.

#if OPLSAA
INSTANTIATE_TEST_SUITE_P(
        Oplsaa,
        Pdb2gmxTest,
        ::testing::Combine(::testing::Values("oplsaa"),
                           ::testing::Values("tip3p", "tip4p", "tip5p"),
                           ::testing::Values("none", "h"),
                           ::testing::Values("id_or_ter"),
                           ::testing::Values("no"),
                           ::testing::Values("A.pdb", "B.pdb", "C.pdb", "D.pdb", "E.pdb"),
                           ::testing::Values(efGRO),
                           ::testing::Values(false)),
        namesOfTests);
#endif

#if GROMOS
INSTANTIATE_TEST_SUITE_P(
        G43a1,
        Pdb2gmxTest,
        ::testing::Combine(::testing::Values("gromos43a1"),
                           ::testing::Values("spc", "spce"),
                           ::testing::Values("none", "h"),
                           ::testing::Values("id_or_ter"),
                           ::testing::Values("no"),
                           ::testing::Values("A.pdb", "B.pdb", "C.pdb", "D.pdb", "E.pdb"),
                           ::testing::Values(efGRO),
                           ::testing::Values(false)),
        namesOfTests);

INSTANTIATE_TEST_SUITE_P(
        G53a6,
        Pdb2gmxTest,
        ::testing::Combine(::testing::Values("gromos53a6"),
                           ::testing::Values("spc", "spce"),
                           ::testing::Values("none", "h"),
                           ::testing::Values("id_or_ter"),
                           ::testing::Values("no"),
                           ::testing::Values("A.pdb", "B.pdb", "C.pdb", "D.pdb", "E.pdb"),
                           ::testing::Values(efGRO),
                           ::testing::Values(false)),
        namesOfTests);
#endif

#if AMBER
INSTANTIATE_TEST_SUITE_P(
        Amber,
        Pdb2gmxTest,
        ::testing::Combine(::testing::Values("amber99sb-ildn"),
                           ::testing::Values("tip3p"),
                           ::testing::Values("none", "h"),
                           ::testing::Values("id_or_ter"),
                           ::testing::Values("no"),
                           ::testing::Values("A.pdb", "B.pdb", "C.pdb", "D.pdb", "E.pdb"),
                           ::testing::Values(efGRO),
                           ::testing::Values(false)),
        namesOfTests);
INSTANTIATE_TEST_SUITE_P(AmberTip4p,
                         Pdb2gmxTest,
                         ::testing::Combine(::testing::Values("amber99sb-ildn"),
                                            ::testing::Values("tip4p"),
                                            ::testing::Values("none"),
                                            ::testing::Values("id_or_ter"),
                                            ::testing::Values("no"),
                                            ::testing::Values("tip4p.pdb"),
                                            ::testing::Values(efGRO),
                                            ::testing::Values(true)),
                         namesOfTests);
#endif

#if CHARMM
INSTANTIATE_TEST_SUITE_P(
        Charmm,
        Pdb2gmxTest,
        ::testing::Combine(
                ::testing::Values("charmm27"),
                ::testing::Values("tip3p"),
                ::testing::Values("none", "h"),
                ::testing::Values("id_or_ter"),
                ::testing::Values("no"),
                ::testing::Values("A.pdb", "B.pdb", "C.pdb", "D.pdb", "E.pdb", "monomer.pdb"),
                ::testing::Values(efGRO),
                ::testing::Values(false)),
        namesOfTests);


INSTANTIATE_TEST_SUITE_P(
        ChainSep,
        Pdb2gmxTest,
        ::testing::Combine(::testing::Values("charmm27"),
                           ::testing::Values("tip3p"),
                           ::testing::Values("none"),
                           ::testing::Values("id", "ter", "id_or_ter", "id_and_ter"),
                           ::testing::Values("all", "no"),
                           ::testing::Values("chainTer.pdb"),
                           ::testing::Values(efGRO),
                           ::testing::Values(false)),
        namesOfTests);

INSTANTIATE_TEST_SUITE_P(
        ChainChanges,
        Pdb2gmxTest,
        ::testing::Combine(::testing::Values("charmm27"),
                           ::testing::Values("tip3p"),
                           ::testing::Values("none"),
                           ::testing::Values("id", "ter", "id_or_ter", "id_and_ter"),
                           ::testing::Values("no"),
                           ::testing::Values("fragments.pdb"),
                           ::testing::Values(efPDB),
                           ::testing::Values(false)),
        namesOfTests);

INSTANTIATE_TEST_SUITE_P(Cyclic,
                         Pdb2gmxTest,
                         ::testing::Combine(::testing::Values("charmm27"),
                                            ::testing::Values("tip3p"),
                                            ::testing::Values("none"),
                                            ::testing::Values("id_or_ter"),
                                            ::testing::Values("no", "all"),
                                            ::testing::Values("cyc-rna.pdb", "cyc-prot.pdb"),
                                            ::testing::Values(efGRO),
                                            ::testing::Values(false)),
                         namesOfTests);
#endif

} // namespace
} // namespace test
} // namespace gmx
