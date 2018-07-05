/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
using test::ConfMatch;

//! Test parameter struct.
using CommandLineOptionParams = std::tuple<std::string, std::string, std::string, std::string,
                                           std::string, std::string>;

class Pdb2gmxTest : public test::CommandLineTestBase,
                    public ::testing::WithParamInterface<CommandLineOptionParams>
{
    public:
        Pdb2gmxTest()
        {
            // TODO It would be preferable to just scrub the content
            // that actually varies, but we don't have enough regular
            // expression support for that yet.

            // Note that the "\n" are needed so these regular
            // expressions match Windows line endings.
            textMatcher_.addRegexToSkip("^;[[:blank:]] *File '.*' was generated.*\n");
            textMatcher_.addRegexToSkip("^;[[:blank:]]*By user:.*\n");
            textMatcher_.addRegexToSkip("^;[[:blank:]]*On host:.*\n");
            textMatcher_.addRegexToSkip("^;[[:blank:]]*At date:.*\n");
            textMatcher_.addRegexToSkip("^;[[:blank:]]*:-\\).*\\(-:.*\n");
            textMatcher_.addRegexToSkip("^;[[:blank:]]*Executable:.*\n");
            textMatcher_.addRegexToSkip("^;[[:blank:]]*Data prefix:.*\n");
            textMatcher_.addRegexToSkip("^;[[:blank:]]*Working dir:.*\n");
            textMatcher_.addRegexToSkip("^;[[:blank:]]*pdb2gmx-test.*\n");
            setOutputFile("-o", "conf.gro", ConfMatch());
            setOutputFile("-p", "topol.top", TextFileMatch(textMatcher_));
        }

        void runTest(const CommandLine &args)
        {
            CommandLine &cmdline = commandLine();
            cmdline.merge(args);

            TestReferenceChecker rootChecker(this->rootChecker());

            ASSERT_EQ(0, CommandLineTestHelper::runModuleFactory(&pdb2gmxInfo::create, &cmdline));

            checkOutputFiles();
        }
        FilteringExactTextMatch textMatcher_;
};

TEST_P(Pdb2gmxTest, ProducesMatchingTopology)
{
    const auto &params    = GetParam();
    std::string cmdline[] = {
        "pdb2gmx", "-ignh", "-ff", std::get<0>(params), "-water", std::get<1>(params), "-vsite", std::get<2>(params),
        "-chainsep", std::get<3>(params), "-merge", std::get<4>(params)
    };
    setInputFile("-f", std::get<5>(params));
    runTest(CommandLine(cmdline));
}

INSTANTIATE_TEST_CASE_P(ForOplsaa, Pdb2gmxTest,
                            ::testing::Combine
                            (::testing::Values("oplsaa"),
                                ::testing::Values("tip3p", "tip4p", "tip5p"),
                                ::testing::Values("none", "h"),
                                ::testing::Values("id_or_ter"),
                                ::testing::Values("no"),
                                ::testing::Values("fragment1.pdb", "fragment2.pdb", "fragment3.pdb", "fragment4.pdb"))
                        );

INSTANTIATE_TEST_CASE_P(ForGromos43a1, Pdb2gmxTest,
                            ::testing::Combine
                            (::testing::Values("gromos43a1"),
                                ::testing::Values("spc", "spce"),
                                ::testing::Values("none", "h"),
                                ::testing::Values("id_or_ter"),
                                ::testing::Values("no"),
                                ::testing::Values("fragment1.pdb", "fragment2.pdb", "fragment3.pdb", "fragment4.pdb"))
                        );

INSTANTIATE_TEST_CASE_P(ForGromos53a6, Pdb2gmxTest,
                            ::testing::Combine
                            (::testing::Values("gromos53a6"),
                                ::testing::Values("spc", "spce"),
                                ::testing::Values("none", "h"),
                                ::testing::Values("id_or_ter"),
                                ::testing::Values("no"),
                                ::testing::Values("fragment1.pdb", "fragment2.pdb", "fragment3.pdb", "fragment4.pdb"))
                        );

INSTANTIATE_TEST_CASE_P(ForAmber99sb_ildn, Pdb2gmxTest,
                            ::testing::Combine
                            (::testing::Values("amber99sb-ildn"),
                                ::testing::Values("tip3p"),
                                ::testing::Values("none", "h"),
                                ::testing::Values("id_or_ter"),
                                ::testing::Values("no"),
                                ::testing::Values("fragment1.pdb", "fragment2.pdb", "fragment3.pdb", "fragment4.pdb"))
                        );

INSTANTIATE_TEST_CASE_P(ForCharmm27, Pdb2gmxTest,
                            ::testing::Combine
                            (::testing::Values("charmm27"),
                                ::testing::Values("tip3p"),
                                ::testing::Values("none", "h"),
                                ::testing::Values("id_or_ter"),
                                ::testing::Values("no"),
                                ::testing::Values("fragment1.pdb", "fragment2.pdb", "fragment3.pdb", "fragment4.pdb"))
                        );


INSTANTIATE_TEST_CASE_P(ChainSep, Pdb2gmxTest,
                            ::testing::Combine
                            (::testing::Values("charmm27"),
                                ::testing::Values("tip3p"),
                                ::testing::Values("none"),
                                ::testing::Values("id", "ter", "id_or_ter", "id_and_ter"),
                                ::testing::Values("all", "no"),
                                ::testing::Values("chainTer.pdb"))
                        );

} // namespace
} // namespace test
} // namespace gmx
