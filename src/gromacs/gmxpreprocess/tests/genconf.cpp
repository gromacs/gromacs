/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Tests for genconf.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/genconf.h"

#include <filesystem>
#include <optional>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::CommandLine;
using gmx::test::ExactTextMatch;

class GenconfTest : public gmx::test::CommandLineTestBase
{
public:
    GenconfTest()
    {
        std::string confFileName =
                gmx::test::TestFileManager::getInputFilePath("spc-and-methanol.gro").string();
        commandLine().addOption("-f", confFileName);
        commandLine().addOption("-seed", "1993"); // make random operations reproducible
        setOutputFile("-o", "out.gro", ExactTextMatch());
    }

    void runTest(const CommandLine& args, const std::optional<CommandLine>& fnArgs = std::nullopt)
    {
        CommandLine& cmdline = commandLine();
        cmdline.merge(args);

        // Path arguments are not build independent, so we keep them as a separate argument
        // from regular arguments. Otherwise the rootChecker test below fails.
        if (fnArgs.has_value())
        {
            cmdline.merge(fnArgs.value());
        }

        gmx::test::TestReferenceChecker rootChecker(this->rootChecker());
        rootChecker.checkString(args.toString(), "CommandLine");

        ASSERT_EQ(0, gmx_genconf(cmdline.argc(), cmdline.argv()));

        checkOutputFiles();
    }
};

TEST_F(GenconfTest, nbox_Works)
{
    const char* const cmdline[] = { "genconf", "-nbox", "2", "1", "1" };
    runTest(CommandLine(cmdline));
}

TEST_F(GenconfTest, nbox_norenumber_Works)
{
    const char* const cmdline[] = { "genconf", "-nbox", "2", "1", "1", "-norenumber" };
    runTest(CommandLine(cmdline));
}

TEST_F(GenconfTest, nbox_dist_Works)
{
    const char* const cmdline[] = { "genconf", "-nbox", "2", "2", "3", "-dist", "0.1" };
    runTest(CommandLine(cmdline));
}

TEST_F(GenconfTest, nbox_rot_Works)
{
    const char* const cmdline[] = { "genconf", "-nbox", "2", "2", "3", "-rot" };
    runTest(CommandLine(cmdline));
}

TEST_F(GenconfTest, trj_Works)
{
    const char* const cmdline[] = { "genconf", "-nbox", "2", "2", "1" };
    CommandLine       commandLine(cmdline);

    const std::string trajFileName =
            gmx::test::TestFileManager::getInputFilePath("spc-and-methanol-traj.gro").string();
    CommandLine commandLineFnArgs;
    commandLineFnArgs.addOption("-trj", trajFileName);

    runTest(commandLine, commandLineFnArgs);
}

} // namespace
} // namespace test
} // namespace gmx
