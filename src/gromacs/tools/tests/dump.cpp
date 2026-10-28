/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Tests for functionality of the "dump" tool.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/tools/dump.h"

#include <filesystem>
#include <functional>
#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"
#include "testutils/tprfilegenerator.h"

namespace gmx
{

namespace test
{

class DumpTest : public ::testing::Test
{
public:
    //! Run test case.
    static void runTest(CommandLine* cmdline);

    //! Prepare tpr resource shared over all test cases (for efficiency)
    static void SetUpTestSuite()
    {
        // This requires earlier setup work to decide what the name of
        // the temporary directory is.
        s_tprFileManager = std::make_unique<TprAndFileManager>("spc-and-methanol");
    }
    //! Storage for tpr file manager
    static std::unique_ptr<TprAndFileManager> s_tprFileManager;
};

std::unique_ptr<TprAndFileManager> DumpTest::s_tprFileManager;

void DumpTest::runTest(CommandLine* cmdline)
{
    EXPECT_EQ(0, gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::DumpInfo::create, cmdline));
}

TEST_F(DumpTest, WorksWithTpr)
{
    const std::string tprName   = s_tprFileManager->tprName().string();
    const char* const command[] = { "dump", "-s", tprName.c_str() };
    CommandLine       cmdline(command);
    runTest(&cmdline);
}

TEST_F(DumpTest, WorksWithTprAndMdpWriting)
{
    TestFileManager   fileManager;
    const std::string tprName   = s_tprFileManager->tprName().string();
    const std::string mdpName   = fileManager.getTemporaryFilePath("output.mdp").string();
    const char* const command[] = { "dump", "-s", tprName.c_str(), "-om", mdpName.c_str() };
    CommandLine       cmdline(command);
    runTest(&cmdline);
}

} // namespace test

} // namespace gmx
