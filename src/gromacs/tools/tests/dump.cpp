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

protected:
    // TODO this is changed in newer googletest versions
    //! Prepare shared resources.
    static void SetUpTestSuite() { s_tprFileHandle = new TprAndFileManager("lysozyme"); }
    //! Clean up shared resources.
    static void TearDownTestSuite()
    {
        delete s_tprFileHandle;
        s_tprFileHandle = nullptr;
    }
    //! Storage for opened file handles.
    static TprAndFileManager* s_tprFileHandle;
};

TprAndFileManager* DumpTest::s_tprFileHandle = nullptr;

void DumpTest::runTest(CommandLine* cmdline)
{
    EXPECT_EQ(0, gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::DumpInfo::create, cmdline));
}

TEST_F(DumpTest, WorksWithTpr)
{
    const char* const command[] = { "dump", "-s", s_tprFileHandle->tprName().c_str() };
    CommandLine       cmdline(command);
    runTest(&cmdline);
}

TEST_F(DumpTest, WorksWithTprAndMdpWriting)
{
    TestFileManager fileManager;
    std::string     mdpName = fileManager.getTemporaryFilePath("output.mdp").string();
    const char* const command[] = { "dump", "-s", s_tprFileHandle->tprName().c_str(), "-om", mdpName.c_str() };
    CommandLine cmdline(command);
    runTest(&cmdline);
}

} // namespace test

} // namespace gmx
