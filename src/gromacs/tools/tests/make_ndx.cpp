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
 * Tests for gmx make_ndx.
 *
 * \author R. Thomas Ullmann <tullman@gwdg.de>
 */
#include "gmxpre.h"

#include "gromacs/tools/make_ndx.h"

#include <string>

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/stdiohelper.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

class GmxMakeNdx : public gmx::test::CommandLineTestBase
{
public:
    /*! \brief runs the test for the filename prefix \p sysName from the simulation data base
        \param[in]   sysName        file name prefix sysName.g96 must be present in the simulation database
        \param[in]   stdinContent   commands passed to the interactive prompt
        \param[in]   useIndexInputOnly   whether to provide just an index as input (no structure file)
     */
    void runTest(const std::string& sysName,
                 const std::string& stdinContent      = "q\n",
                 bool               useIndexInputOnly = false)
    {
        auto&       cmdline        = commandLine();
        auto        groFileName    = sysName + ".g96";
        std::string ndxFileName    = sysName + ".ndx";
        std::string outputFileName = sysName + "_output.ndx";

        if (useIndexInputOnly)
        {
            setInputFile("-n", ndxFileName);
        }
        else
        {
            setInputFile("-f", groFileName);
        }
        setOutputFile("-o", outputFileName.c_str(), gmx::test::ExactTextMatch());

        gmx::test::StdioTestHelper stdioHelper(&fileManager());
        stdioHelper.redirectStringToStdin(stdinContent.c_str());

        ASSERT_EQ(0, gmx_make_ndx(cmdline.argc(), cmdline.argv()));
        checkOutputFiles();
    }
};

TEST_F(GmxMakeNdx, WritesDefaultProteinIndexGroups)
{
    std::string sysName("villin");
    runTest(sysName);
}

TEST_F(GmxMakeNdx, HandlesNoStructureInput)
{
    std::string sysName("alanine_vacuo");
    runTest(sysName, "1|2\nq\n", true);
}

TEST_F(GmxMakeNdx, HandlesNotProtein)
{
    std::string sysName("spc-dimer");
    runTest(sysName, "q\n", true);
}

TEST_F(GmxMakeNdx, HandlesEmptyIndexResult)
{
    std::string sysName("alanine_vacuo");
    runTest(sysName, "4&8\nq\n", true);
}

TEST_F(GmxMakeNdx, HandlesEmptyIndexFile)
{
    std::string sysName("spc-dimer");
    runTest(sysName, "del 0\nq\n", true);
}

TEST_F(GmxMakeNdx, Splitres)
{
    std::string sysName("spc-dimer");
    runTest(sysName, "splitres 1\nq\n", false);
}

TEST_F(GmxMakeNdx, Splitat)
{
    std::string sysName("spc-dimer");
    runTest(sysName, "splitat 1\nq\n", false);
}

} // namespace
} // namespace test
} // namespace gmx
