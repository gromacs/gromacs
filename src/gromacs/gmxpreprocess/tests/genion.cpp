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
 * Tests for genion.
 *
 * \author Vytautas Gapsys <vgapsys@gwdg.de>
 * \author Christian Blau <blau@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/genion.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/conftest.h"
#include "testutils/refdata.h"
#include "testutils/stdiohelper.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"


namespace gmx
{
namespace test
{
namespace
{

class GenionTest : public CommandLineTestBase
{
public:
    GenionTest()
    {
        CommandLine caller = commandLine();

        const std::string mdpInputFileName(fileManager().getTemporaryFilePath("input.mdp").string());
        TextWriter::writeFileFromString(
                mdpInputFileName,
                "verlet-buffer-tolerance =-1\nrcoulomb=0.5\nrvdw = 0.5\nrlist = 0.5\n");
        caller.addOption("-f", mdpInputFileName);
        caller.addOption("-c", TestFileManager::getInputFilePath("spc216_with_methane.gro").string());
        caller.addOption("-p", TestFileManager::getInputFilePath("spc216_with_methane.top").string());
        caller.addOption("-o", tprFileName_);

        gmx_grompp(caller.argc(), caller.argv());

        setOutputFile("-o", "out.gro", ConfMatch());
    }

    void runTest(const CommandLine& args, const std::string& interactiveCommandLineInput)
    {
        StdioTestHelper stdIoHelper(&fileManager());
        stdIoHelper.redirectStringToStdin(interactiveCommandLineInput.c_str());

        CommandLine& cmdline = commandLine();
        cmdline.addOption("-s", tprFileName_);
        cmdline.merge(args);

        ASSERT_EQ(0, gmx_genion(cmdline.argc(), cmdline.argv()));
        checkOutputFiles();
    }

private:
    const std::string tprFileName_ =
            fileManager().getTemporaryFilePath("spc216_with_methane.tpr").string();
};

TEST_F(GenionTest, HighConcentrationIonPlacement)
{
    const char* const cmdline[] = { "genion", "-seed", "1997", "-conc", "1.0", "-rmin", "0.6" };

    runTest(CommandLine(cmdline), "Water");
}

TEST_F(GenionTest, NoIonPlacement)
{
    const char* const cmdline[] = { "genion", "-seed", "1997", "-conc", "0.0", "-rmin", "0.6" };

    runTest(CommandLine(cmdline), "Water");
}

} // namespace
} // namespace test
} // namespace gmx
