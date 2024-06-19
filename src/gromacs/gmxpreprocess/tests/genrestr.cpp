/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * Tests for gmx genrestr.
 *
 * \author Kevin Boyd <kevin44boyd@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/genrestr.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

#include "testutils/cmdlinetest.h"
#include "testutils/filematchers.h"
#include "testutils/refdata.h"
#include "testutils/stdiohelper.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{
namespace test
{

class GenRestrTest : public CommandLineTestBase
{
public:
    void runTest(const std::string& interactiveCommandLineInput)
    {
        StdioTestHelper stdIoHelper(&fileManager());
        stdIoHelper.redirectStringToStdin(interactiveCommandLineInput.c_str());

        CommandLine& cmdline = commandLine();
        // Provide the name of the module to call
        std::string module[] = { "genrestr" };
        cmdline.merge(CommandLine(module));

        ASSERT_EQ(0, gmx_genrestr(cmdline.argc(), cmdline.argv()));
        checkOutputFiles();
    }
};

TEST_F(GenRestrTest, SimpleRestraintsGenerated)
{
    setInputFile("-f", "lysozyme.pdb");
    ExactTextMatch settings;
    setOutputFile("-o", "restraints.itp", TextFileMatch(settings));
    // Select c-alphas from default index options.
    std::string selection = "3";
    runTest(selection);
}
} // namespace test
} // namespace gmx
