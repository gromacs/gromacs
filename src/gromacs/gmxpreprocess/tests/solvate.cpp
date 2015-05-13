/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
 * Tests for solvation.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/solvate.h"

#include "gromacs/utility/futil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace
{

using gmx::test::CommandLine;

class SolvateTest : public gmx::test::CommandLineTestBase
{
    public:
        SolvateTest()
        {
            setOutputFile("-o", "out.gro");
        }

        void runTest(const CommandLine &args)
        {
            CommandLine &cmdline = commandLine();
            cmdline.merge(args);

            ASSERT_EQ(0, gmx_solvate(cmdline.argc(), cmdline.argv()));
        }
};

TEST_F(SolvateTest, cs_box_Works)
{
    // use default solvent box (-cs without argument)
    const char *const cmdline[] = {
        "solvate", "-cs", "-box", "1.1"
    };
    runTest(CommandLine(cmdline));
}

TEST_F(SolvateTest, cs_cp_Works)
{
    // use default solvent box (-cs without argument)
    const char *const cmdline[] = {
        "solvate", "-cs"
    };
    setInputFile("-cp", "spc-and-methanol.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(SolvateTest, cs_cp_p_Works)
{
    // use default solvent box (-cs without argument)
    const char *const cmdline[] = {
        "solvate", "-cs"
    };
    setInputFile("-cp", "spc-and-methanol.gro");

    // TODO: Consider adding a convenience method for this.
    std::string topFileName           = fileManager().getInputFilePath("spc-and-methanol.top");
    std::string modifiableTopFileName = fileManager().getTemporaryFilePath(".top");
    gmx_file_copy(topFileName.c_str(), modifiableTopFileName.c_str(), true);
    commandLine().addOption("-p", modifiableTopFileName);

    runTest(CommandLine(cmdline));
}

} // namespace
