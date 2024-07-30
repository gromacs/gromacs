/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
 * Tests for solvation.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/solvate.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/textreader.h"

#include "testutils/cmdlinetest.h"
#include "testutils/conftest.h"
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
using gmx::test::ConfMatch;
using gmx::test::ExactTextMatch;

class SolvateTest : public gmx::test::CommandLineTestBase
{
public:
    SolvateTest() { setOutputFile("-o", "out.gro", ConfMatch()); }

    void runTest(const CommandLine& args)
    {
        CommandLine& cmdline = commandLine();
        cmdline.merge(args);

        ASSERT_EQ(0, gmx_solvate(cmdline.argc(), cmdline.argv()));
        checkOutputFiles();
    }
};

TEST_F(SolvateTest, cs_box_Works)
{
    // use default solvent box (-cs without argument)
    const char* const cmdline[] = { "solvate", "-cs", "-box", "1.1" };
    runTest(CommandLine(cmdline));
}

TEST_F(SolvateTest, cs_cp_Works)
{
    // use default solvent box (-cs without argument)
    const char* const cmdline[] = { "solvate", "-cs" };
    setInputFile("-cp", "spc-and-methanol.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(SolvateTest, cs_cp_p_Works)
{
    // use default solvent box (-cs without argument)
    const char* const cmdline[] = { "solvate", "-cs" };
    setInputFile("-cp", "spc-and-methanol.gro");
    setModifiableInputFile("-p", "spc-and-methanol.top");

    runTest(CommandLine(cmdline));
}

TEST_F(SolvateTest, shell_Works)
{
    // use default solvent box (-cs without argument)
    const char* const cmdline[] = { "solvate", "-cs" };
    setInputFile("-cp", "spc-and-methanol.gro");
    commandLine().addOption("-shell", 1.0);

    runTest(CommandLine(cmdline));
}

TEST_F(SolvateTest, update_Topology_Works)
{
    // use solvent box with 2 solvents, check that topology has been updated
    const char* const cmdline[] = { "solvate" };
    setInputFile("-cs", "mixed_solvent.gro");
    setInputFile("-cp", "simple.gro");
    setInputAndOutputFile("-p", "simple.top", ExactTextMatch());

    runTest(CommandLine(cmdline));
}

TEST_F(SolvateTest, cs_pdb_big_box_Works)
{
    // use SPC216 solvent, but in PDB format
    const char* const cmdline[] = { "solvate", "-box", "2" };
    setInputFile("-cs", "spc216.pdb");
    runTest(CommandLine(cmdline));
}

} // namespace
} // namespace test
} // namespace gmx
