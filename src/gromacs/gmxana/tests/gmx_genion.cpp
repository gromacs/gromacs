/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2017,2018, by the GROMACS development team, led by
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
 * Tests for adding ions
 *
 * \author Joe Jordan <e.jjordan12@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/gmxana/gmx_ana.h"

#include "gromacs/utility/futil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/conftest.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace
{

using gmx::test::CommandLine;
using gmx::test::ConfMatch;
using gmx::test::ExactTextMatch;

class GenionTest : public gmx::test::CommandLineTestBase
{
    public:
        GenionTest()
        {
            tprFilename_ = fileManager().getTemporaryFilePath("spc216.tpr");
            auto inputMdpFilename  = fileManager().getTemporaryFilePath("input.mdp");
            auto outputMdpFilename = fileManager().getTemporaryFilePath("output.mdp");
            gmx::TextWriter::writeFileFromString(inputMdpFilename, "");

            CommandLine caller = commandLine();

            caller.addOption("-f", inputMdpFilename);
            caller.addOption("-c", fileManager().getInputFilePath("spc216.gro"));
            caller.addOption("-p", fileManager().getInputFilePath("spc216.top"));
            caller.addOption("-o", tprFilename_);
            gmx_grompp(caller.argc(), caller.argv());
        }

        std::string getTprFilename() {return tprFilename_; }

        void runTest(const CommandLine &args)
        {
            CommandLine &cmdline = commandLine();
            cmdline.merge(args);

            ASSERT_EQ(0, gmx_genion(cmdline.argc(), cmdline.argv()));
            checkOutputFiles();
        }

    private:
        std::string tprFilename_;
};

TEST_F(GenionTest, add_Cation_Works)
{
    const char *const cmdline[] = {
        "genion -np 1"
    };
    commandLine().addOption("-s", getTprFilename());
    setOutputFile("-o", "out.gro", ConfMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(GenionTest, add_Anion_Works)
{
    const char *const cmdline[] = {
        "genion -nn 1"
    };
    commandLine().addOption("-s", getTprFilename());
    setOutputFile("-o", "out.gro", ConfMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(GenionTest, conc_Works)
{
    const char *const cmdline[] = {
        "genion -conc 0.1"
    };
    commandLine().addOption("-s", getTprFilename());
    setOutputFile("-o", "out.gro", ConfMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(GenionTest, neutral_Works)
{
    const char *const cmdline[] = {
        "genion -nn 1 -neutral"
    };
    commandLine().addOption("-s", getTprFilename());
    setOutputFile("-o", "out.gro", ConfMatch());
    runTest(CommandLine(cmdline));
}

TEST_F(GenionTest, update_Topology_Works)
{
    const char *const cmdline[] = {
        "genion -nn 1"
    };
    commandLine().addOption("-s", getTprFilename());
    setOutputFile("-o", "out.gro", ConfMatch());

    // TODO: Consider adding a convenience method for this.
    // Copies topology file to where it would be found as an output file, so the copied
    // .top file is used as both input and output
    std::string topFileName           = gmx::test::TestFileManager::getInputFilePath("spc216.top");
    std::string modifiableTopFileName = fileManager().getTemporaryFilePath("ions.top");
    gmx_file_copy(topFileName.c_str(), modifiableTopFileName.c_str(), true);
    setOutputFile("-p", "ions.top", ExactTextMatch());
    runTest(CommandLine(cmdline));
}

} // namespace
