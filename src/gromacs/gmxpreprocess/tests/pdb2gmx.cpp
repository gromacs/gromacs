/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Integration tests for pdb2gmx.
 *
 * For now, only covers non-interactive -merge and -chainsep
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/gmxpreprocess/pdb2gmx.h"

#include "testutils/cmdlinetest.h"

namespace
{

using gmx::test::CommandLine;

class Pdb2gmxTest : public gmx::test::CommandLineTestBase
{
    public:
        Pdb2gmxTest()
        {
            setInputFile("-f", "lysozyme-chainid-and-ter.pdb");
            setOutputFile("-o", "conf.gro");
            setOutputFile("-p", "topol.top");
            /* We ignore a bunch of extraneous output files */
        }

        void runTest(const CommandLine &args)
        {
            CommandLine &cmdline = commandLine();
            cmdline.merge(args);

            ASSERT_EQ(0, gmx_pdb2gmx(cmdline.argc(), cmdline.argv()));

            /* Checks that they were written and that their content
               matches the reference copies */
            checkOutputFiles();
        }
};

TEST_F(Pdb2gmxTest, MergeNo)
{
    const char *const cmdline[] = {
        "pdb2gmx", "-testmode", "-ff", "oplsaa", "-water", "tip3p", "-merge", "no"
    };
    runTest(CommandLine(cmdline));
}

TEST_F(Pdb2gmxTest, MergeAll)
{
    const char *const cmdline[] = {
        "pdb2gmx", "-testmode", "-ff", "oplsaa", "-water", "tip3p", "-merge", "all"
    };
    runTest(CommandLine(cmdline));
}

TEST_F(Pdb2gmxTest, ChainsepIdOrTer)
{
    const char *const cmdline[] = {
        "pdb2gmx", "-testmode", "-ff", "oplsaa", "-water", "tip3p", "-chainsep", "id_or_ter"
    };
    runTest(CommandLine(cmdline));
}

TEST_F(Pdb2gmxTest, ChainsepIdAndTer)
{
    const char *const cmdline[] = {
        "pdb2gmx", "-testmode", "-ff", "oplsaa", "-water", "tip3p", "-chainsep", "id_and_ter"
    };
    runTest(CommandLine(cmdline));
}

TEST_F(Pdb2gmxTest, ChainsepId)
{
    const char *const cmdline[] = {
        "pdb2gmx", "-testmode", "-ff", "oplsaa", "-water", "tip3p", "-chainsep", "id"
    };
    runTest(CommandLine(cmdline));
}

TEST_F(Pdb2gmxTest, ChainsepTer)
{
    const char *const cmdline[] = {
        "pdb2gmx", "-testmode", "-ff", "oplsaa", "-water", "tip3p", "-chainsep", "ter"
    };
    runTest(CommandLine(cmdline));
}

} // namespace
