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
 * Tests for insertion of molecules.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/gmxpreprocess/insert-molecules.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"

namespace
{

using gmx::test::CommandLine;

class InsertMoleculesTest : public gmx::test::CommandLineTestBase
{
    public:
        InsertMoleculesTest()
        {
            setOutputFile("-o", "out.gro");
        }

        void runTest(const CommandLine &args)
        {
            CommandLine &cmdline = commandLine();
            cmdline.merge(args);

            gmx::test::TestReferenceChecker rootChecker(this->rootChecker());
            rootChecker.checkString(args.toString(), "CommandLine");

            ASSERT_EQ(0, gmx::test::CommandLineTestHelper::runModule(
                              &gmx::InsertMoleculesInfo::create, &cmdline));

            checkOutputFiles();
        }
};

TEST_F(InsertMoleculesTest, InsertsMoleculesIntoExistingConfiguration)
{
    const char *const cmdline[] = {
        "insert-molecules", "-nmol", "1"
    };
    setInputFile("-f", "spc-and-methanol.gro");
    setInputFile("-ci", "x2.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(InsertMoleculesTest, InsertsMoleculesIntoEmptyBox)
{
    const char *const cmdline[] = {
        "insert-molecules", "-box", "4", "-nmol", "5"
    };
    setInputFile("-ci", "x2.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(InsertMoleculesTest, InsertsMoleculesIntoEnlargedBox)
{
    const char *const cmdline[] = {
        "insert-molecules", "-box", "4", "-nmol", "2"
    };
    setInputFile("-f", "spc-and-methanol.gro");
    setInputFile("-ci", "x.gro");
    runTest(CommandLine(cmdline));
}

TEST_F(InsertMoleculesTest, InsertsMoleculesIntoFixedPositions)
{
    const char *const cmdline[] = {
        "insert-molecules", "-box", "4"
    };
    const char *const positions[] = {
        "0.0  0.0  0.0",
        "1.0  2.0  3.0",
        "0.99 2.01 3.0",
        "2.0  1.0  2.0"
    };
    setInputFile("-ci", "x0.gro");
    setInputFileContents("-ip", "dat", positions);
    runTest(CommandLine(cmdline));
}

} // namespace
