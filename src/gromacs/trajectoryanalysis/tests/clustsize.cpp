/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests for gmx clustsize
 *
 * \todo These will be superseded by tests of the new style analysis
 * modules.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include <cstring>

#include <string>

#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/path.h"

#include "testutils/cmdlinetest.h"
#include "testutils/integrationtests.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/xvgtest.h"

namespace gmx
{

namespace test
{

namespace
{

class ClustsizeTest : public CommandLineTestBase
{
    public:
        ClustsizeTest()
        {
            setInputFile("-f", "clustsize.pdb");
        }

        void runTest(bool        mol,
                     bool        cutoff)
        {
            const char *const command[] = { "clustsize" };
            CommandLine       args      = CommandLine(command);
            double            tolerance = 1e-4;
            test::XvgMatch    xvg;

            setOutputFile("-mc", ".xvg",
                          xvg.tolerance(gmx::test::relativeToleranceAsFloatingPoint(1, tolerance)));
            setOutputFile("-nc", ".xvg",
                          xvg.tolerance(gmx::test::relativeToleranceAsFloatingPoint(1, tolerance)));
            setOutputFile("-ac", ".xvg",
                          xvg.tolerance(gmx::test::relativeToleranceAsFloatingPoint(1, tolerance)));
            setOutputFile("-hc", ".xvg",
                          xvg.tolerance(gmx::test::relativeToleranceAsFloatingPoint(1, tolerance)));

            if (mol)
            {
                setInputFile("-s", "clustsize.tpr");
                args.addOption("-mol");
            }
            else
            {
                setInputFile("-n", "clustsize.ndx");
            }
            if (cutoff)
            {
                args.addOption("-cut", "0.3");
            }
            rootChecker().checkString(args.toString(), "CommandLine");

            CommandLine &cmdline = commandLine();
            cmdline.merge(args);

            ASSERT_EQ(0, gmx_clustsize(cmdline.argc(), cmdline.argv()));

            checkOutputFiles();
        }
};

TEST_F(ClustsizeTest, NoMolDefaultCutoff)
{
    runTest(false, false);
}

TEST_F(ClustsizeTest, NoMolShortCutoff)
{
    runTest(false, true);
}

TEST_F(ClustsizeTest, MolDefaultCutoff)
{
    runTest(true, false);
}

TEST_F(ClustsizeTest, MolShortCutoff)
{
    runTest(true, true);
}

} // namespace

} // namespace

} // namespace
