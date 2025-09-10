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
 * Tests for gmx mindist.
 *
 * \author Kevin Boyd <kevin.boyd@uconn.edu>
 */

#include "gmxpre.h"

#include <cstdio>
#include <cstdlib>

#include <string>

#include <gtest/gtest.h>

#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/textreader.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/stdiohelper.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"
#include "testutils/xvgtest.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::CommandLine;
using gmx::test::StdioTestHelper;
using gmx::test::XvgMatch;

class MindistTest : public gmx::test::CommandLineTestBase
{
public:
    MindistTest()
    {
        setInputFile("-f", "mindist_coords.gro");
        setInputFile("-s", "mindist_coords.gro");
        setInputFile("-n", "mindist.ndx");
    }

    void runTest(const CommandLine& args, const char* stringForStdin)
    {
        StdioTestHelper stdioHelper(&fileManager());
        stdioHelper.redirectStringToStdin(stringForStdin);

        CommandLine& cmdline = commandLine();
        cmdline.merge(args);
        ASSERT_EQ(0, gmx_mindist(cmdline.argc(), cmdline.argv()));
        checkOutputFiles();
    }
};

/* mindist_coords.pdb has 3 beads spaced out in a 5 nm box, with the same yz coordinates
 * and x coordinates of 1, 4, and 4.5. Indices are as follows
 * index 0 : atom 1
 * index 1 : atom 2
 * index 2 : atom 3
 * index 3 : atoms (1 ,2)
 * index 4 : atoms (2, 3)
 * index 5 : atoms (1, 2, 3)
 */

// Mindist between beads 0 and 1 should = 2 (across periodic boundaries)
TEST_F(MindistTest, mindistWorksWithSingleAtoms)
{
    setOutputFile("-od", "mindist.xvg", XvgMatch());
    const char* const cmdline[] = { "mindist" };
    const char* const stdIn     = "0 1";
    runTest(CommandLine(cmdline), stdIn);
}

// Mindist between group (0, 1) and bead 2 should = 0.5
TEST_F(MindistTest, mindistWorksWithMultipleAtoms)
{
    setOutputFile("-od", "mindist.xvg", XvgMatch());
    const char* const cmdline[] = { "mindist" };
    const char* const stdIn     = "2 3";
    runTest(CommandLine(cmdline), stdIn);
}

/* Should have 0 contacts with default cutoff */
TEST_F(MindistTest, mindistDoesNotPickUpContacts)
{
    setOutputFile("-on", "ncontacts.xvg", XvgMatch());
    const char* const cmdline[] = { "mindist" };
    const char* const stdIn     = "0 1";
    runTest(CommandLine(cmdline), stdIn);
}

/* Should pick up one contact with 2.5 nm cutoff */
TEST_F(MindistTest, mindistPicksUpContacts)
{
    setOutputFile("-on", "ncontacts.xvg", XvgMatch());
    const char* const cmdline[] = {
        "mindist",
        "-d",
        "2.5",
    };
    const char* const stdIn = "0 1";
    runTest(CommandLine(cmdline), stdIn);
}

TEST_F(MindistTest, ngWorks)
{
    setOutputFile("-od", "mindist.xvg", XvgMatch());
    const char* const cmdline[] = {
        "mindist",
        "-ng",
        "2",
    };
    const char* const stdIn = "0 1 2";
    runTest(CommandLine(cmdline), stdIn);
}

// 2 contacts within this cutoff, but only one should be reported
TEST_F(MindistTest, groupWorks)
{
    setOutputFile("-on", "ncontacts.xvg", XvgMatch());
    const char* const cmdline[] = { "mindist", "-group", "-d", "3" };
    const char* const stdIn     = "3, 2";
    runTest(CommandLine(cmdline), stdIn);
}

// Maximum distance between group (1, 2) and atom 3
TEST_F(MindistTest, maxDistWorks)
{
    setOutputFile("-od", "mindist.xvg", XvgMatch());
    const char* const cmdline[] = { "mindist", "-max" };
    const char* const stdIn     = "2 3";
    runTest(CommandLine(cmdline), stdIn);
}

/* Particles 1 and 2 are 2 nm away through pbc, but
   should be 3 nm away with no pbc */
TEST_F(MindistTest, noPbcWorks)
{
    setOutputFile("-od", "mindist.xvg", XvgMatch());
    const char* const cmdline[] = { "mindist", "-nopbc" };
    const char* const stdIn     = "0 1";
    runTest(CommandLine(cmdline), stdIn);
}

// Group (1, 2), each res compared to particle 3
TEST_F(MindistTest, resPerTimeWorks)
{
    setOutputFile("-or", "respertime.xvg", XvgMatch());
    const char* const cmdline[] = { "mindist", "-respertime" };
    const char* const stdIn     = "3 2";
    runTest(CommandLine(cmdline), stdIn);
}

TEST_F(MindistTest, matrixWorks)
{
    setOutputFile("-od", "mindist.xvg", XvgMatch());
    const char* const cmdline[] = { "mindist", "-matrix" };
    const char* const stdIn     = "5";
    runTest(CommandLine(cmdline), stdIn);
}

// TODO test periodic image - needs a tpr?

} // namespace
} // namespace test
} // namespace gmx
