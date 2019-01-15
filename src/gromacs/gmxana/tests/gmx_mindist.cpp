/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Tests for gmx mindist.
 *
 * \author Kevin Boyd <kevin.boyd@uconn.edu>
 */

#include "gmxpre.h"

#include <cstdio>
#include <cstdlib>

#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/textreader.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/stdiohelper.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"
#include "testutils/xvgtest.h"

namespace
{

using gmx::test::CommandLine;
using gmx::test::XvgMatch;
using gmx::test::StdioTestHelper;

class MindistTest : public gmx::test::CommandLineTestBase
{
    public:
        MindistTest()
        {
            setInputFile("-f", "mindist_coords.gro");
            setInputFile("-s", "mindist_coords.gro");
            setInputFile("-n", "mindist.ndx");
        }

        void runTest(const CommandLine &args, const char * stringForStdin)
        {
            StdioTestHelper stdioHelper(&fileManager());
            stdioHelper.redirectStringToStdin(stringForStdin);

            CommandLine &cmdline = commandLine();
            cmdline.merge(args);
            ASSERT_EQ(0, gmx_mindist(cmdline.argc(), cmdline.argv()));
            checkOutputFiles();
        }
};
/* mindist_coords.pdb has 3 beads spaced out in a 5 nm box, with the same yz coordinates
 * and x coordinates of 1, 4, and 4.5. Index groups 0, 1, and 2 correspond to the individual beads.
 * Index 3 contains the first 2 beads */

// Mindist between beads 0 and 1 should = 2 (across periodic boundaries)
TEST_F(MindistTest, mindistWorksWithSingleAtoms)
{
    setOutputFile("-od", "mindist.xvg", XvgMatch());
    const char *const  cmdline[] = {
        "mindist"
    };
    const char * const stdIn = "0 1";
    runTest(CommandLine(cmdline), stdIn);
}

// Mindist between group (0, 1) and bead 2 should = 0.5
TEST_F(MindistTest, mindistWorksWithMultipleAtoms)
{
    setOutputFile("-od", "mindist.xvg", XvgMatch());
    // TODO index file reading fails on this test without the -ng flag
    const char *const  cmdline[] = {
        "mindist", "-ng", "1"
    };
    const char * const stdIn = "2 3";
    runTest(CommandLine(cmdline), stdIn);
}

/* Should have 0 contacts with default cutoff */
TEST_F(MindistTest, mindistDoesNotPickUpContacts)
{
    // TODO Default distance changes if used after mindistPicksUpContacts
    setOutputFile("-on", "ncontacts.xvg", XvgMatch());
    const char * const cmdline[] = {
        "mindist", "-ng", "1"
    };
    const char * const stdIn = "0 1";
    runTest(CommandLine(cmdline), stdIn);
}

/* Should pick up one contact with 2.5 nm cutoff */
TEST_F(MindistTest, mindistPicksUpContacts)
{
    setOutputFile("-on", "ncontacts.xvg", XvgMatch());
    const char *const  cmdline[] = {
        "mindist", "-d", "2.5", "-ng", "1"
    };
    const char * const stdIn = "0 1";
    runTest(CommandLine(cmdline), stdIn);
}

TEST_F(MindistTest, matrixWorks)
{

}

TEST_F(MindistTest, maxDistWorks)
{

}

TEST_F(MindistTest, distanceCutoffWorks)
{

}
TEST_F(MindistTest, ngWorks)
{

}

TEST_F(MindistTest, noPbcWorks)
{

}

TEST_F(MindistTest, resPerTimeWorks)
{

}

} //namespace
