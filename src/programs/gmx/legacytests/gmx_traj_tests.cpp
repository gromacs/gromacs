/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Tests for gmx traj
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_integration_tests
 */

#include "testutils/integrationtests.h"
#include "gromacs/options/filenameoption.h"
#include "programs/gmx/gmx.h"

namespace
{

//! Test fixture for gmx traj
typedef gmx::test::IntegrationTestFixture GmxTrajTest;

class GmxTrajTestCanRead : public GmxTrajTest
{
    public:
        GmxTrajTestCanRead()
        {
            createProgramCaller("gmx traj", gmx_cmain);
            redirectStringToStdin("0\n");
            setProgramOption<gmx::FileNameOption>("gmx traj", "s", fileManager.getInputFilePath("spc2.gro"));
            setProgramOption<gmx::FileNameOption>("gmx traj", "ox", fileManager.getTemporaryFilePath("spc2.xvg"));
        }

        int runTest(const char *fileName)
        {
            setProgramOption<gmx::FileNameOption>("gmx traj", "f", fileManager.getInputFilePath(fileName));
            return runProgram("gmx traj");
        }
};

/* There are fancy ways to use "value-parameterized tests" to reduce
 * the duplication below, but for the handful of file formats we care
 * about, the code would probably get longer (and need fancy template
 * stuff, too).
 */

TEST_F(GmxTrajTestCanRead, TRRFile)
{
    runTest("spc2-traj.trr");
}

TEST_F(GmxTrajTestCanRead, XTCFile)
{
    runTest("spc2-traj.xtc");
}

TEST_F(GmxTrajTestCanRead, GROFile)
{
    runTest("spc2-traj.gro");
}

TEST_F(GmxTrajTestCanRead, PDBFile)
{
    runTest("spc2-traj.pdb");
}

TEST_F(GmxTrajTestCanRead, G96File)
{
    runTest("spc2-traj.g96");
}

// .g87 file reading has been broken since at least v4.5
// proposed on gmx-developers for removing that support
TEST_F(GmxTrajTestCanRead, DISABLED_G87File)
{
    redirectStringToStdin("0\n4\n6\n0.0\n");
    runTest("spc2-traj.g87");
}

} // namespace
